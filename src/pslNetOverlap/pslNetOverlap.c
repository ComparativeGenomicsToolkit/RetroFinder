#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "memalloc.h"
#include "psl.h"
#include "obscure.h"
#include "bed.h"
#include "axt.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "nib.h"
#include "fa.h"
#include "dlist.h"
#include "binRange.h"
#include "options.h"
#include "genePred.h"
#include "genePredReader.h"
#include "dystring.h"
#include "retroMrnaInfo.h"
#include "twoBit.h"
#include "hdb.h"
#include "gapCalc.h"
#include "chainNet.h"
#include "chainConnect.h"
#include "chainNetDbLoad.h"
#include "genbank.h"
#include "verbose.h"

#define BEDCOUNT 3
/* max number of nets */
#define MAXNETS 32 
/* max number of levels in a net - gaps count as a level */
#define MAX_NET_LEVELS  15

/* smallest gap in a net that we consider a break in orthology */
#define MINGAP 50
#define MINNETSIZE 150
#define MAXNETSIZE 150000
/* structure for reading net files - read subset of full net*/
struct netSummary {
    int level ;
    int qSize ;
    int qN ;
    int tR ;
    int qFar ;
    char *type ;
    int chainId;
};


/* command line */
static struct optionSpec optionSpecs[] = {
    {"maxNetSize", OPTION_INT},
    {"species1", OPTION_STRING},
    {"species2", OPTION_STRING},
    {"nibdir1", OPTION_STRING},
    {"nibdir2", OPTION_STRING},
    {"mrna1", OPTION_STRING},
    {"mrna2", OPTION_STRING},
    {"bedOut", OPTION_STRING},
    {"score", OPTION_STRING},
    {"snp", OPTION_STRING},
    {"cdsFile", OPTION_STRING},
    {"minDiff", OPTION_INT},
    {"skipBlatMerge", OPTION_BOOLEAN},
    {"computeSS", OPTION_BOOLEAN},
    {"nohead", OPTION_BOOLEAN},
    {"ignoreSize", OPTION_BOOLEAN},
    {"noIntrons", OPTION_BOOLEAN},
    {"minCover", OPTION_FLOAT},
    {"minCoverPseudo", OPTION_FLOAT},
    {"minAli", OPTION_FLOAT},
    {"minAliPseudo", OPTION_FLOAT},
    {"splicedOverlapRatio", OPTION_FLOAT},
    {"minIntronSize", OPTION_INT},
    {"intronSlop", OPTION_INT},
    {"nearTop", OPTION_FLOAT},
    {"minNearTopSize", OPTION_INT},
    {"maxBlockGap", OPTION_INT},
    {"maxRep", OPTION_FLOAT},
    {"maxTrf", OPTION_FLOAT},
    {"stripVersion", OPTION_BOOLEAN},
    {"showall", OPTION_BOOLEAN},
    {NULL, 0}
};

struct hash *rmskHash = NULL, *synHash[MAXNETS] ;
char database[64];
int maxNetSize = MAXNETSIZE;


struct hash *readNetToBinKeeper(char *sizeFileName, char *netFileName)
/* read a truncated net table and return results in hash of binKeeper structure for fast query*/
/* free net in binKeeper to save memory only start/end coord */
/*  select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR from netMm6 */
{
struct binKeeper *bk; 
struct lineFile *lf = lineFileOpen(sizeFileName, TRUE);
struct lineFile *nf = lineFileOpen(netFileName , TRUE);
struct hash *hash = newHash(0);
char *chromRow[2];
char *row[12] ;

while (lineFileRow(lf, chromRow))
    {
    char *name = chromRow[0];
    int size = lineFileNeedNum(lf, chromRow, 1);

    if (hashLookup(hash, name) != NULL)
        warn("Duplicate %s, ignoring all but first\n", name);
    else
        {
        bk = binKeeperNew(0, size);
        assert(size > 1);
	hashAdd(hash, name, bk);
        }
    }
while (lineFileNextRow(nf, row, ArraySize(row)))
    {
    struct netSummary *np = NULL;
    AllocVar(np);
    np->level = sqlUnsigned(row[3]);
    np->qSize = sqlUnsigned(row[6])-sqlUnsigned(row[5]);
    np->qN = sqlUnsigned(row[8]);
    np->tR = sqlUnsigned(row[9]);
    np->qFar = sqlSigned(row[10]);
    np->chainId = sqlSigned(row[11]);
    np->type = cloneString(row[7]);
    char *chrom = cloneString(row[0]);
    int chromStart = sqlUnsigned(row[1]);
    int chromEnd = sqlUnsigned(row[2]);
   
    bk = hashMustFindVal(hash, chrom);
    binKeeperAdd(bk, chromStart, chromEnd, np);
    }
lineFileClose(&nf);
lineFileClose(&lf);
return hash;
}

void binKeeperHashFree(struct hash **hash)
{
if (*hash != NULL)
    {
    struct hashEl *hashEl = NULL;
    struct hashCookie cookie = hashFirst(*hash);
    while ((hashEl = hashNext(&cookie)) != NULL)
        {
        struct binKeeper *bk = hashEl->val;
        binKeeperFree(&bk);
        }
    hashFree(hash);
    }
}
struct hash *readBedCoordToBinKeeper(char *sizeFileName, char *bedFileName, int wordCount)
/* read a list of beds and return results in hash of binKeeper structure for fast query*/
/* free bed in binKeeper to save memory only start/end coord */
{
struct binKeeper *bk; 
struct bed *bed;
struct lineFile *lf = lineFileOpen(sizeFileName, TRUE);
struct lineFile *bf = lineFileOpen(bedFileName , TRUE);
struct hash *hash = newHash(0);
char *chromRow[2];
char *row[3] ;

assert (wordCount == 3);
while (lineFileRow(lf, chromRow))
    {
    char *name = chromRow[0];
    int size = lineFileNeedNum(lf, chromRow, 1);

    if (hashLookup(hash, name) != NULL)
        warn("Duplicate %s, ignoring all but first\n", name);
    else
        {
        bk = binKeeperNew(0, size);
        assert(size > 1);
	hashAdd(hash, name, bk);
        }
    }
while (lineFileNextRow(bf, row, ArraySize(row)))
    {
    bed = bedLoadN(row, wordCount);
    bk = hashMustFindVal(hash, bed->chrom);
    binKeeperAdd(bk, bed->chromStart, bed->chromEnd, bed);
    bedFree(&bed);
    }
lineFileClose(&bf);
lineFileClose(&lf);
return hash;
}

double getChainScore(int chainId, char *seqName, int start, int end, 
        char *database, char *track, char *otherDb, int *overlap)
{
struct chain *chain = NULL, *subChain = NULL, *toFree = NULL;
boolean nullSubset = FALSE;
chain = chainLoadIdRange(database, track, seqName, start, end, chainId);
chainSubsetOnT(chain, start, end, &subChain, &toFree);
if (subChain != NULL)
    {
    *overlap = positiveRangeIntersection(start, end, subChain->tStart, subChain->tEnd);
    verbose(5, "new chain overlap is %d \n", *overlap);
    }
else
    verbose(5, "no chain found %s %d %d %s\n", seqName, start, end, track);
double subSetScore = 0;

if (subChain == NULL)
    nullSubset = TRUE;
else if (hDbIsActive(otherDb) && subChain != chain)
    {
    //char *linearGap = trackDbSettingOrDefault(tdb, "chainLinearGap", "loose");
    struct gapCalc *gapCalc = gapCalcFromFile("medium");
    struct axtScoreScheme *scoreScheme = axtScoreSchemeDefault();
    int qStart = subChain->qStart;
    int qEnd   = subChain->qEnd  ;
    struct dnaSeq *tSeq = hDnaFromSeq(database, subChain->tName, subChain->tStart, subChain->tEnd, dnaLower);
    struct dnaSeq *qSeq = NULL;
    char *matrix = NULL;
    if (matrix != NULL)
        {
        char *words[64];
        int size = chopByWhite(matrix, words, 64) ;
        if (size == 2 && atoi(words[0]) == 16)
            {
            scoreScheme = axtScoreSchemeFromBlastzMatrix(words[1], 400, 30);
            }
        else
            {
            if (size != 2)
                errAbort("error parsing matrix entry in trackDb, expecting 2 word got %d ",
                        size);
            else
                errAbort("error parsing matrix entry in trackDb, size 16 matrix, got %d ",
                        atoi(words[0]));
            }
        }

    if (subChain->qStrand == '-')
        reverseIntRange(&qStart, &qEnd, subChain->qSize);
    qSeq = hChromSeq(otherDb, subChain->qName, qStart, qEnd);
    if (subChain->qStrand == '-')
        reverseComplement(qSeq->dna, qSeq->size);
    subChain->score = chainCalcScoreSubChain(subChain, scoreScheme, gapCalc,
        qSeq, tSeq);
    verbose(5, "%s chainScore %4.0g\n",subChain->qName, subChain->score);
    subSetScore = subChain->score;
    }
else if (subChain == chain)
    subSetScore = chain->score;
chainFree(&toFree);
return subSetScore;
}

int netOverlap(struct psl *psl, struct hash *nHash, char *name)
/* calculate if pseudogene overlaps the syntenic diagonal with another species */
{
//  case 1: [no orthlo DNA in net]
//     net-1                  XXXXXXXXXXXXXXXXXX
//    Repeat                                XXXXXXXXXXX
//                     retro  XXXXXXXXXXXXXXXXXXXXXXXXXX
//                            <----netPart1----> 
//                            <----netSize[1]-->
//                            <--gapsize[2]-------------------------->
//                            <-ovrSize[1]----->  
//  case 2: [no break in net]
//     net-1 XXXXXXXXXXXXXXXXX                                XXXXXXXXXXXXXXXXXX
//     Gap-2                  XXXXXXXXXXXXXXXXXXXXXXXXXNNNNNNN
//     net-3                  XXX(inv)XXXXXXXXXX
//    Repeat                                XXXXXXXXXXX
//                     retro  XXXXXXXXXXXXXX           XXXXXXXXXXXX
//           <----netPart1--->                                <------netPart1-->
//                            <-----------netPart2----------->                  
//                            <----netPart3--->
//           <---------------------------------netSize[1]---------------------->
//                            <--gapsize[2]--->
//                            <-ovrSize[2]-->         <-oSize[2]-->
//                            <----netSize[3]->
//  case 3: [break in net]
//     net-1 XXXXXXXXXXXXXXXXX                                     XXXXXXXXXXXXX
//     Gap-2                  XXXXXXXXXXXXXXXXXXXXXXXXXNNNNNNNNNNN
//     net-3                  XXXX(qFAR>200K)X
//    Repeat                                  XXXXXXXXX
//                     retro  XXXXXXXXXXXXXX           XXXXXXXXXXXX
//           <----netPart1--->                                <------netPart1-->
//                            <-----------netPart2--------------->                  
//                            <---netPart3-->
//           <---------------------------------netSize[1]---------------------->
//                            <-gapsize[2]--->
//                            <-ovrSize[2]-->         <-oSize[2]-->
//                            <---netSize[3]->
//  case 4: [inversion]
int overlap = 0;
int percentBreak = 0;
if (nHash != NULL)
    {
    int maxlevel = 0;
    int maxChainScore = 0;
    int netSize[MAX_NET_LEVELS];
    int gapSize[MAX_NET_LEVELS];
    int overlapSize[MAX_NET_LEVELS];
    double chainScore[MAX_NET_LEVELS];
    char *netType[MAX_NET_LEVELS];
    char netName[128];
    char chainName[128];
    char otherDb[128];
    int rptSize = 0;
    //float coverage = 0.0;
    struct binKeeper *bk = hashFindVal(nHash, psl->tName);
    struct binElement *el, *elist, *rptlist;
    struct binKeeper *rptbk = hashFindVal(rmskHash, psl->tName);
    int i;
    splitPath(name, NULL, netName, NULL);
    safef(otherDb, sizeof(otherDb), "%s",netName);
    strSwapStrs(otherDb, sizeof(otherDb), "net", "");
    otherDb[0] = tolower(otherDb[0]);
    safef(chainName, sizeof(chainName), "%s",netName);
    strSwapStrs(chainName, sizeof(chainName), "net", "chain");
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) netSize[i] = 0;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) gapSize[i] = 0;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) overlapSize[i] = 0;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) netType[i] = NULL;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) chainScore[i] = 0;
    elist = binKeeperFindSorted(bk, psl->tStart , psl->tEnd ) ;
    rptlist = binKeeperFindSorted(rptbk, psl->tStart , psl->tEnd ) ;
    verbose(2, "\n## %s %s %s %s\n",psl->qName, netName, chainName, name);
    for (el = rptlist; el != NULL ; el = el->next)
        {
        rptSize += min(el->end,psl->tEnd)-max(el->start,psl->tStart) ;
        verbose(5, "** subtract repeat %d cum rpt size %d \n", 
                (el->end)-(el->start), rptSize) ;
        }
    int retroSize = psl->tEnd - psl->tStart - rptSize;
    verbose(5, "retrosize %d orig size %d %s\n",retroSize,(psl->tEnd)-(psl->tStart), psl->qName);
    if (retroSize < 0) retroSize = 1;
    /* calculate size gap and net at each level */
    for (el = elist; el != NULL ; el = el->next)
        {
        struct netSummary *np = el->val;
        int netPart = (el->end)-(el->start);
        float overlapThreshold = 0.01;
        int currentLevel = np->level;
        netType[currentLevel] = np->type;
        if (sameString(np->type, "gap"))
            {
            if (((float)netPart/(float)retroSize) > overlapThreshold)
                {
                gapSize[currentLevel] += netPart ;
                verbose(5, "%s gapSize %d netpart %d\n",
                        psl->qName, gapSize[currentLevel], netPart);
                if (np->qN > 0)
                    {
                    int delta = 0;
                    /*treat gaps as syntenic sequence */
                    verbose(5, "%s gap reduced from %d ", 
                            psl->qName, gapSize[currentLevel]); 
                    delta = (float)gapSize[currentLevel]*(float)(np->qN)/(float)np->qSize; 
                    gapSize[currentLevel] -= delta;
                    netPart -= delta;
                    verbose(5, " to %d by netPart*qN %d*%d delta %d\n", 
                            gapSize[currentLevel], netPart, np->qN, delta); 
                    netSize[currentLevel-1] += delta;
                    verbose(5, "%s netSize[%d] increased to %d\n",
                            psl->qName, currentLevel-1, netSize[currentLevel-1]);
                    }
                else
                    {
                    verbose(5, "%s gap reduced from %d to %d by tR %d ", 
                            psl->qName, gapSize[currentLevel], gapSize[currentLevel]-np->tR, np->tR); 
                    gapSize[currentLevel] -= np->tR;
                    }
                if (gapSize[currentLevel] <= 0)
                    gapSize[currentLevel] = 1;
                }
            assert(gapSize[currentLevel]>=0);
            verbose(5, "NET %s gapSize[%d]=%d qN %d tR %d netPart %d\n",
                    psl->qName, np->level, gapSize[currentLevel], np->qN, np->tR, netPart);
            if (currentLevel > 1)
                {
                overlapSize[currentLevel-1] -= netPart;
                if (overlapSize[currentLevel-1] < 0)
                    overlapSize[currentLevel-1] = 0;
                verbose(4, "%s overlap on level %d reduced by netPart %d to %d\n",
                        psl->qName, currentLevel-1, netPart, overlapSize[currentLevel-1]);
                }
            }
        else /* not a gap */
            {
            int overlapPart = positiveRangeIntersection(psl->tStart, psl->tEnd, 
                el->start, el->end);
            overlapSize[np->level] += overlapPart;
            netSize[np->level] += netPart;
            verbose(4,"NET NonGap %s %s:%d netSize[%d]=%d retroSize %d part %d tRep %d qN %d overSize %d %s chain %d\n",
                    psl->qName, psl->tName, psl->tStart, np->level, netSize[np->level], 
                    retroSize, el->end-el->start, np->tR, np->qN, overlapSize[np->level], np->type, np->chainId );
            if (np->level > maxlevel && netSize[np->level] > MINNETSIZE)
                {
                verbose(5,"%s level %d > maxlevel %d , netsize[%d] %d > MINNETSIZE %d\n", 
                        psl->qName, np->level, maxlevel, np->level, netSize[np->level], MINNETSIZE);
                maxlevel = np->level;
                verbose(5, "%s max level = %d\n",psl->qName, maxlevel);
                chainScore[maxlevel] = getChainScore(np->chainId, psl->tName, psl->tStart, psl->tEnd, 
                        database, chainName, otherDb, &overlap);
                if (overlap > overlapSize[maxlevel])
                    overlapSize[maxlevel] = overlap;
                if (chainScore[maxlevel] > maxChainScore)
                    maxChainScore = chainScore[maxlevel];

                verbose(4, "%s max level = %d chainscore %4.0f overlapSize[] %d \n",
                        psl->qName, maxlevel, chainScore[maxlevel], overlapSize[maxlevel]);
                if (sameString(np->type , "inv") && maxlevel > 2 && netSize[maxlevel] < maxNetSize && np->qFar > 0)
                    {
                    overlapSize[maxlevel-2] = overlapSize[maxlevel];
                    verbose(2, "%s INVSERION new level %d, qFar %d overlapSize %d\n",
                            psl->qName, maxlevel, np->qFar, overlapSize[maxlevel]);
                    }
                if (sameString(np->type , "syn") && maxlevel > 2 && 
                        netSize[maxlevel] < maxNetSize && np->qFar > 0  && np->qFar < 10000000)
                    {
                    overlapSize[maxlevel-2] = overlapSize[maxlevel];
                    chainScore[maxlevel-2] = max(chainScore[maxlevel-2], chainScore[maxlevel]+1);
                    verbose(2, "%s Syntenic rearrangement new level %d, qFar %d overlapSize %d cScore %4.0f\n",
                            psl->qName, maxlevel, np->qFar, overlapSize[maxlevel], chainScore[maxlevel-2]);
                    }
//                if (maxlevel > 2)
//                    {
//                    verbose(4, "%s checking gapSize[%d] %d < MINGAP %d\n",
//                            psl->qName, maxlevel-1, gapSize[maxlevel-1], MINGAP);
//                    if (gapSize[maxlevel-1] < MINGAP)
//                        {
//                        maxlevel -= 2;
//                        verbose(4,"%s small gap, NET level bumped down from %d to %d \n",
//                            psl->qName, maxlevel+2, maxlevel);
//                        }
//                    }
                }
            }
        verbose(2, "%s type %s type %s qFar %d level %d gap %d\n",
                psl->qName, np->type, netType[np->level], np->qFar, np->level, gapSize[np->level]);
        }

    for (i = 0 ; i <= maxlevel ; i++) 
        verbose(2, "#### %s %%break %d retroSize %d gapSize[%d] %d overlapSize %d maxlevel %d netSize[%d] %d type %s chainSc %4.0f \n",
            psl->qName, percentBreak, retroSize, i, gapSize[i], 
            overlapSize[i], maxlevel, i , netSize[i], netType[i], chainScore[i]);
    if (maxlevel > 1 && rptSize > 0)
        {
        if (gapSize[maxlevel-1] > netSize[maxlevel])
            {
            gapSize[maxlevel-1] -= rptSize;
            if (gapSize[maxlevel-1] < netSize[maxlevel])
                gapSize[maxlevel-1] = netSize[maxlevel]+1;
            verbose(4, "%s NET gap reduced by %d to %d\n",
                    psl->qName, rptSize, gapSize[maxlevel-1]);
            assert(gapSize[maxlevel-1]);
            }
        }
    int maxChain = 0;
    for (i = maxlevel ; i >= 0 ; i--) 
        /* find largest scoring chain */
        {
        if (chainScore[i] > maxChain)
            {
            maxChain = chainScore[i];
            maxlevel = i;
            }
        }
    if (maxlevel > 2)
        if (sameString(netType[maxlevel] , "inv") && maxlevel > 2 && netSize[maxlevel] < maxNetSize)
            {
            maxlevel-=2;
            verbose(2, "%s INVSERION new level %d, overlapSize %d\n",
                    psl->qName, maxlevel, overlapSize[maxlevel]);
            }

    /* no orthologous DNA, treat as break */
    if (maxlevel == 0)
        {
        percentBreak = 99;
        verbose(4,"%s case1: no ortho DNA maxlevel == %d, percentBreak=%d\n", psl->qName, maxlevel, percentBreak);
        }
    else if (maxlevel==1 )
        {
        if (overlapSize[maxlevel] == 0)
            {
            percentBreak = 99; /* no orthologous DNA, treat as break */
            verbose(4,"%s case1: no ortho DNA maxlevel == %d, percentBreak=%d\n", psl->qName, maxlevel, percentBreak);
            }
        else
            {
            percentBreak = (float)overlapSize[maxlevel]*100/(float)netSize[maxlevel];
            verbose(4,"%s case2: top level maxlevel == %d, percentBreak=%d overlapSize %d\n", 
                psl->qName, maxlevel, percentBreak, overlapSize[maxlevel]);
            }
        }
    /* orthologous DNA, no break */
    else if (netSize[maxlevel] > maxNetSize)
        {
        percentBreak = 0;
        verbose(4,"%s case2: %d net no break maxlevel == %d, percentBreak=%d netSize %d\n", 
                psl->qName, maxNetSize, maxlevel, percentBreak, netSize[maxlevel]);
        }
    /* sequence gaps should not be treated as breaks in orthology */
//    else if (2*gapSize[maxlevel-1] < netSize[maxlevel])
//        {
//        verbose(4, "NET gap %d < net %d zero break\n",gapSize[maxlevel-1], netSize[maxlevel]);
//        percentBreak = 0;
//        }
    else if (gapSize[maxlevel-1] < maxNetSize)
        {
        percentBreak = (float)retroSize*100/(float)gapSize[maxlevel-1];
        verbose(4, "%s case3: gap<%d retro %d / gap %d  netSize[%d]=%d overlapSize %d break %d \n",
                psl->qName, maxNetSize, retroSize, gapSize[maxlevel-1], maxlevel, 
                netSize[maxlevel], overlapSize[maxlevel], percentBreak );
        if (percentBreak >= 130) percentBreak = 120;
        }
    else if (gapSize[maxlevel-1] >= maxNetSize)
        {
        if (overlapSize[maxlevel] == 0)
            percentBreak = 99;
        else
            percentBreak = (float)overlapSize[maxlevel]*100/(float)retroSize;
        if (percentBreak > 100)
            percentBreak = 101;
        verbose(4, "%s case3: gap>%d retro %d / gap %d  break %d netSize[%d]=%d overlapSize %d\n",
                psl->qName, maxNetSize, retroSize, gapSize[maxlevel-1], percentBreak, maxlevel, 
                netSize[maxlevel], overlapSize[maxlevel]);
        }
    //verbose(3,"NET AFTER #score %s %s:%d maxlevel %d  overlap %d percentBreak=%d netSize[max] = %d \n",
            //psl->qName, psl->tName, psl->tStart, maxlevel, overlapSize[maxlevel] , 
            //percentBreak, netSize[maxlevel]);
    //slFreeList(&elist);
    }
assert (percentBreak < 130);
return percentBreak;
}
void usage(int count)
/* Print usage instructions and exit. */
{
errAbort(
    "pslNetOverlap - calculate break in orthlogy based on nets\n"
    "and mRNA alignments \n"
    "usage:\n"
    "    pslPseudo database input.psl chrom.sizes rmsk.bed net.txt(s) ...\n\n"
    "where \n"
    "\tinput.psl - alignment of mRNAs sorted by pslSort, qName should be unique\n"
    "\tchrom.sizes - list of chromosome followed by size\n"
    "\trmsk.bed - bed file with coordinates of repeats\n" 
    "\tnet.txt - subset of net tab delimited with 12 columns, select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf, qFar, chainId from net\n"
    "\noptions:\n"
    "\t-verbose=N - higher N results in more debugging output\n"
    "\t-maxNetSize=N - largeest net that can still be considered a break in orthlogy default = %d\n\n, ",MAXNETSIZE
    );
}
int main(int argc, char *argv[])
/* Process command line. */
{
struct psl *psl = NULL, *pslList = NULL;
int lineSize;
char *line;
char *words[32];
int wordCount;
int aliCount = 0;
int i;

optionInit(&argc, argv, optionSpecs);
if (argc < 6)
    usage(argc);
//verboseSetLogFile("stdout");
//verbosity = optionInt("verbose", 1);
//verboseSetLevel(verbosity);
//verbose(1,"version is %s\n",rcsid);
maxNetSize = optionInt("maxNetSize", MAXNETSIZE);
safef(database, sizeof(database), "%s",argv[1]);
struct lineFile *in = pslFileOpen(argv[2]);
verbose(2,"Reading Repeats from %s\n",argv[4]);
rmskHash = readBedCoordToBinKeeper(argv[3], argv[4], BEDCOUNT);
printf("#qName\t");
for (i = 5 ; i < argc ; i++)
    {
    verbose(2,"Loading net %s\n",argv[i]);
    printf("%s\t",argv[i]);
    synHash[i] = readNetToBinKeeper(argv[3], argv[i]);
    }
printf("\n");
while (lineFileNext(in, &line, &lineSize))
    {
    if ((++aliCount & 0x1ffff) == 0)
        {
	verboseDot();
	}
    wordCount = chopTabs(line, words);
    if (wordCount != 21)
	errAbort("Bad line %d of %s\n", in->lineIx, in->fileName);
    psl = pslLoad(words);
    printf("%s\t",psl->qName);
    for (i = 5 ; i < argc ; i++)
        {
        int overlapPercent = netOverlap(psl, synHash[i], argv[i]);
        printf("%d\t", overlapPercent);
        }
    printf("\n");
    slAddHead(&pslList, psl);
    }
slReverse(&pslList);


//binKeeperHashFree(&synHash);
return 0;
}
