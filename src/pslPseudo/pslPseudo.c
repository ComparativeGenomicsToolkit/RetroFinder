/* pslPseudo - analyse repeats and generate list of processed pseudogenes
 * from a genome wide, sorted by mRNA .psl alignment file.
 * 
 * If you are looking at this for the first time, start by reading the initWeights function. It describes the features and weights that make up
 * the score function.
 */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "memalloc.h"
#include "jksql.h"
#include "psl.h"
#include "obscure.h"
#include "bed.h"
#include "axt.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "dlist.h"
#include "binRange.h"
#include "options.h"
#include "genePred.h"
#include "genePredReader.h"
#include "dystring.h"
#include "ucscRetroInfo.h"
#include "ucscRetroOrtho.h"
#include "twoBit.h"
#include "chainNet.h"
#include "genbank.h"
#include "verbose.h"

#define ScoreNorm 7
#define BEDCOUNT 3
#define POLYASLIDINGWINDOW 10
#define POLYAREGION 150
#define INTRONMAGIC 10 /* allow more introns if lots of exons covered - if (exonCover - intronCount > INTRONMAGIC) */

/* label for classification stored in ucscRetroInfo table */
#define PSEUDO 1
#define NOTPSEUDO -1
#define EXPRESSED -2

/* max number of levels in a net - gaps count as a level */
#define MAX_NET_LEVELS  15

/* smallest gap in a net that we consider a break in orthology */
#define MINGAP 50
#define MINNETSIZE 150

static char const rcsid[] = "$Id: pslPseudo.c,v 1.18 2012/09/22 15:17:05 baertsch Exp $";

char *db;
char *mrnaSeq;
int verbosity = 1;
float minAli = 0.98;
float maxRep = 0.60;
float maxTrf = 0.50;
float minAliPseudo = 0.60;
float nearTop = 0.005;
float repsPerIntron = 0.7;
float splicedOverlapRatio = 0.1;
float minCover = 0.50;
float minCoverPseudo = 0.01;
int maxBlockGap = 50;
int intronSlop = 35;
int minIntronSize = 80;
int spliceDrift = 15;
int scoreThreshold = 425;
float intronRatio = 1.5;
bool skipBlatMerge = FALSE;
bool ignoreSize = FALSE;
bool noIntrons = FALSE;
bool noHead = FALSE;
bool skipExciseRepeats = FALSE;
bool stripVersion = FALSE;
bool showall = FALSE;
double  wt[12];     /* weights on score function*/
int minNearTopSize = 10;
struct genePred *gpList1 = NULL, *gpList2 = NULL, *kgList = NULL, *mrnaGene = NULL;
FILE *bestFile, *pseudoFile, *linkFile, *axtFile, *orthoFile;
struct twoBitFile *genomeSeqFile = NULL;
struct twoBitFile *mrnaFile = NULL;
struct axtScoreScheme *ss = NULL; /* blastz scoring matrix */
struct slName *mrnaList = NULL;/* list of all input mrna sequence names  */
struct hash *fileHash = NULL;  
char mrnaOverlap[255];
struct hash *rmskHash = NULL, *trfHash = NULL, *synHash = NULL, *syn2Hash, *syn3Hash, *exprHash = NULL, *qNameHash = NULL;
/* hash table of accession to CDS */
static struct hash *gCdsTable = NULL;
bool abortAtEnd = FALSE;
char orthoNet1[128];
char orthoNet2[128];
char orthoNet3[128];

/* structure for reading net files - read subset of full net*/
struct netSummary {
    int level ;
    int qSize ;
    int qN ;
    int tR ;
    char *type ;
};

/* command line */
static struct optionSpec optionSpecs[] = {
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

void usage(int count)
/* Print usage instructions and exit. */
{
errAbort(
    "pslPseudo - analyse repeats and generate genome wide best\n"
    "alignments from a sorted set of local alignments\n"
    "Usage:\n"
    "    pslPseudo db in.psl sizes.lst rmsk.bed mouseNet.txt dogNet.txt trf.bed all_mrna.psl out.psl outPseudo.psl outPseudoLink.txt out.axt genome.2bit mrna.2bit refGene.tab mgcGene.tab kglist.tab rheMac2Net.txt orthoFileName \n\n"
    "where in.psl is an blat alignment of mrnas sorted by pslSort\n"
    "blastz.psl is an blastz alignment of mrnas sorted by pslSort\n"
    "sizes.lst is a list of chromosome followed by size\n"
    "rmsk.bed is the repeat masker bed file\n"
    "mouseNet.txt = select tName, tStart, tEnd, level, qName, qStart, qEnd, type from netMm6\n"
    "dogNet.txt = select tName, tStart, tEnd, level, qName, qStart, qEnd, type from netCanFam2\n"
    "trf.bed is the simple repeat (trf) bed file\n"
    "all_mrna.psl is the blat best mrna alignments\n"
    "out.psl is the best mrna alignment for the gene \n"
    "outPseudo.psl contains pseudogenes\n"
    "outPseudoLink.txt will have the link between gene and pseudogene\n"
    "out.axt contains the pseudogene aligned to the gene (can be /dev/null)\n"
    "genome.2bit genome sequence in 2bit format\n"
    "mrna.2bit  sequence data for all aligned mrnas using lastz\n"
    "refGene.tab refseq annotations for finding parent gene (genePred format)\n"
    "mgcGene.tab mgc annotations for finding parent gene (genePred format)\n"
    "kglist.tab knownGene annotations for finding parent gene (genePred format)\n"
    "mrnaGene.tab mrna cds  annotations (genePred format)\n"
    "rhesusNet.txt = select tName, tStart, tEnd, level, qName, qStart, qEnd, type from netRheMac2\n"
    "orthoFile.txt = \n"
    "options:\n"
    "    -nohead don't add PSL header\n"
    "    -skipBlatMerge Do not merge alignments from all_mrna into input in.psl. Avoids size mismatches. \n"
    "    -ignoreSize Will not weigh in favor of larger alignments so much\n"
    "    -noIntrons Will not penalize for not having introns when calculating size factor\n"
    "    -minCover=0.N minimum coverage for parent mrna , default %4.3f\n"
    "    -minCoverPseudo=0.N minimum coverage of pseudogene to output, default %4.3f\n"
    "    -minAli=0.N minimum alignment ratio for mrna, default %4.3f\n"
    "    -minAliPseudo=0.N minimum alignment ratio for pseudogenes, default %4.3f\n"
    "    -minIntronSize=N minium size alignment gap to be counted as an intron, default %d \n"
    "    -splicedOverlapRatio=0.N max overlap with spliced mrna,  default %4.3f\n"
    "    -intronSlop=N max delta of intron position on q side alignment, default %d bp\n"
    "    -nearTop=0.N how much can deviate from top and be taken, default %4.3f\n"
    "    -minNearTopSize=N  Minimum size of alignment that is near top for aligmnent to be kept,  default %d.\n"
    "    -maxBlockGap=N  Max gap size between adjacent blocks that are combined, default %d.\n"
    "    -maxRep=N  max ratio of overlap with repeat masker\n"
    "               for aligmnent to be kept,  default %4.3f\n"
    "    -maxTrf=N  max ratio of overlap with simple repeats\n"
    "               for aligmnent to be kept,  default %4.3f\n"
    "    -stripVersion  ignore version number of mRNA in input file \n"
    "    -cdsFile=file - get CDS from this database, psl is a file.\n"
    "    -showall  do not eliminate low scoring hits from output file , default=false\n"
    " Arg count is %d , should be 20\n"
       , minCover, minCoverPseudo, minAli, minAliPseudo, minIntronSize, splicedOverlapRatio, 
       intronSlop, nearTop, minNearTopSize, maxBlockGap, maxRep, maxTrf, count );
}

bool samePrefix(char *string1, char *string2, char sep)
/* check if accessions match ignoring suffix and version */
{
char buf1[256], buf2[256];
char *name1[3];
char *name2[3];
safef(buf1, sizeof(buf1), "%s",string1);
safef(buf2, sizeof(buf2), "%s",string2);
chopString(buf1, &sep, name1, ArraySize(name1));
chopString(buf2, &sep, name2, ArraySize(name2));
if (name1[0] == NULL)
    name1[0] = string1;
if (name2[0] == NULL)
    name2[0] = string2;
if (sameString(name1[0], name2[0]))
    {
    return TRUE;
    }
return FALSE;
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
char *row[10] ;

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

struct axt *axtCreate(char *q, char *t, int size, struct psl *psl)
/* create axt */
{
int qs = psl->qStart, qe = psl->qEnd;
int ts = psl->tStart, te = psl->tEnd;
int symCount = 0;
struct axt *axt = NULL;

AllocVar(axt);
if (psl->strand[0] == '-')
    reverseIntRange(&qs, &qe, psl->qSize);

if (psl->strand[1] == '-')
    reverseIntRange(&ts, &te, psl->tSize);

axt->qName = cloneString(psl->qName);
axt->tName = cloneString(psl->tName);
axt->qStart = qs;
axt->qEnd = qe;
axt->qStrand = psl->strand[0];
axt->tStrand = '+';
if (psl->strand[1] != 0)
    {
    axt->tStart = ts;
    axt->tEnd = te;
    }
else
    {
    axt->tStart = psl->tStart;
    axt->tEnd = psl->tEnd;
    }
axt->symCount = symCount = strlen(t);
axt->tSym = cloneString(t);
if (strlen(q) != symCount)
    warn("Symbol count %d != %d inconsistent at t %s:%d and qName %s\n%s\n%s\n",
    	symCount, (int)strlen(q), psl->tName, psl->tStart+1, psl->qName, t, q);
axt->qSym = cloneString(q);
axt->score = axtScoreFilterRepeats(axt, ss);
verbose(1,"axt score = %d\n",axt->score);
return axt;
}

void loadCdsFile(char *cdsFile)
/* read a CDS file into the global hash */
{
struct lineFile *lf = lineFileOpen(cdsFile, TRUE);
char *row[4];

gCdsTable = hashNew(20);
while (lineFileNextRowTab(lf, row, 4))
    hashAdd(gCdsTable, row[0], 
            lmCloneString(gCdsTable->lm, row[2]));

lineFileClose(&lf);
}

char *cdsQuery(char *acc)
/* query for a CDS, either in the hash table */
{
if (gCdsTable != NULL)
    return hashFindVal(gCdsTable, acc);
else
    return NULL;
}

int getGbVersionIdx(char *acc)
/* determine if a accession appears to have a genbank version followed by optional dash delimited extension.  If
 * so, return index of the dot, otherwise -1 */
{
char *verPtr = strrchr(acc, '.');
int dotIdx;
if (verPtr == NULL)
    return -1;
dotIdx = verPtr - acc;
verPtr++;
while (*verPtr != '\0')
    {
    if (!isdigit(*verPtr) && *verPtr != '-')
        return -1;
    verPtr++;
    }
return dotIdx;
}

char *getCdsForAcc(char *acc)
/* look up a cds, trying with and without version */
{
/* get the CDS for acc */
char *cdsStr = cdsQuery(acc);
/* if it is not found then try using the accession without version number */
if (cdsStr == NULL)
    {
    int dotIdx = getGbVersionIdx(acc);
    if (dotIdx >= 0)
        {
        acc[dotIdx] = '\0';
        cdsStr = cdsQuery(acc);
        acc[dotIdx] = '.';
        }
    }
if (cdsStr != NULL)
    verbose(5, "cds str %s %s\n",acc, cdsStr);
else
    verbose(5, "cds str %s NULL\n",acc);
return cdsStr;
}

void getCdsInMrna(struct genePred *gp, int *retCdsStart, int *retCdsEnd)
/* Given a gene prediction, figure out the
 * CDS start and end in mRNA coordinates. */
{
int missingStart = 0, missingEnd = 0;
int exonStart, exonEnd, exonSize, exonIx;
int totalSize = 0;

for (exonIx = 0; exonIx < gp->exonCount; ++exonIx)
    {
    exonStart = gp->exonStarts[exonIx];
    exonEnd = gp->exonEnds[exonIx];
    exonSize = exonEnd - exonStart;
    totalSize += exonSize;
    missingStart += exonSize - positiveRangeIntersection(exonStart, exonEnd, gp->cdsStart, exonEnd);
    missingEnd += exonSize - positiveRangeIntersection(exonStart, exonEnd, exonStart, gp->cdsEnd);
    }
*retCdsStart = missingStart;
*retCdsEnd = totalSize - missingEnd;
}

void writeInsert(struct dyString *aRes, struct dyString *bRes, char *aSeq, int gapSize)
/* Write out gap, possibly shortened, to aRes, bRes. */
{
dyStringAppendN(aRes, aSeq, gapSize);
dyStringAppendMultiC(bRes, '-', gapSize);
}

void writeGap(struct dyString *aRes, int aGap, char *aSeq, struct dyString *bRes, int bGap, char *bSeq)
/* Write double - gap.  Something like:
 *         --c
 *         ag-  */
{
/* append bGap '-' characters to end of aRes */
dyStringAppendMultiC(aRes, '-', bGap);
/* bSeq is appended to bRes */
dyStringAppendN(bRes, bSeq, bGap);
/* aSeq is appended to aRes */
dyStringAppendN(aRes, aSeq, aGap);
/* append aGap '-' characters to end of bRes */
dyStringAppendMultiC(bRes, '-', aGap);
}

struct genePred *getOverlappingGene2(struct genePred **list, char *table, char *chrom, int cStart, int cEnd, char *name, int *retOverlap)
{
/* read all genes from a table find the gene with the biggest overlap. 
   Cache the list of genes to so we only read it once */

struct genePred *el = NULL, *bestMatch = NULL, *gp = NULL;
int overlap = 0 , bestOverlap = 0, i;
int *eFrames;

for (el = *list; el != NULL; el = el->next)
    {
    if (chrom != NULL && el->chrom != NULL)
        {
        overlap = 0;
        if ( sameString(chrom, el->chrom))
            {
            for (i = 0 ; i<(el->exonCount); i++)
                {
                overlap += positiveRangeIntersection(cStart,cEnd, el->exonStarts[i], el->exonEnds[i]) ;
                }
            if (overlap > 20 && sameString(name, el->name))
                {
                bestMatch = el;
                bestOverlap = overlap;
                *retOverlap = bestOverlap;
                }
            if (overlap > bestOverlap)
                {
                bestMatch = el;
                bestOverlap = overlap;
                *retOverlap = bestOverlap;
                }
            }
        }
    }
if (bestMatch != NULL)
    {
    /* Allocate genePred and fill in values. */
    AllocVar(gp);
    gp->name = cloneString(bestMatch->name);
    gp->chrom = cloneString(bestMatch->chrom);
    gp->strand[1] = bestMatch->strand[1];
    gp->strand[0] = bestMatch->strand[0];
    gp->txStart = bestMatch->txStart;
    gp->txEnd = bestMatch->txEnd;
    gp->cdsStart = bestMatch->cdsStart;
    gp->cdsEnd = bestMatch->cdsEnd;
    gp->exonCount = bestMatch->exonCount;
    AllocArray(gp->exonStarts, bestMatch->exonCount);
    AllocArray(gp->exonEnds, bestMatch->exonCount);
    for (i=0; i<bestMatch->exonCount; ++i)
        {
        gp->exonStarts[i] = bestMatch->exonStarts[i] ;
        gp->exonEnds[i] = bestMatch->exonEnds[i] ;
        }
    gp->optFields = bestMatch->optFields;
    gp->score = bestMatch->score;

    if (bestMatch->optFields & genePredName2Fld)
        gp->name2 = cloneString(bestMatch->name2);
    else
        gp->name2 = NULL;
    if (bestMatch->optFields & genePredCdsStatFld)
        {
        gp->cdsStartStat = bestMatch->cdsStartStat;
        gp->cdsEndStat = bestMatch->cdsEndStat;
        }
    if (bestMatch->optFields & genePredExonFramesFld)
        {
        gp->exonFrames = AllocArray(eFrames, bestMatch->exonCount);
        for (i = 0; i < bestMatch->exonCount; i++)
            gp->exonFrames[i] = bestMatch->exonFrames[i];
        }
    eFrames = gp->exonFrames;
    }

return gp;
}

struct axt *pslToAxt(struct psl *psl, struct genbankCds *cds)
{
static struct dnaSeq *tSeq = NULL, *qSeq = NULL;
static struct slName *mrna;
struct dyString *q = newDyString(16*1024);
struct dyString *t = newDyString(16*1024);
int blockIx;
int qs, ts ;
int lastQ = 0, lastT = 0, size;
int qOffset = 0;
int tOffset = psl->tStart;
struct axt *axt = NULL;
char name[512];

freeDnaSeq(&qSeq);	
assert(mrnaList != NULL);
for (mrna = mrnaList; mrna != NULL ; mrna = mrna->next)
    {
    safef(name, sizeof(name), "%s", psl->qName);
    chopSuffixAt(name, '-');
    assert(mrna != NULL);
    if (sameString(mrna->name, name))
        {
        qSeq = twoBitReadSeqFrag(mrnaFile, name, 0, 0);
        toLowerN(qSeq->dna, qSeq->size);
        if (cds != NULL)
            {
            int i;
//            assert(cds->end <= qSeq->size);
            for (i = cds->start ; i < min(cds->end,qSeq->size) ; i++)
                {
                qSeq->dna[i] = toupper(qSeq->dna[i]);
                }
            verbose(5, "AFTER cds %d-%d sz %d %s\n",cds->start, cds->end, qSeq->size, qSeq->dna);
            }
        assert(qSeq != NULL);
        if(abs((qSeq->size)-psl->qSize) >= 3)
        {
            warn("Error: psl %s qSize = %d and sequence len is %d\n",
                    name, psl->qSize, qSeq->size);
            verbose(1,"Error: psl %s qSize = %d and sequence len is %d\n",
                    name, psl->qSize, qSeq->size);
        }
        assert(abs((qSeq->size)-psl->qSize) < 3); 
        break;
        }
    }
if (qSeq == NULL)
    {
    verbose(5,"mrna sequence data not found %s %s:%d-%d\n",psl->qName, psl->tName, psl->tStart+1,psl->tEnd);
    dyStringFree(&q);
    dyStringFree(&t);
    dnaSeqFree(&tSeq);
    dnaSeqFree(&qSeq);
    return NULL;
    }
if (qSeq->size != psl->qSize)
    {
    warn("sequence %s aligned is different size %d from mrna.fa file %d \n",name,psl->qSize,qSeq->size);
    verbose(2,"sequence %s aligned is different size %d from mrna.fa file %d \n",name,psl->qSize,qSeq->size);
    dyStringFree(&q);
    dyStringFree(&t);
    dnaSeqFree(&tSeq);
    dnaSeqFree(&qSeq);
    return NULL;
    }
freeDnaSeq(&tSeq);
// Read region of genome sequence.
tSeq = twoBitReadSeqFrag(genomeSeqFile, psl->tName, psl->tStart, psl->tEnd);

if (psl->strand[0] == '-')
   {
    reverseComplement(qSeq->dna, qSeq->size);
   }
if (psl->strand[1] == '-')
   {
    reverseComplement(tSeq->dna, tSeq->size);
   }
for (blockIx=0; blockIx < psl->blockCount; ++blockIx)
    {
    qs = psl->qStarts[blockIx] - qOffset;
    ts = psl->tStarts[blockIx] - tOffset;

    if (blockIx != 0)
        {
	int qGap, tGap, minGap;
	qGap = qs - lastQ;
	tGap = ts - lastT;
	minGap = min(qGap, tGap);
	if (minGap > 0)
	    {
	    writeGap(q, qGap, qSeq->dna + lastQ, t, tGap, tSeq->dna + lastT);
	    }
	else if (qGap > 0)
	    {
	    writeInsert(q, t, qSeq->dna + lastQ, qGap);
	    }
	else if (tGap > 0)
	    {
	    writeInsert(t, q, tSeq->dna + lastT, tGap);
	    }
	}
    size = psl->blockSizes[blockIx];
    assert(qSeq != NULL);
    dyStringAppendN(q, qSeq->dna + qs, size);
    lastQ = qs + size;
    dyStringAppendN(t, tSeq->dna + ts, size);
    lastT = ts + size;
    }

if (strlen(q->string) != strlen(t->string))
    warn("Symbol count(t) %d != %d inconsistent at t %s:%d and qName %s\n%s\n%s\n",
    	(int)strlen(t->string), (int)strlen(q->string), psl->tName, psl->tStart+1, name, t->string, q->string);
//if (psl->strand[0] == '-')
//    {
//    reverseComplement(q->string, q->stringSize);
//    reverseComplement(t->string, t->stringSize);
//    }
axt = axtCreate(q->string, t->string, min(q->stringSize,t->stringSize), psl);
dyStringFree(&q);
dyStringFree(&t);
dnaSeqFree(&tSeq);
dnaSeqFree(&qSeq);
return axt;
}

void binKeeperPslHashFree(struct hash **hash)
{
if (*hash != NULL)
    {
    struct hashEl *hashEl = NULL;
    struct hashCookie cookie = hashFirst(*hash);
    while ((hashEl = hashNext(&cookie)) != NULL)
        {
        struct binKeeper *bk = hashEl->val;
        struct binElement *elist = NULL, *el = NULL;;
        elist = binKeeperFindAll(bk) ;
        for (el = elist; el != NULL ; el = el->next)
            {
            struct psl *psl = el->val;
            pslFreeList(&psl);
            }
        binKeeperFree(&bk);
        }
    hashFree(hash);
    }
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

bool checkForGene(struct genePred **list, char *table, char *chrom, int cStart, int cEnd, char *name, int *retOverlap, struct genePred **gp)
{
/* check if this coor/name matches something in gene list
   Cache the list of genes to so we only read it once */

struct genePred *el = NULL, *bestMatch = NULL;
int overlap = 0 , bestOverlap = 0, i;
int cdsStart = -1;
int cdsEnd = -1;

if (*list == NULL)
    return FALSE;
/* go through the gene list and count overlapping bases between exons of the 
 * genePreds in the list and cStart to cEnd range if on the same chrom 
 * and if coordinates and name have a match then return TRUE and the genePred
 * and overlap */
for (el = *list; el != NULL; el = el->next)
    {
    if (chrom != NULL && el->chrom != NULL)
        {
        overlap = 0;
        if ( sameString(chrom, el->chrom))
            {
            for (i = 0 ; i<(el->exonCount); i++)
                {
                overlap += positiveRangeIntersection(cStart,cEnd, el->exonStarts[i], el->exonEnds[i]) ;
                }
            if (overlap > 20 && sameString(name, el->name))
                {
                bestMatch = el;
                bestOverlap = overlap;
                cdsStart = el->cdsStart;
                cdsEnd = el->cdsEnd;
                *retOverlap = bestOverlap;
                if (gp != NULL)
                    *gp = bestMatch;
                verbose(5, "genePred found %s\n",bestMatch->name);
                return TRUE;
                }
            if (overlap > bestOverlap)
                {
                bestMatch = el;
                cdsStart = el->cdsStart;
                cdsEnd = el->cdsEnd;
                bestOverlap = overlap;
                *retOverlap = bestOverlap;
                }
            }
        }
    }

verbose(5, "genePred not found %s\n",name);
return FALSE;
}

int getGenePred(struct ucscRetroInfo *pg, struct psl *psl, struct genePred **retGp)
/* return 0 if no gene pred or positive number to bump the score to 
 * favor parents with CDS annotation */
{
int geneOverlap = 0;
char name[512];
if (psl == NULL) 
    return 0;
safef(name, sizeof(name), "%s", psl->qName);
chopSuffix(name);
if (checkForGene(&gpList2, "mgcGene", psl->tName, psl->tStart, 
                        psl->tEnd , name, &geneOverlap, retGp))
    return 100;
if (checkForGene(&gpList1, "refGene", psl->tName, psl->tStart, 
                        psl->tEnd , name, &geneOverlap, retGp))
    return 70;
if (checkForGene(&kgList, "knownGene", psl->tName, psl->tStart, 
                        psl->tEnd , name, &geneOverlap, retGp))
    return 50;
return 0;
}

int calcMilliScore(struct psl *psl)
/* Figure out percentage score. */
/* Is this really a percentage? for BLAT, use 100.0 - milliBad * 0.1 */
{
return 1000-pslCalcMilliBad(psl, TRUE);
}

void outputNoLinkScore(struct psl *psl, struct ucscRetroInfo *pg, 
                int overlapOrtho1, int overlapOrtho2, int overlapOrtho3)
/* output bed record with pseudogene details */
{
struct ucscRetroOrtho *puro1 = NULL;
struct ucscRetroOrtho *puro2 = NULL;
struct ucscRetroOrtho *puro3 = NULL;
int blockCount;
int i, *chromStarts, chromStart;
AllocVar(puro1);
AllocVar(puro2);
AllocVar(puro3);

pg->chrom = cloneString(psl->tName);
pg->chromStart = pg->thickStart = chromStart = psl->tStart;
pg->chromEnd = pg->thickEnd = psl->tEnd;
strncpy(pg->strand,  psl->strand, sizeof(pg->strand));
pg->blockCount = blockCount = psl->blockCount;
pg->blockSizes = (int *)cloneMem(psl->blockSizes,(sizeof(int)*psl->blockCount));
pg->chromStarts = chromStarts = (int *)cloneMem(psl->tStarts, (sizeof(int)*psl->blockCount));

/* Switch minus target strand to plus strand. */
if (psl->strand[1] == '-')
    {
    int chromSize = psl->tSize;
    reverseInts(pg->blockSizes, blockCount);
    reverseInts(chromStarts, blockCount);
    for (i=0; i<blockCount; ++i)
	chromStarts[i] = chromSize - chromStarts[i];
    }

/* Convert coordinates to relative. */
for (i=0; i<blockCount; ++i)
    chromStarts[i] -= chromStart;

if (pg->type == NULL)
    pg->type = cloneString("NONE");
if (strlen(pg->type)<=1)
    pg->type = cloneString("NONE");
if (pg->gChrom==NULL)
    pg->gChrom = cloneString("NONE");
pg->polyAstart = (psl->strand[0] == '+') ? pg->polyAstart - psl->tEnd : psl->tStart - pg->polyAstart ;
pg->matches = psl->match+psl->repMatch ;
pg->qSize = psl->qSize;
pg->qEnd = psl->qEnd;
if (pg->overName == NULL)
    pg->overName = cloneString("none");
if (strlen(pg->overStrand) < 1)
    {
    pg->overStrand[0] = '0'; 
    pg->overStrand[1] = '\0'; 
    }
verbose(3, "## %s score = %d\n",pg->name, pg->score);
ucscRetroInfoOutput(pg, linkFile, '\t','\n');
puro1->name = pg->name;
puro2->name = pg->name;
puro3->name = pg->name;
puro1->db = orthoNet1;
puro2->db = orthoNet2;
puro3->db = orthoNet3;
puro1->overlap = overlapOrtho1;
puro2->overlap = overlapOrtho2;
puro3->overlap = overlapOrtho3;
ucscRetroOrthoOutput(puro1, orthoFile, '\t','\n');
ucscRetroOrthoOutput(puro2, orthoFile, '\t','\n');
ucscRetroOrthoOutput(puro3, orthoFile, '\t','\n');
//ucscRetroOrthoFree(&puro1);
//ucscRetroOrthoFree(&puro2);
//ucscRetroOrthoFree(&puro3);
}

void initWeights()
{
/*
    0 = + milliBad - measures the divergence from between the parent gene and the retro (1000=perfect match, 600 is very low homology)
    1 = + exon Coverage - counts how many exons from the parent gene are "covered" by the alignment to the retrogene.
    2 = + log axtScore - lastz score between the parent and the retrogene
    3 = + log polyAlen - length of the polyA tail inserted into the retro location
    4 = + max(overlapMouse overlapDog) - max(overlap with break in synteny with other two species )
    5 = + processedIntrons - count of the number of introns from the parent that were removed in the retro location
    6 = + intronCount ^.5
    7 = - log maxOverlap - downweight retros that have evidence of expression - 
		 percentage of bases from retro that overlap all_mnra.psl (arg #8) 
    8 = + coverage *((qSize-qEnd)/qSize)
    9 = - repeats
    10 =- alignGapCount

pseudoScore = ( wt[0]*scaledMilliBad
                + wt[1]*(log(pg->exonCover+1)/log(2))*600       
                + wt[2]*(((log(pg->axtScore>0?pg->axtScore:1)/log(2))*170)-1000)   
                + wt[3]*(log(pg->polyAlen+2)*200)                               
                + wt[4]*(overlapOrtholog*10)                            
                + wt[5]*(((log(pg->processedIntrons > 0 ? (pg->processedIntrons)+1 : 1))/log(2))*1200 ) 
                - wt[6]*pow(pg->intronCount,0.5)*1200                   
                - wt[7]*(maxOverlap*300)                       
                + wt[8]*scaledCoverage 
                - wt[9]*(pg->tReps*10)
                - wt[10]*(pg->alignGapCount) 
                ) / ScoreNorm + bump;

 */
wt[0] = 1; wt[1] = 1; wt[2] = 0.1; wt[3] = 0.2; wt[4] = 0.7; 
wt[5] = 1; wt[6] = 1; wt[7] = 0.5; wt[8] = 1; wt[9] = 1; wt[10] = 1;

}

int netOverlap(struct psl *psl, struct hash *nHash)
/* calculate if pseudogene overlaps the syntenic diagonal with another species */
{
int percentBreak = 0;
if (nHash != NULL)
    {
    int maxlevel = 0;
    int netSize[MAX_NET_LEVELS];
    int gapSize[MAX_NET_LEVELS];
    int overlapSize[MAX_NET_LEVELS];
    int rptSize = 0;
    struct binKeeper *bk = hashFindVal(nHash, psl->tName);
    struct binElement *el;
    struct binElement *rptlist = NULL;
    struct binElement *elist = NULL;
    struct binKeeper *rptbk = hashFindVal(rmskHash, psl->tName);
    int i;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) netSize[i] = 0;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) gapSize[i] = 0;
    for (i = 0 ; i < MAX_NET_LEVELS ; i++) overlapSize[i] = 0;
    if (bk != NULL)
        elist = binKeeperFindSorted(bk, psl->tStart , psl->tEnd ) ;
    if (rptbk != NULL)
        rptlist = binKeeperFindSorted(rptbk, psl->tStart , psl->tEnd ) ;
    for (el = rptlist; el != NULL ; el = el->next)
        {
        rptSize += min(el->end,psl->tEnd)-max(el->start,psl->tStart) ;
        verbose(5, "NET add rpt size %d to cum rpt %d retro %d-%d\n", 
                (el->end)-(el->start), rptSize, psl->tStart, psl->tEnd) ;
        }
    int retroSize = psl->tEnd - psl->tStart - rptSize;
    if (retroSize < 0) retroSize = 1;
    /* calculate size gap and net at each level */
    for (el = elist; el != NULL ; el = el->next)
        {
        struct netSummary *np = el->val;
        int netPart = (el->end)-(el->start);
        float overlapThreshold = 0.01;
        if (sameString(np->type, "gap"))
            {
            if (((float)netPart/(float)retroSize) > overlapThreshold)
                {
                gapSize[np->level] += netPart ;
                if (np->qN > 0)
                    /*treat gaps as syntenic sequence */
                    gapSize[np->level] -= (float)netPart*(float)(np->qN)/(float)np->qSize; 
                else
                    gapSize[np->level] -= np->tR;
                if (gapSize[np->level] < 0)
                    gapSize[np->level] = 0;
                }
            assert(gapSize[np->level]>=0);
            verbose(5, "NET %s gapSize[%d]=%d qN %d tR %d netPart %d\n",
                    psl->qName, np->level, gapSize[np->level], np->qN, np->tR, netPart);
            }
        else 
            {
            int overlapPart = positiveRangeIntersection(psl->tStart, psl->tEnd, 
                el->start, el->end);
            netSize[np->level] += netPart;
            overlapSize[np->level] += overlapPart;
            verbose(4,"NET %s %s:%d netSize[%d]=%d retroSize %d part %d tRep %d qN %d overSize %d %s\n",
                    psl->qName, psl->tName, psl->tStart+1, np->level, netSize[np->level], 
                    retroSize, el->end-el->start, np->tR, np->qN, overlapSize[np->level], np->type );
            if (np->level > maxlevel && netSize[np->level] > MINNETSIZE)
                {
                maxlevel = np->level;
                if (maxlevel > 2)
                    {
                    if (gapSize[maxlevel-1] < MINGAP)
                        maxlevel -= 2;
                    verbose(4,"NET level bumped down from %d to %d gapSize[%d]=%d \n",
                            maxlevel+2, maxlevel, maxlevel-1,gapSize[maxlevel-1]);
                    }
                }
            }
        }
    if (maxlevel > 1 && rptSize > 0)
        if (gapSize[maxlevel-1] > netSize[maxlevel])
            {
            gapSize[maxlevel-1] -= rptSize;
            if (gapSize[maxlevel-1] < netSize[maxlevel])
                gapSize[maxlevel-1] = netSize[maxlevel]+1;
            verbose(4, "NET gap reduced by %d to %d\n",rptSize, gapSize[maxlevel-1]);
            assert(gapSize[maxlevel-1]);
            }
    /* no orthologous DNA, treat as break */
    if (maxlevel == 0)
        percentBreak = 99;
    /* orthologous DNA, no break */
    else if (netSize[maxlevel] > 100000)
        percentBreak = 0;
    else if (maxlevel==1)
        percentBreak = (float)overlapSize[maxlevel]*100/(float)netSize[maxlevel];
    /* sequence gaps should not be treated as breaks in orthology */
//    else if (2*gapSize[maxlevel-1] < netSize[maxlevel])
//        {
//        verbose(4, "NET gap %d < net %d zero break\n",gapSize[maxlevel-1], netSize[maxlevel]);
//        percentBreak = 0;
//        }
    else if (gapSize[maxlevel-1] < 100000)
        {
        percentBreak = (float)retroSize*100/(float)gapSize[maxlevel-1];
        verbose(4, "NET retro %d / gap %d  break %d netSize[%d]=%d\n",
                retroSize, gapSize[maxlevel-1], percentBreak, maxlevel, netSize[maxlevel]);
        if (percentBreak >= 130) percentBreak = 120;
        }
    else if (gapSize[maxlevel-1] >= 100000)
        {
        percentBreak = (float)overlapSize[maxlevel]*100/(float)retroSize;
        if (percentBreak > 100)
            percentBreak = 101;
        }
    verbose(3,"NET AFTER #score %s %s:%d maxlevel %d  overlap %d percentBreak=%d netSize[max] = %d \n",
            psl->qName, psl->tName, psl->tStart+1, maxlevel, overlapSize[maxlevel] , 
            percentBreak, netSize[maxlevel]);
    //slFreeList(&elist);
    }
assert (percentBreak < 130);
return percentBreak;
}

void outputLink(struct psl *psl, struct ucscRetroInfo *pg , struct dyString *reason, struct psl *bestParentPsl)
   /* char *type, char *bestqName, char *besttName, 
                int besttStart, int besttEnd, int maxExons, int geneOverlap, 
                char *bestStrand, int polyA, int polyAstart, int label, 
                int exonCover, int intronCount, int bestAliCount, 
                int tReps, int qReps, int overlapMouse, 
                struct genePred *kg, struct genePred *mgc, struct genePred *gp, 
                int alignGapCount, struct dyString *iString, int conservedIntrons, int maxOverlap) */
/* output bed record with pseudogene details and link to gene*/
{
struct axt *axt = NULL;
int pseudoScore = 0;
int bump = 0;
float maxOverlap = (float)pg->maxOverlap/(float)(psl->match+psl->misMatch+psl->repMatch)  ;
struct genePred *gp = NULL;
char name[512];
static struct genbankCds cds ;
/* calculate if pseudogene overlaps the syntenic diagonal with another species */
verbose(2,"NETFar1\n");
int overlapOrtho1 = netOverlap(psl, synHash);
verbose(2,"NETFar2\n");
int overlapOrtho2 = netOverlap(psl, syn2Hash);
verbose(2,"NETNear1\n");
int overlapOrtho3 = netOverlap(psl, syn3Hash);
int overlapOrtholog = max(overlapOrtho1, overlapOrtho2);

/* if retroGene has no overlap with orthologous DNA, assume it is a new insertion
   and give it a proxy orthlog score with 70% overlap with parent */
if (overlapOrtho1 == -10 && overlapOrtho2 == -10)
    overlapOrtholog = 70;
verbose(3,"%s netFar1 %d netFar2 %d composite %d netNear %d \n",
        psl->qName, overlapOrtho1, overlapOrtho2, overlapOrtholog, overlapOrtho3);
pg->milliBad = calcMilliScore(psl);
pg->axtScore = -1;
pg->type = reason->string;
safef(name, sizeof(name), "%s", psl->qName);
chopSuffixAt(name,'-');
cds.start = 0;
cds.end = 0;

/* map orf from parent and write genePred */
/* bump score more for knownGene, mgc and refSeq parents or parents valid ORF */
bump = getGenePred(pg, bestParentPsl, &gp);
if (gp == NULL)
    {
    verbose(5,"%s no genePred %s ",
            pg->name, pg->refSeq );
    if (bestParentPsl != NULL)
        verbose(5,"%s %s:%d-%d",
            bestParentPsl->qName, bestParentPsl->tName, bestParentPsl->tStart+1, bestParentPsl->tEnd);
    verbose(5,"\n");
    }

char *cdsStr = getCdsForAcc(name);
/* parses the CDS from cdsStr and stores the start (0-based coord) and end in 
 * the cds struct */ 
if (cdsStr != NULL)
    {
    genbankCdsParse(cdsStr, &cds);
    verbose(5,"cds %s %s cds %d-%d\n",
        name, pg->name, cds.start, cds.end);
    }
else
    verbose(5,"no cds %s\n",name);

/* convert to axt and compute axt score */
if (pg->label == PSEUDO || pg->label == EXPRESSED || pg->label == NOTPSEUDO) 
    {
    verbose(2,"pslToAxt %s \n",psl->qName);
    axt = pslToAxt(psl, &cds);
    
    if (axt != NULL)
        {
        pg->axtScore = axtScoreFilterRepeats(axt, ss);
        verbose(6,"qseq %s %s\n",psl->qName, axt->qSym);
        }
    else
        verbose(2,"no axt Score %d qname %s\n",pg->axtScore, psl->qName);
    }

/* all weighted features are scaled to range from 0 to 1000 */
assert(psl->qSize > 0);
float scaledMilliBad = (max(pg->milliBad,600)-600)*2.5;
int scaledCoverage =(pg->coverage/100.0)*(1.0-((float)(psl->qSize-psl->qEnd)/(float)psl->qSize))*log(psl->qSize)*100.0;
pseudoScore = ( wt[0]*scaledMilliBad
                + wt[1]*(log(pg->exonCover+1)/log(2))*600       // *1
                //+ wt[2]*log(pg->axtScore>0?pg->axtScore:1)*70 
                + wt[2]*(((log(pg->axtScore>0?pg->axtScore:1)/log(2))*170)-1000)   // *.10
                + wt[3]*(log(pg->polyAlen+2)*200)                               // *.20
                + wt[4]*(overlapOrtholog*10)                            // *.70
                + wt[5]*(((log(pg->processedIntrons > 0 ? (pg->processedIntrons)+1 : 1))/log(2))*1200 ) // *1.00
                //- wt[5]*pow(pg->conservedIntrons,0.5)*2000 
                - wt[6]*pow(pg->intronCount,0.5)*1200                    // *-1.00
                //- wt[6]*pow(pg->intronCount,0.5)*2000 
                - wt[7]*(maxOverlap*300)                        // *-.50
                + wt[8]*scaledCoverage //0.50
                - wt[9]*(pg->tReps*10)  // *-1.00
                - wt[10]*(pg->alignGapCount)  //*-1.00
                ) / ScoreNorm + bump;
//wt[0] = 0; wt[1] = 0.85; wt[2] = 0.1; wt[3] = 0.2; wt[4] = 0.7; 
//wt[5] = 1; wt[6] = 1  ; wt[7] = 0.5; wt[8] = 0.5; wt[9] = 1; wt[10] = 1;
verbose(1, "##polyA Len: %d, ScoreNorm: %4.1f \n", pg->polyAlen, ScoreNorm);
verbose(1,"##score %d %s %s:%d-%d scr %d milbad %d +%4.1f xon %d +%4.1f retainS %d ax +%4.1f pA +%4.1f net +%4.1f max (%d, %d) procIntrons %d +%4.1f in.cnt %d -%4.1f ov -%4.1f  cov %4.2f*qCov %4.2f log(qSize)%4.1f qSize %d +%4.1f tRep -%4.1f alignGapCount %d -%4.1f norm/%d %s bump+%d \n", 
                pg->label, psl->qName, psl->tName, psl->tStart+1, psl->tEnd, pseudoScore, 
                pg->milliBad, scaledMilliBad/ScoreNorm,
                pg->exonCover,
                wt[1]*(log(pg->exonCover+1)/log(2))*600/ScoreNorm , 
                pg->conservedSpliceSites,
                wt[2]*(((log(pg->axtScore>0?pg->axtScore:1)/log(2))*170)-1000)/ScoreNorm,
                wt[3]*(log(pg->polyAlen+2)*200)/ScoreNorm ,
                wt[4]*overlapOrtholog*10/ScoreNorm , overlapOrtho1, overlapOrtho2,
                pg->processedIntrons,
                wt[5]*(((log(pg->processedIntrons > 0 ? (pg->processedIntrons)+1: 1))/log(2))*1200)/ScoreNorm ,
                pg->intronCount, 
                wt[6]*pow(pg->intronCount,0.5)*1200/ScoreNorm ,
                wt[7]*(maxOverlap*300)/ScoreNorm,
                (float)pg->coverage/100.0, 1.0-((float)(psl->qSize-psl->qEnd)/(float)psl->qSize), log(psl->qSize), psl->qSize,
                wt[8]*scaledCoverage/ScoreNorm,
                wt[9]*(pg->tReps*10)/ScoreNorm, 
                pg->alignGapCount,
                wt[10]*pg->alignGapCount/ScoreNorm,
                ScoreNorm, pg->type, bump
                ) ;
verbose(1,"###score %d %s %s:%d-%d scr %d milbad +%4.1f ParentExonCover +%4.1f axtScore +%4.1f pA +%4.1f net +%4.1f PInt +%4.1f ic -%4.1f ov -%4.1f cov +%4.1f tRep -%4.1f alignGapCount -%4.1f norm/%d %s bump+%d \n", 
                pg->label, psl->qName, psl->tName, psl->tStart+1, psl->tEnd, pseudoScore, 
                wt[0]*scaledMilliBad/ScoreNorm,
                wt[1]*(log(pg->exonCover+1)/log(2))*600/ScoreNorm , 
                wt[2]*(((log(pg->axtScore>0?pg->axtScore:1)/log(2))*170)-1000)/ScoreNorm,
                wt[3]*(log(pg->polyAlen+2)*200)/ScoreNorm ,
                wt[4]*overlapOrtholog*10/ScoreNorm , 
                wt[5]*(((log(pg->processedIntrons > 0 ? (pg->processedIntrons)+1.00 : 1.00))/log(2))*1200.0)/ScoreNorm ,
                wt[6]*pow(pg->intronCount,0.5)*1200/ScoreNorm ,
                wt[7]*(maxOverlap*300)/ScoreNorm,
                wt[8]*scaledCoverage/ScoreNorm,
                wt[9]*(pg->tReps*10)/ScoreNorm, 
                wt[10]*pg->alignGapCount/ScoreNorm,
                ScoreNorm, pg->type, bump
                ) ;
if (pseudoScore > 0)
    pg->score = pseudoScore;
else 
    pg->score = 0;

/* write pseudogene to gene alignment in psl and axt format */
if ((pg->label == PSEUDO || pg->label == EXPRESSED) ) 
    {
    pslTabOut(psl, pseudoFile);
    if (axt != NULL && pg->score > scoreThreshold)
        {
        axtWrite(axt, axtFile);
        axtFree(&axt);
        }
    else
        verbose(3, "no axt, type = %d score = %d\n",pg->label, pg->score);
    }

//wt[0] = 0.3; wt[1] = 0.85; wt[2] = 0.7; wt[3] = 0.4; wt[4] = 0.3; 
//wt[5] = 0;   wt[6] = 1  ;  wt[7] = 0.5; wt[8] = 1;   wt[9] = 1;
outputNoLinkScore(psl, pg, overlapOrtho1, overlapOrtho2, overlapOrtho3);
}

int intronFactor(struct psl *psl, struct hash *rmskHash, struct hash *trfHash, int *baseCount)
/* Figure out the approximate number of introns.  
 * An intron in this case is just a gap of 0 bases in query and
 * maxBlockGap or more in target with repeats masked. */
{
int i, blockCount = psl->blockCount  ;
int ts, qs, te, qe, sz, tsRegion, teRegion;
struct binKeeper *bk;
int intronCount = 0, tGap;

*baseCount = 0;
/* psl should not be NULL and if there are <=1 blocks then return 0 */
assert (psl != NULL);
if (blockCount <= 1)
    return 0;
/* get size of first block and the mRNA and genomic coords for the end of block */
sz = psl->blockSizes[0];
qe = psl->qStarts[0] + sz;
te = psl->tStarts[0] + sz;
tsRegion = psl->tStarts[0];
for (i=1; i<blockCount; ++i)
    {
    int trf = 0;
    int reps = 0, regionReps = 0;
    struct binElement *el, *elist;
    qs = psl->qStarts[i];
    ts = psl->tStarts[i];
    teRegion = psl->tStarts[i] + psl->blockSizes[i];

    /* if repeatMasker hash and trf hash are not NULL then find the 
     binKeeper for tName (chrom) and if this is not NULL then get the 
     list of trf repeats for the range tStart to End and count the total
     number of bases of trf repeats the intersect with gap. */
    if (rmskHash != NULL)
        {
        if (trfHash != NULL)
            {
            bk = hashFindVal(trfHash, psl->tName);
            if (bk != NULL)
                {
                elist = binKeeperFindSorted(bk, psl->tStart , psl->tEnd ) ;
                trf = 0;
                for (el = elist; el != NULL ; el = el->next)
                    {
                    //bed = el->val;
                    trf += positiveRangeIntersection(te, ts, el->start, el->end);
                    }
                slFreeList(&elist);
                }
            }
        /* then do the same for RepeatMasker repeats counting the total 
           number of repeat bases that intersect with the gap. */
        bk = hashFindVal(rmskHash, psl->tName);
        if (bk != NULL)
            {
            elist = binKeeperFindSorted(bk, tsRegion , teRegion ) ;
            for (el = elist; el != NULL ; el = el->next)
                {
                /* count repeats in the gap */
                reps += positiveRangeIntersection(te, ts, el->start, el->end);
                /* count repeats in the gap  and surrounding exons  */
                regionReps += positiveRangeIntersection(tsRegion, teRegion, el->start, el->end);
                }
            if (elist != NULL)
                slFreeList(&elist);
            }
        }
    /* don't subtract repeats if the entire region is masked */
    if ( regionReps + 10 >= teRegion - tsRegion )
        reps = 0;
    /* why multiple repeat bases by 1.5? */
    tGap = ts - te - reps*1.5 - trf;
    /* why is 2*abs(qs-qe) used */
    /* only counts it if qGap small compared to tGap, if same then some
       other insertion */
    if (2*abs(qs - qe) <= tGap && tGap >= maxBlockGap ) 
        /* don't count q gaps in mrna as introns , if they are similar in size to tGap*/
        {
        /* increment intron count and baseCount */
        intronCount++;
        *baseCount += tGap - (2*abs(qs-qe));
        verbose(6, "YES ");
        }
    verbose(6,"%s:%d (%d < tGap %d) and > %d  qs %d qe %d   ts %d te %d reps %d trf %d \n",
            psl->tName, psl->tStart+1,2*abs(qs-qe),tGap, maxBlockGap,  qs,qe, ts, te, reps, trf);
    assert(psl != NULL);
    /* values for next round */
    sz = psl->blockSizes[i];
    qe = qs + sz;
    te = ts + sz;
    tsRegion = psl->tStarts[i];
    }
verbose(6, "old intronCount %d %s %s bases %d\n",intronCount, psl->qName, psl->tName, *baseCount);
return intronCount;
}

int sizeFactor(struct psl *psl, struct hash *rmskHash, struct hash *trfHash)
/* Return a factor that will favor longer alignments. An intron is worth 3 bases...  */
{
int score;
int bases;
if (ignoreSize) return 0;
assert(psl != NULL);
/* how did you decide on this score? */
score = 4*round(sqrt(psl->match + psl->repMatch/4));
/* if noIntrons not specified then take these into account */
if (!noIntrons)
    {
    int bonus = intronFactor(psl, rmskHash, trfHash, &bases) * 3;
    if (bonus > 10) bonus = 10;
    score += bonus;
    }
return score;
}

int calcSizedScore(struct psl *psl, struct hash *rmskHash, struct hash *trfHash)
/* Return score that includes base matches and size. */
{
int score = calcMilliScore(psl) + sizeFactor(psl, rmskHash, trfHash);
return score;
}

boolean closeToTop(struct psl *psl, int *scoreTrack , int milliScore) /*struct hash *rmskHash, struct hash *trfHash)*/
/* Returns TRUE if psl is near the top scorer for at least minNearTopSize bases. */
{
int threshold = round(milliScore * (1.0+nearTop));
int i, blockIx;
int start, size, end;
int topCount = 0;
char strand = psl->strand[0];

verbose(5,"%s:%d milliScore %d, threshold %d\n", psl->tName, psl->tStart+1, milliScore, threshold);
for (blockIx = 0; blockIx < psl->blockCount; ++blockIx)
    {
    start = psl->qStarts[blockIx];
    size = psl->blockSizes[blockIx];
    end = start+size;
    if (strand == '-')
	reverseIntRange(&start, &end, psl->qSize);
    verbose(5,"block %d threshold=%d s-e:%d-%d. ",blockIx, threshold, start, end);
    for (i=start; i<end; ++i)
	{
        verbose(6,"s=%d tc=%d ",scoreTrack[i],topCount);
        /* not sure I see how this gives values close to the top scorer */
        /* also means that the threshold >= topscorer that has that base */
	if (scoreTrack[i] <= threshold)
	    {
	    if (++topCount >= minNearTopSize)
                {
		return TRUE;
                }
	    }
	}
    }
return FALSE;
}

int scoreWindow(char c, char *s, int size, int *score, int *start, int *end, int match, int misMatch, float threshold)
/* calculate score array with score at each position in s, match to char c adds 1 to score, mismatch adds -1 */
/* index of max score is returned , size is size of s */
/* number char c within window must have density greater than threshold */
{
int i=0, j=0, max=2, count = 0; 

*end = 0;

assert(match >= 0);
assert(misMatch <= 0);
/* go through the sequence one base at a time and score each one, adding 1 to
 * the previous scores if it matches char c and -1 if it doesn't. this will 
 * find high scoring runs of the character c. */
for (i=0 ; i<size ; i++)
    {
    int prevScore = (i > 0) ? score[i-1] : 0;

    if (toupper(s[i]) == toupper(c) )
        score[i] = prevScore+match;
    else
        score[i] = prevScore+misMatch;
    /* if the score is greater than max (2 initially), then
     * reset max to this score and set the end of the run of char c to 
     * this position and go back to find the beginning of the run of c chars */
    if (score[i] >= max)
        {
        max = score[i];
        *end = i;
        for (j=i ; j>=0 ; j--)
            if (score[j] == 0)
                {
                *start = j+1;
                break;
                }
        }
    verbose(6,"max is %d %d - %d score=%d i=%d\n",max, *start, *end, score[i],i);
    if (score[i] < 0) 
        score[i] = 0;
    }
assert (*end < size);
/* traceback to find start of the run of char c*/
for (i=*end ; i>=0 ; i--)
    if (score[i] == 0)
        {
        *start = i+1;
        break;
        }

for (i=*start ; i<=*end ; i++)
    {
    assert (i < size);
    if (toupper(s[i]) == toupper(c) )
        count++;
    }
if (*end != 0 )
    {
    if (count/(*end-*start) < threshold)
        verbose(4,"below threshold count %d %d - %d %5.2f < %5.2f \n",count, *start, *end, (float)count/(float)(*end-*start), threshold);
    }
else
    count = 0;
return count;
}
   
int polyACalc(int start, int end, char *strand, int tSize, char *chrom, int region, 
        int *polyAstart, int *polyAend, float divergence)
/* get size of polyA tail in genomic dna , count bases in a 
 * sliding window in region of the end of the sequence*/
/* use the divergence from the gene as the scoring threshold */
{
int seqSize;
struct dnaSeq *seq = NULL;
int count = 0;
int length = 0;
int score[POLYAREGION+1], pStart = 0; 
int match = 1;
int misMatch = -1;
float threshold = divergence/100;
int tOffset = 0;
int seqStart = strand[0] == '+' ? end - region/2 : start - region/2;
/* use mixed case for sequence, repeats in lower case */
boolean doMask = TRUE;
int retFullSeqSize;

seqSize = twoBitSeqSize(genomeSeqFile, chrom);
assert(region > 0);
assert(end != 0);
*polyAstart = 0 , *polyAend = 0;
assert (seqSize == tSize);
if (seqStart < 0) seqStart = 0;
if (seqStart + region > seqSize) region = seqSize - seqStart;

verbose(2,"search for polyA at %s:%d-%d len %d start %d tSize %d match %d misMatch %d %f\n",
        chrom, seqStart+1, seqStart+region,region, start, tSize,match,misMatch,divergence);
if (region == 0)
    return 0;
assert(region > 0);
seq = twoBitReadSeqFrag(genomeSeqFile, chrom, seqStart, seqStart+region);

if (strand[0] == '+')
    {
    assert (seq->size <= POLYAREGION);
verbose(4,"\n + range (0-based start)=%d %d %s %s\n",seqStart, seqStart+region, seq->dna, chrom );
    count = scoreWindow('A',seq->dna,seq->size, score, polyAstart, polyAend, match, misMatch, threshold);
    }
else
    {
    assert (seq->size <= POLYAREGION);
verbose(4,"\n - range=%d %d %s %s\n",seqStart, seqStart+region, seq->dna, chrom );
    count = scoreWindow('T',seq->dna,seq->size, score, polyAend, polyAstart, match, misMatch, threshold);
    }
/* add genomic base position of start of window to search for polyA tail */
pStart += seqStart;
/* get genomic coordinate of polyAstart */
*polyAstart += seqStart;
/* get genomic coordinate of polyAEnd */
*polyAend += seqStart;
/* calculate length of polyA tail */
length = (strand[0]=='+'?*polyAend-*polyAstart:(*polyAstart)-(*polyAend))+1;
verbose(4,"\npolyA is %d from end. seqStart=%d %s polyS/E %d-%d exact matches %d length %d \n",
        strand[0]=='+'?*polyAstart-seqStart:*polyAend-seqStart+seq->size, 
	seqStart, seq->dna+(*polyAstart)-seqStart, 
	*polyAstart, *polyAend, count, length);
freeDnaSeq(&seq);
return count;
}

void pslMergeBlocks(struct psl *psl, struct psl *outPsl,
                       int insertMergeSize)
/* merge together blocks separated by small inserts. */
{
int iBlk, iExon = -1, blockIx;
int startIdx, stopIdx, idxIncr;
int ii = 0;

assert(outPsl!=NULL);
outPsl->qStarts = needMem(psl->blockCount*sizeof(unsigned));
outPsl->tStarts = needMem(psl->blockCount*sizeof(unsigned));
outPsl->blockSizes = needMem(psl->blockCount*sizeof(unsigned));

verbose(4,"MERGE %d blocks on %c match %d misMatch %d %s %s \n",
        psl->blockCount, psl->strand[0], psl->match, psl->misMatch, psl->qName, psl->tName);
verbose(4, "before merge %s %c %s %d blocks  ",
        psl->qName, psl->strand[0], psl->tName, psl->blockCount);
for (blockIx = 0; blockIx < psl->blockCount; ++blockIx)
    verbose(4, "%d-%d, ",psl->qStarts[blockIx], 
            psl->qStarts[blockIx]+psl->blockSizes[blockIx]);
verbose(4, "\n");
verbose(4, "before merge %s %c %s %d blocks  ",
        psl->qName, psl->strand[0], psl->tName, psl->blockCount);
for (blockIx = 0; blockIx < psl->blockCount; ++blockIx)
    verbose(4, "%d-%d, ",psl->tStarts[blockIx], 
            psl->tStarts[blockIx]+psl->blockSizes[blockIx]);
verbose(4, "\n");

    startIdx = 0;
    stopIdx = psl->blockCount;
    idxIncr = 1;

/* go through the aligned blocks of the PSL */
for (iBlk = startIdx; iBlk != stopIdx; iBlk += idxIncr)
    {
    int  tStart = psl->tStarts[iBlk];
    int  qStart = psl->qStarts[iBlk];
    int  size = psl->blockSizes[iBlk];
    int  qEnd = psl->qStarts[iBlk]+size;
    int  outQStart ;
    int  outQEnd ;
    int tdiff = 0;
    int qdiff = 0;
    outQStart = outPsl->qStarts[iExon >= 0 ? iExon :0];
    outQEnd = outPsl->qStarts[iExon >= 0 ? iExon : 0]+outPsl->blockSizes[iExon >= 0 ? iExon : 0];
    if (psl->strand[0] == '-')
        {
        reverseIntRange(&qStart, &qEnd, psl->qSize);
        reverseIntRange(&outQStart, &outQEnd, psl->qSize);
        }
    if (iExon >= 0)
        {
        tdiff = abs(tStart - (outPsl->tStarts[iExon]+outPsl->blockSizes[iExon]));
        qdiff = abs(qStart - (outQEnd));
        }
    if (iExon < 0 || (tdiff > insertMergeSize))
        {
        iExon++;
        verbose(5, " init or Not merge tdiff %d > %d qdiff %d=%d-%d tStart %d out tStarts[%d] %d outPsl->size %d\n",
                tdiff, insertMergeSize, qdiff, qStart, outQEnd, tStart, iExon, outPsl->tStarts[iExon], outPsl->blockSizes[iExon]);
        verbose(5,"  init or Not merge %s[%d] new q %d t %s %d %d size %u to %d \n", 
            psl->qName, iExon, qStart,
            psl->tName, psl->tStarts[iBlk], tStart, outPsl->blockSizes[iBlk], size);
        outPsl->tStarts[iExon] = tStart;
        if (psl->strand[0] == '-')
            reverseIntRange(&qStart, &qEnd, psl->qSize);
        outPsl->qStarts[iExon] = qStart;
        if (size > psl->qSize)
            assert(size <= psl->qSize);
        outPsl->blockSizes[iExon] = size; 
	}
    else
        {
        /* why is this qdiff not tdiff */
        outPsl->blockSizes[iExon] += size + qdiff ;
        verbose(5, " MERGE tdiff %d < %d tStart %d out.tStarts[%d] %d size %d outPsl->blkSmize %d tdiff %d qdiff %d tStart %d\n",
                tdiff, insertMergeSize, tStart, iExon, outPsl->tStarts[iExon], size, outPsl->blockSizes[iExon], tdiff, qdiff, psl->tStart);
        }
verbose(5, "during merge %s %c %s %d blocks  ",
        psl->qName, psl->strand[0], psl->tName, iExon);
for (blockIx = 0; blockIx < iExon; ++blockIx)
    verbose(5, "%d-%d, ",outPsl->tStarts[blockIx], 
            outPsl->tStarts[blockIx]+outPsl->blockSizes[blockIx]);
verbose(5, "\n");
    }

verbose(5, "before last step merge %s %d ",psl->qName, outPsl->blockCount);
/* calculate ii as the sum of block sizes in outPsl this variable is not used
 * elsewhere */
for (blockIx = 0; blockIx < outPsl->blockCount; ++blockIx)
    {
    verbose(5, "%d-%d, ",outPsl->qStarts[blockIx], 
            outPsl->qStarts[blockIx]+outPsl->blockSizes[blockIx]);
    ii+=outPsl->blockSizes[blockIx];
    }
verbose(5, "\n");

/* do you need to revese the int range for a PSL? */
if (psl->strand[0] == '-')
    for (iBlk = 0; iBlk != outPsl->blockCount; iBlk += 1)
        {
        int qStart = outPsl->qStarts[iBlk];
        int qEnd = outPsl->qStarts[iBlk]+outPsl->blockSizes[iBlk];
        reverseIntRange(&qStart, &qEnd, psl->qSize);
        outPsl->qStarts[iBlk] = qStart;
        if (qStart < 0)
            assert(qStart >= 0);
        if (qStart+outPsl->blockSizes[iBlk] > outPsl->qSize)
            assert(qStart+outPsl->blockSizes[iBlk] <= outPsl->qSize);
        }     
/* set other fields in the outPsl, mostly just copying over from psl */ 
outPsl->blockCount = iExon+1;
outPsl->match = psl->match;
outPsl->misMatch = psl->misMatch;
outPsl->repMatch = psl->repMatch;
outPsl->nCount = psl->nCount;
outPsl->qNumInsert = psl->qNumInsert;
outPsl->qBaseInsert = psl->qBaseInsert;
outPsl->tNumInsert = psl->tNumInsert;
outPsl->tBaseInsert = psl->tBaseInsert;
strcpy(outPsl->strand, psl->strand);
outPsl->qName = cloneString(psl->qName);
outPsl->qSize = psl->qSize;
outPsl->qStart = psl->qStart;
outPsl->qEnd = psl->qEnd;
outPsl->tName = cloneString(psl->tName);
outPsl->tSize = psl->tSize;
outPsl->tStart = psl->tStart;
outPsl->tEnd = psl->tEnd;
for (iBlk = 0; iBlk != psl->blockCount; iBlk += 1)
    verbose(5, "%d,",psl->qStarts[iBlk]);
verbose(5,"  ");
for (iBlk = 0; iBlk != psl->blockCount; iBlk += 1)
    verbose(5, "%d,",psl->blockSizes[iBlk]);
verbose(5," to ");
for (iBlk = 0; iBlk != outPsl->blockCount; iBlk += 1)
    verbose(5, "%d,",outPsl->qStarts[iBlk]);
verbose(5,"  ");
for (iBlk = 0; iBlk != outPsl->blockCount; iBlk += 1)
    verbose(5, "%d,",outPsl->blockSizes[iBlk]);
verbose(5,"\n");

verbose(5,"AFTER MERGE %d blocks on %c match %d misMatch %d %s %s sum of blocks %d\n",
        outPsl->blockCount, outPsl->strand[0], outPsl->match, outPsl->misMatch, outPsl->qName, outPsl->tName, ii);
verbose(5, "AFTER merge %s %c %s %d blocks  ",
        outPsl->qName, outPsl->strand[0], outPsl->tName, outPsl->blockCount);
for (blockIx = 0; blockIx < outPsl->blockCount; ++blockIx)
    verbose(5, "%d-%d, ",outPsl->tStarts[blockIx], 
            outPsl->tStarts[blockIx]+outPsl->blockSizes[blockIx]);
verbose(5, "\n");
}

bool isRepeat (char *chrom, int is, int ie, struct hash *rmskHash, int *retRepCount)
/* check for repeats from intron*/
{
struct binKeeper *bk = NULL;
struct binElement *elist = NULL, *el = NULL;
int reps = 0;
if(ie < is) return FALSE;

if (rmskHash != NULL)
    {
    bk = hashFindVal(rmskHash, chrom);
    if (bk != NULL)
        {
        elist = binKeeperFindSorted(bk, is, ie) ;
        for (el = elist; el != NULL ; el = el->next)
            {
            reps += positiveRangeIntersection(is, ie, el->start, el->end);
            verbose(5,"    isRep? chkRep %s:%d-%d reps %d ratio %f < 0.7\n",chrom,is+1, ie, reps,(float)reps/(float)(ie-is) );
            }
        slFreeList(&elist);
        if (reps > ie-is)
            {
            reps = ie-is;
            verbose(5,"Warning: too many reps %d in %s:%d-%d size %d\n",
                            reps, chrom, is+1, ie, ie-is);
            }
        verbose(4,"     check if Intron is repeat: ratio = %f rep %d intron %d \n",
                (float)reps/(float)(ie-is),reps, ie-is);
        if ((float)reps/(float)(ie-is) > repsPerIntron)
            {
            verbose(3,"     Intron is repeat: ratio = %f \n",(float)reps/(float)(ie-is));
            if (retRepCount != NULL)
                *retRepCount = reps;
            return TRUE;
            }
        }
    }
if (retRepCount != NULL)
    *retRepCount = reps;
return FALSE;
}

/* get overlap of q blocks and set the number of overlapping blocks to numBlocks */
int getQOverlap(struct psl *psl, int start, int end, int *numBlocks)
{
int total = 0;
int i;
int blocks = 0;
for (i = 0 ; i < psl->blockCount ; i++)
    {
    int qs = psl->qStarts[i];
    int qe = psl->qStarts[i] + psl->blockSizes[i];
    int te = psl->tStarts[i] + psl->blockSizes[i];
    int overlap = 0;
    int tdiff = 0;
    int tsNext = 999999;
    float coverage = 0;
    if (i < psl->blockCount -1)
        {
        tsNext = psl->tStarts[i+1];
        tdiff = tsNext - te;
        }
    else
        tdiff = 9999999;
    /* merges blocks together if target gap is less than a certain size */
    while (tdiff < maxBlockGap && i< (psl->blockCount)-1)
        {
        i++;
        te = psl->tStarts[i] + psl->blockSizes[i];
        qe = psl->qStarts[i] + psl->blockSizes[i];
        if (i < psl->blockCount -1)
            {
            tsNext = psl->tStarts[i+1];
            tdiff = tsNext - te;
            }
        else
            tdiff = 9999999;
        }
    if (psl->strand[0] == '-')
        reverseIntRange(&qs, &qe, psl->qSize);
    overlap = positiveRangeIntersection(start, end, qs, qe);
    coverage = (float)overlap/(float)(qs-qe);
    total += overlap;
    if (overlap > 25 || coverage > 0.4)
        {       
        blocks++;
        verbose(4,"%s QBLOCK %d-%d overlaps %d-%d %4.2f > 0.40 %d bp > 25\n",
                psl->qName, start, end, qs, qe, coverage, overlap);
        }
    }
*numBlocks = blocks;
return total;
}

/* count Spliced introns and retain Introns */
void calcIntrons(struct psl *psl, int maxBlockGap, struct psl *nestedPsl, 
        int *retProcessedIntrons, int *retIntronCount, int *retBlockCover, int *pseudoExonCount)
{
int i, qsStart = 0;
int blocksCovered = 0;
int totalBases = 0;
int totalExons = 0;
int intronsPresent = 0;
int intronsPresentBases = 0;
int intronsSpliced = 0;
int exonCount = 0;
/* go through the aligned blocks of the PSL */
for (i = 0 ; i < psl->blockCount ; i++)
    {
    int qs = psl->qStarts[i];
    int qe = psl->qStarts[i] + psl->blockSizes[i];
    int te = psl->tStarts[i] + psl->blockSizes[i];
    int tdiff = 0;
    int qdiff = 0;
    int cumTdiff = 0;
    int oqs, oqe;
    int tsNext = 999999;
    int qsNext = 0;
    float coverage = 0;
    int bases = 0;
    int oldte = te;
    int gapRatio = 0;
    int repCnt = 0;
    bool isRep = FALSE;
    if (i < psl->blockCount -1)
        {
        tsNext = psl->tStarts[i+1];
        qsNext = psl->qStarts[i+1];
        /* get size of gap between current block and next for target */
        tdiff = tsNext - te;
        verbose(6, "%s tdiff %d = tsNext %d - te %d\n",psl->qName, tdiff, tsNext, te);
        /* get size of gap between current block and the next for query */
        qdiff = qsNext - qe;
        }
    else
        tdiff = 9999999;
    /* set cumulative tdiff to tidff */
    cumTdiff = tdiff;
    qsStart = qs;
    /* increment exon count */
    exonCount++;
/* combine blocks that are close together */
    while (tdiff < maxBlockGap && i< (psl->blockCount)-1)
        {
        i++;
        te = psl->tStarts[i] + psl->blockSizes[i];
        qe = psl->qStarts[i] + psl->blockSizes[i];
        if (i < psl->blockCount -1)
            {
            tsNext = psl->tStarts[i+1];
            qsNext = psl->qStarts[i+1];
            tdiff = tsNext - te;
            qdiff = qsNext - qe;
            verbose(6, "%s combining block %d out of %d tgap removed %d max allowed is %d\n",
                        psl->qName, i, psl->blockCount, tdiff, maxBlockGap);
            verbose(6, "%s tdiff %d = tsNext %d - te %d\n",psl->qName, tdiff, tsNext, te);
            }
        else
            tdiff = 9999999;
        cumTdiff += tdiff;
        oldte = te;
        }
    /* get number of repeats bases in this gap region and calculate the
       number of non-repeat bases in the gap */
    isRep = isRepeat(psl->tName, te,tsNext,rmskHash, &repCnt);
    verbose(6, "%s %s:%d-%d exon %d reps = %d tdiff %d rep ratio %f isRep %d\n",
            psl->qName, psl->tName, te, tsNext, i, repCnt, tdiff, (float)repCnt/tdiff, isRep);
    tdiff -= repCnt;
    verbose(6, "%s tdiff %d -= repCnt  %d \n",psl->qName, tdiff, repCnt);
    //assert (repCnt >= 0);
    /* if the psl query is on the - strand, then get query start and end
       relative to the - strand */ 
    oqs = qs; oqe = qe;
    if (psl->strand[0] == '-')
        reverseIntRange(&oqs, &oqe, psl->qSize);
    /* calculate the gap ratio as a percentage of the proportion of 
     * query gap over target gap */
    if (tdiff > 0)
        gapRatio = qdiff*100/tdiff;
    else
        gapRatio = 100;
    verbose(5, "%d-%d/%d  ",oqs,oqe,gapRatio); 
    /* If this gap is an intron increment the intron count and add the 
     * number of bases to the cumulative total of bases */
    if (gapRatio < 30 && tdiff > minIntronSize && 
            i < (psl->blockCount)-1 && !isRep)
        {
        int bases = (tdiff-qdiff > 0) ? tdiff - qdiff : 0;
        verbose(5, "INTRON FOUND %s [q/t gapRatio %d < 30%%] qdiff %d [tdiff %d > %d]  q %d-%d oqs/e %d-%d t %d-%d block %d bases %d  ",
                psl->qName, gapRatio, qdiff, tdiff, minIntronSize, qs, qe, oqs, oqe, psl->tStarts[i], te, i, bases);
        intronsPresent++;
        intronsPresentBases += bases;
        }
    /* if a nestedPsl is provided get the number of overlapping bases 
       and keep a cumulative count of the bases and keep a count of the total
       number of exons and spliced introns */
    if (nestedPsl != NULL)
        {
        int numBlocks = 0;
        bases = getQOverlap(nestedPsl, oqs, oqe, &numBlocks); 
        coverage = (float)bases/(float)(oqe-oqs);
        verbose(5, " qBlocks overlap %d - %d = %d %4.2f ", oqs, oqe, bases, coverage);
        totalBases += bases;
        totalExons++;
        if (coverage > 0.20)
            blocksCovered++;
        verbose(5,"/%d, ", numBlocks-1);
        intronsSpliced += numBlocks-1;
        }
    verbose (4, "%s intronsSpliced %d introns %d exons Covering alignment %d\n",
            psl->qName, intronsSpliced, intronsPresent, blocksCovered);
    }
if (retProcessedIntrons != NULL)
    *retProcessedIntrons = intronsSpliced;
if (retIntronCount != NULL)
    {
    *retIntronCount = intronsPresent;
verbose (4, "%s FINAL intronsSpliced %d introns %d exons Covering alignment %d\n",
            psl->qName, *retProcessedIntrons, *retIntronCount, *retBlockCover);
    }
if (pseudoExonCount != NULL)
    *pseudoExonCount = exonCount;
if (retBlockCover != NULL)
    *retBlockCover = blocksCovered;
}

bool pslInRange(struct psl *psl, int i)
/* check if block index is out of range of psl */
{
if (i < 0) 
    return FALSE;
if (i >= psl->blockCount)
    return FALSE;
return TRUE;
}

void printPslQBlocks(struct psl *psl, bool printEnd)
{
int blockIx;

verbose(4, ":%s: ",psl->qName);
for (blockIx = 0; blockIx < psl->blockCount; ++blockIx)
    {
    verbose(4, "%d",psl->qStarts[blockIx]); 
    if (printEnd)
        verbose(4, "-%d",psl->qStarts[blockIx]+psl->blockSizes[blockIx]);
    verbose(4,", ");
    }
verbose(4, "\n");
}

int countRepeatOverlap(char *name, int start, int end)
/* count repeats inside a range */
{
int teReps = 0;
struct binKeeper *bk = NULL;
struct binElement *elist = NULL, *el = NULL;
if (end < start)
    return 0;
if (rmskHash != NULL)
    {
    bk = hashFindVal(rmskHash, name);
    if (end < start)
        errAbort("name %s end %d < start %d\n",
                name, end, start);
    assert(start <= end);
    if (bk != NULL)
        {
        elist = binKeeperFindSorted(bk, start, end) ;
        for (el = elist; el != NULL ; el = el->next)
            {
            teReps += positiveRangeIntersection(start, end, el->start, el->end);
            }
        slFreeList(&elist);
        }
    }
return teReps;
}

int countRetainedSpliceSites(struct psl *target, struct psl *query, int spliceDrift)
/* count number of splice sites from parent gene that are within spliceDrift 
 * bases of splice site in retro. splice site is measured relative to the 
 * query coordinates */
{
int count = 0, i;
bool swapT = FALSE, swapQ = FALSE;
//struct psl *targetM = NULL, *queryM = NULL;

if (target == NULL || query == NULL)
    return 0;

verbose(5, "countRetainSS t strand %s %s ",target->strand, target->tName);
printPslQBlocks(target, FALSE);
verbose(5, "countRetainSS q strand %s %s ",query->strand, query->tName);
printPslQBlocks(query, FALSE);
/* if either the target or query PSL strand is - then reverse coordinates so relative to - strand */
if (target->strand[0] == '-')
    {
    pslRc(target);
    swapT = TRUE;
    }
if (query->strand[0] == '-')
    {
    pslRc(query);
    swapQ = TRUE;
    }
verbose(4, "countRetainSS target blocks %d %s ",target->blockCount, target->strand);
printPslQBlocks(target, FALSE);
verbose(4, "countRetainSS query blocks %d %s ",query->blockCount,query->strand);
printPslQBlocks(query, FALSE);
/* query size for target PSL and query PSL should be the same else abort with warning message */
if (target->qSize != query->qSize)
    {
    abortAtEnd = TRUE;
    warn("size mismatch parent %s %s:%d-%d %d != %s %d %s:%d-%d \n",
            target->qName, 
            target->tName, target->tStart+1, target->tEnd,target->qSize,
            query->qName, query->qSize, 
            query->tName, query->tStart+1, query->tEnd);
    }

for (i = 0 ; i < target->blockCount ; i++)
    {
    int qe = target->qStarts[i] + target->blockSizes[i];
    int qs = target->qStarts[i];
    int te = target->tStarts[i] + target->blockSizes[i];
    int ts = target->tStarts[i];
    int j;
    bool negStrand = query->strand[1] == '-';
    if (target->strand[1] == '-')
        reverseIntRange(&ts, &te, target->tSize);
    for (j = 0 ; j < query->blockCount ; j++)
        {
        int offset = (negStrand) ? -1 : 1;
        int firstBlock = (negStrand) ? query->blockCount-1 : 0;
        int lastBlock = (negStrand) ? 0 : query->blockCount-1 ;
        int qqe = query->qStarts[j] + query->blockSizes[j];
        int qqs = query->qStarts[j] ;
        int qte = query->tStarts[j] + query->blockSizes[j];
        int qts = query->tStarts[j] ;
        int qpts = pslInRange(query, j-offset) ? query->tStarts[j-offset] : query->tStarts[firstBlock];
        int qpte = pslInRange(query, j-offset) ? query->tStarts[j-offset]+query->blockSizes[j-offset] : query->tStarts[firstBlock]+query->blockSizes[lastBlock];
        int qsts = pslInRange(query, j+offset) ? query->tStarts[j+offset] : query->tStarts[lastBlock];
        int qste = pslInRange(query, j+offset) ? query->tStarts[j+offset]+query->blockSizes[j+offset] : query->tStarts[0]+query->blockSizes[lastBlock];
        int distLeftExon = 0;
        int distRightExon = 0;
        verbose(6, "ps j-offset %d before %d %d prev %d-%d succ %d-%d ",
                        j-offset, qts, qte, qpts, qpte, qsts, qste);
        if (negStrand)
            {
            reverseIntRange(&qts, &qte, query->tSize);
            reverseIntRange(&qpts, &qpte, query->tSize);
            reverseIntRange(&qsts, &qste, query->tSize);
            }
        verbose(6, "t after %d %d prev %d-%d succ %d-%d diff %d %d\n",qts, qte, qpts, qpte, qsts, qste,
                        abs(qpte-qts), abs(qte-qsts));
        distLeftExon =  qts-qpte;
        distLeftExon -= countRepeatOverlap(query->tName, qpte, qts);
        distRightExon = qsts-qte;
        distRightExon -= countRepeatOverlap(query->tName, qte, qsts);
        if (distLeftExon < 0)
                distLeftExon = 0;
        if (distRightExon < 0)
                distRightExon = 0;
        if (negStrand)
                {
                int temp = distLeftExon;
                distLeftExon = distRightExon;
                distRightExon = temp;
                }
        verbose(6, "DIST %d %d prev %d-%d succ %d-%d left %d rep %d %s:%d-%d rt %d rep %d %s:%d-%d\n    i=%d/%d j=%d/%d\n",
                        qts, qte, qpts, qpte, qsts, qste,
                        distLeftExon, countRepeatOverlap(query->tName, qpte, qts), 
                        query->tName, qpte, qts,
                        distRightExon, countRepeatOverlap(query->tName, qte, qsts), 
                        query->tName, qte, qsts,
                        i,target->blockCount,j,query->blockCount);
        if ((abs(qqs - qs) <= spliceDrift) && qqs != target->qStart && 
                //((negStrand) ? abs(qsts-qte) : abs(qpte-qts)) > spliceDrift &&
                (distLeftExon > intronSlop || (j == 0) /* first exon has no chance for intron*/) &&
                i > 0 /* skip beginning of gene*/) 
            {
            count++;
            verbose(4, "   COUNT LEFT dist=%d qdist %d-%d=%d ss %d %s %d %d %s %d %d\n", 
                        distLeftExon, qqs, qs, abs(qqs-qs), count, 
                        target->tName, ts, te, query->tName, qts, qte);
            }
        else
            verbose(6, "   no left dist=%d qdist %d-%d=%d ss %d i > 0 %d qqs %d != target->qStart %d\n", 
                        distLeftExon, qqs, qs, abs(qqs-qs), count, i, qqs, target->qStart);
        if ((abs(qqe - qe) <= spliceDrift) && qqe != target->qEnd   && 
                //((negStrand) ? abs(qts-qpte) : abs(qte-qsts)) > spliceDrift  &&
                (distRightExon > intronSlop || (j == target->blockCount-1))  &&
                i < target->blockCount -1 /* skip end of gene */)
            {
            count++;
            verbose(4, "   COUNT RIGHT dist=%d qdist %d-%d=%d ss %d\n", 
                        distRightExon, qqe, qe , abs(qqe-qe), count);
            verbose(4, "   COUNT RIGHT i=%d j=%d %c q %d-%d parent q %d-%d \
t %c %s:%d-%d q %s %s:%d-%d target start %d End %d left %d %d-%d right %d %d-%d\n",
                 i,j,query->strand[0],qqs,qqe, qs,qe, 
                 target->strand[0],target->tName, ts, te, 
                 query->qName, query->tName, qts,qte, target->qStart, target->qEnd,
                 qpte-qts, qpte, qts, qte-qsts, qte, qsts
                 );
            }
        else
            verbose(6, "   NO i=%d j=%d %c qqs/e %d-%d parent qs/e %d-%d \
t %c %s:%d-%d q %s:%d-%d \
qqe-qe %d qqe %d qEnd %d qte-qsts %d qte %d qsts %d\n",
                 i,j,query->strand[0], qqs, qqe, qs, qe,
                 target->strand[0], target->tName, ts, te, 
                 query->tName, qts,qte , 
                 abs(qqe - qe), qqe, target->qEnd  , abs(qte-qsts), qte, qsts);
        }
    }
verbose(2, "FINAL COUNT Cons Splices %d %s:%d-%d q %s:%d-%d\n",
     count,
     target->tName,  target->tStart+1, target->tEnd,
     query->tName, query->tStart+1, query->tEnd);
/* if target or query PSL coords were swapped then revrse them again */
if (swapT)
    pslRc(target);
if (swapQ)
    pslRc(query);
/* if second strand is '-' then set to null */
if (query->strand[1] == '+')
    query->strand[1] = '\0';
if (target->strand[1] == '+')
    target->strand[1] = '\0';
return count;
}

float calcTrfRatio(struct psl *psl, struct hash *trfHash)
/* calc ratio of simple repeats overlap to aligning bases */
{
int trf = 0;
if (trfHash != NULL)
    {
    int i;
    for (i = 0 ; i < psl->blockCount ; i++)
        {
        int ts = psl->tStarts[i] ;
        int te = psl->tStarts[i] + psl->blockSizes[i];
        struct binKeeper *bk = hashFindVal(trfHash, psl->tName);
        if (bk != NULL)
            {
            struct binElement *el, *elist = binKeeperFindSorted(bk, ts, te ) ;
            for (el = elist; el != NULL ; el = el->next)
                {
                trf += positiveRangeIntersection(ts, te, el->start, el->end);
                verbose(7,"%s trf overlap %d psl match+mismatch %d\n",
                        psl->qName,trf,(psl->match+psl->misMatch));
                }
            slFreeList(&elist);
            }
        }
    }
verbose(5,"%s trf overlap %d psl match+mismatch %d ratio %4.3f\n",
        psl->qName,trf,(psl->match+psl->misMatch), (float)trf/(float)(psl->match+psl->misMatch));
return (float)trf/(float)(psl->match+psl->misMatch);
}

int overlapMrna(struct psl *psl, int *exonOverlapCount, struct psl **overlapPsl, int *exonCount)
/* Count bases that mRNA (exprHash) overlaps with pseudogenes. If self match 
 * then don't filter it.*/
/* exonOverlapCount has number of exons in matched mRNA */
{
int maxOverlap = 0;
safef(mrnaOverlap,255,"NONE");
if (exprHash != NULL)
    {
    int mrnaBases = 0;
    struct psl *mPsl = NULL , *mPslMerge = NULL;
    struct binKeeper *bk = hashFindVal(exprHash, psl->tName);
    if (bk != NULL)
        {
        struct binElement *el, *elist = binKeeperFindSorted(bk, psl->tStart , psl->tEnd ) ;
        int blockIx;
        for (el = elist; el != NULL ; el = el->next)
            {
            mrnaBases = 0;
            mPsl = el->val;
            if (mPsl != NULL)
                {
                assert (psl != NULL);
                assert (mPsl != NULL);
                assert (psl->tName != NULL);
                if (differentString(psl->qName, mPsl->qName))
                    {
                    mPslMerge = mPsl;
                    if (exonCount != NULL)
                        *exonCount = mPsl->blockCount;
                    assert(mPslMerge != NULL);
                    verbose(6,"blk %d %s %d ec %d\n",mPslMerge->blockCount, mPsl->qName, mrnaBases, *exonOverlapCount);

                    if (mPslMerge->blockCount > 0) /* mPslMerge contains overlap of retro with all_mrna.psl. 
						no need to loop through intersection blocks, if only one */
                        {
                        for (blockIx = 0; blockIx < mPslMerge->blockCount; ++blockIx)
                            {
                            mrnaBases += positiveRangeIntersection(psl->tStart, psl->tEnd, 
                                mPslMerge->tStarts[blockIx], mPslMerge->tStarts[blockIx]+mPslMerge->blockSizes[blockIx]);
                            verbose(6,"%s overlapMrna %s %d %s %d-%d\n", mPslMerge->qName, psl->qName, mrnaBases, 
                                    mPslMerge->tName, mPslMerge->tStarts[blockIx], mPslMerge->tStarts[blockIx]+mPslMerge->blockSizes[blockIx]);
                            }
                        }
                    else
                        {
                        mrnaBases += positiveRangeIntersection(psl->tStart, psl->tEnd, mPsl->tStart, mPsl->tEnd);
    //                    verbose(6,"blk merge %d %s %d ec %d\n",mPslMerge->blockCount, mPsl->qName, mrnaBases, *exonOverlapCount);
                        }
		    /* mrnaBases now contains overlap between retro and all_mrna.psl */
                    verbose(6,"MRNABASES %d block cnt %d maxOverlap %d exonOverlapCount %d \
                            %s %s %d-%d best so far %s\n",
                            mrnaBases, mPslMerge->blockCount, maxOverlap, *exonOverlapCount, 
                            mPslMerge->qName, mPslMerge->qName, mPslMerge->tStart, mPslMerge->tEnd, mrnaOverlap);
		    /* must have at least 50 bases overlap to be consider 
                      * "expressed" retro */
		    /* pick mRNA with greatest overlap with retro */
 		    /* Refseq mRNA take priority over non-refseq mRNA */
                    if (mrnaBases > 50 && (mPslMerge->blockCount > 0) && 
				/* throw out small alignment blocks less than 50bp */
                            (((int)mPslMerge->blockCount > *exonOverlapCount) || 
                             (((int)mPslMerge->blockCount == *exonOverlapCount) && 
                              ((mrnaBases > maxOverlap) || (startsWith("NM",mPslMerge->qName) && !startsWith("NM",mrnaOverlap)) 
                               )))
                       )
                        {
			/* return # of                      overlapping exons, bases and qName 
                         * of best overlapping mRNA*/
                            *exonOverlapCount = (int)mPslMerge->blockCount;
                            safef(mrnaOverlap,255,"%s",mPslMerge->qName);
                            maxOverlap = mrnaBases;
                            *overlapPsl = mPslMerge;
                        }
                    }
                }
            }
        slFreeList(&elist);
        }
    }
return maxOverlap ;
}

void pseudoFeaturesCalc(struct psl *psl, struct psl *bestParentPsl, int maxExons, int bestAliCount, 
        char *bestChrom, int bestStart, int bestEnd, char *bestStrand) 
/* calculate features of retroGene */
{
struct ucscRetroInfo *pg = NULL;
struct dyString *iString = newDyString(16*1024);
struct dyString *reason = newDyString(255);
struct genePred *gp = NULL, *kg = NULL, *mgc = NULL;
int milliMinPseudo = round(1000*minAliPseudo);
int processedIntrons = 0;    
int intronCount = 0;    
int exonCover = 0;
int qBlockCover = 0;
//int conservedSpliceSites = 0;    
int geneOverlap = -1;
int polyAstart = 0;
int polyAend = 0;
//int tReps , qReps;
int trf = 0, rep = 0;
float trfRatio = 0;
bool keepChecking = TRUE;
int intronBases;
int pseudoExonCount = 0;

/* set and initalise the ucscRetroInfo struct */
AllocVar(pg);
pg->name = cloneString(psl->qName);
pg->parentSpliceCount = (maxExons*2)-2; /* really splice sites */
pg->bestAliCount = bestAliCount;
pg->alignGapCount = intronFactor(psl, rmskHash, trfHash, &intronBases);
pg->alignGapCount *= (intronBases/10);
pg->gChrom = cloneString(bestChrom);
pg->gStart = bestStart;
pg->gEnd = bestEnd;
safef(pg->gStrand, sizeof(pg->gStrand), bestStrand);
pg->milliBad = calcMilliScore(psl);
pg->coverage = ((psl->match+psl->misMatch+psl->repMatch)*100)/psl->qSize;
verbose(1,"\nchecking new %s:%d-%d %s best %s:%d-%d milli %d cover %d parent exons %d strand %c \n", 
         psl->tName, psl->tStart+1, psl->tEnd, psl->qName,  bestChrom, bestStart+1, bestEnd, 
         pg->milliBad, pg->coverage, pg->parentSpliceCount, psl->strand[0]);
pg->overStart = pg->overExonCover = pg->kStart = pg->kEnd = pg->rStart = pg->rEnd = pg->mStart = pg->mEnd = -1;
pg->polyA = polyACalc(psl->tStart, psl->tEnd, psl->strand, psl->tSize, psl->tName, 
                POLYAREGION, &polyAstart, &polyAend, pg->milliBad/10);
pg->polyAlen = abs(polyAend-polyAstart)+1;
pg->polyAstart = polyAstart;
/* calc introns processed in retro in two ways and blocks and exons covered into RETRO*/
calcIntrons(psl, maxBlockGap, bestParentPsl, &processedIntrons, &intronCount, &qBlockCover, &pseudoExonCount);
/* reverse the paramters to calculate number of exons covered in parent */
if (bestParentPsl != NULL)
    calcIntrons(bestParentPsl, maxBlockGap, psl, NULL, NULL, &exonCover, NULL);
/* if bestParentPsl is NULL and exonCover not calculated this will be set to 0 */
pg->exonCover = exonCover;
pg->retroExonCount = pseudoExonCount;
pg->processedIntrons = processedIntrons;
pg->intronCount = intronCount;
/* target is bestParentPsl, query is psl */
pg->conservedSpliceSites = countRetainedSpliceSites(bestParentPsl, psl , spliceDrift);
/* reduced intronsProcessed by conserved splice site count */
pg->processedIntrons -= pg->conservedSpliceSites;
if (pg->processedIntrons < 0)
    pg->processedIntrons = 0;
/* Calculate the ratio of simple repeat based to aligned non-repeat bases 
 * (matches + mismatches) */
trfRatio = calcTrfRatio(psl, trfHash);
verbose(4, "%s calcIntrons.ExonsSpliced_exon_covered - conserved_SS -> %d-%d=%d calcIntrons.intronCount %d trfRatio %4.2f\n", 
        psl->qName, pg->exonCover,pg->conservedSpliceSites, 
        pg->exonCover-pg->conservedSpliceSites, pg->intronCount, trfRatio);
//if (bestParentPsl == NULL)
//    pg->intronCount = pg->alignGapCount;

/* find overlapping gene annotations */
geneOverlap = 0;
genePredFree(&kg); 
genePredFree(&gp); 
genePredFree(&mgc);
/* if bestParentPsl (parent) is not NULL then find best overlapping genePred from
 * three genesets that are input to the program */ 
if (bestParentPsl != NULL)
    {
    kg = getOverlappingGene2(&kgList, "knownGene", bestParentPsl->tName, bestParentPsl->tStart, 
                        bestParentPsl->tEnd , bestParentPsl->qName, &geneOverlap);
    if (kg != NULL)
        {
        pg->kStart = kg->txStart;
        pg->kEnd = kg->txEnd;
        pg->kgName = cloneString(kg->name);
        }
    else
        {
        pg->kgName = cloneString("noKg");
        }
    gp = getOverlappingGene2(&gpList1, "refGene", bestParentPsl->tName, bestParentPsl->tStart, 
                        bestParentPsl->tEnd , bestParentPsl->qName, &geneOverlap);
    if (gp != NULL)
        {
        pg->refSeq = cloneString(gp->name);
        pg->rStart = gp->txStart;
        pg->rEnd = gp->txEnd;
        }
    else
        {
        pg->refSeq = cloneString("noRefSeq");
        }
    mgc = getOverlappingGene2(&gpList2, "mgcGenes", bestParentPsl->tName, bestParentPsl->tStart, 
                        bestParentPsl->tEnd , bestParentPsl->qName, &geneOverlap);
    if (mgc != NULL)
        {
        pg->mgc = cloneString(mgc->name);
        pg->mStart = mgc->txStart;
        pg->mEnd = mgc->txEnd;
        }
    else
        {
        pg->mgc = cloneString("noMgc");
        }
    }
else
    {
    pg->refSeq = cloneString("noRefSeq");
    pg->kgName = cloneString("noKg");
    pg->mgc = cloneString("noMgc");
    }

if (trfRatio > maxTrf)
    {
    verbose(1,"NO. %s trf overlap %f > %4.2f %s %d \n",
            psl->qName, (float)trf/(float)(psl->match+psl->misMatch) , maxTrf,
            psl->tName, psl->tStart);
    keepChecking = FALSE;
    }
/* blat sometimes overlaps parts of the same mrna , filter these */

/* if this is quite repetitive and there is an intersection between the parent
 * start and end and the retro start and end and they are on the same chrom */
if ( keepChecking && positiveRangeIntersection(bestStart, bestEnd, psl->tStart, psl->tEnd) && 
            sameString(psl->tName, bestChrom))
   {
   verbose(2,"NO. self overlap %s %d %d parent %s %d %d\n",
            psl->tName, psl->tStart, psl->tEnd,
            bestChrom, bestStart, bestEnd
            );
   dyStringAppend(reason,"self;");
   keepChecking = FALSE;
   }

/* count repeat overlap with pseudogenes and skip ones with more than maxRep% overlap*/
if (keepChecking && rmskHash != NULL)
    {
    int i;
    struct binElement *el, *elist;
    rep = 0;
    assert (psl != NULL);
    for (i = 0 ; i < psl->blockCount ; i++)
        {
        int ts = psl->tStarts[i] ;
        int te = psl->tStarts[i] + psl->blockSizes[i];
        struct binKeeper *bk;
        bk = hashFindVal(rmskHash, psl->tName);
        if (bk != NULL)
            {
            elist = binKeeperFindSorted(bk, ts, te) ;
            for (el = elist; el != NULL ; el = el->next)
                {
                rep += positiveRangeIntersection(ts, te, el->start, el->end);
                }
            slFreeList(&elist);
            }
        }

    }
pg->tReps = round((float)(rep*100)/(float)(psl->match+(psl->misMatch)));
/* label predictions with high repeat coverage as NOTPSEUDO */
if ((float)rep/(float)(psl->match+(psl->misMatch)) > maxRep )
    {
    verbose(1,"NO %s reps %.3f %.3f\n",psl->tName,(float)rep/(float)(psl->match+(psl->misMatch)) , maxRep);
    dyStringAppend(reason,"maxRep;");
    pg->label = NOTPSEUDO;
    keepChecking = FALSE;
    }

/* for retros with low number of introns and certain score and low repeats, 
 * check mRNAs aligned by BLAT for overlap and therefore signs of the retro
 * being expressed */
if (keepChecking && (pg->intronCount <= 2 /*|| (pg->exonCover - pg->intronCount > INTRONMAGIC)*/) && 
    /*maxExons > 1 && */ pg->bestAliCount > 0 && bestChrom != NULL &&
    (calcMilliScore(psl) >= milliMinPseudo && trfRatio < maxTrf && (pg->exonCover-pg->conservedSpliceSites) > 0 &&
    psl->match + psl->misMatch + psl->repMatch >= minCoverPseudo * (float)psl->qSize))
    {
    struct psl *mPsl = NULL;
    int exonOverlapCount = -1;
    int exonCount = -1;
    int maxOverlap = overlapMrna(psl, &exonOverlapCount, &mPsl, &exonCount);
    pg->maxOverlap = maxOverlap;
    if ((float)maxOverlap/(float)(psl->match+psl->misMatch+psl->repMatch) > splicedOverlapRatio 
            && maxOverlap > 10 ) 
        /* if overlap > 50 bases  and 10% overlap with pseudogene, then skip */
        {
        verbose(1,"NO %s:%d-%d %s expressed blat mrna %s %d bases overlap %f %%\n",
                psl->tName, psl->tStart+1, psl->tEnd, psl->qName,mrnaOverlap, 
                maxOverlap, (float)maxOverlap/(float)psl->qSize);
        dyStringAppend(reason,"expressed");
        pg->overName = cloneString(mPsl->qName); 
        pg->overStart = mPsl->tStart;
        pg->overExonCover = exonOverlapCount;
        strncpy(pg->overStrand, mPsl->strand , sizeof(pg->overStrand));
        //pslFree(&mPsl);
        pg->label = EXPRESSED;
        outputLink(psl, pg, reason, bestParentPsl);
        keepChecking = FALSE;
        }

    if (keepChecking)
       {
       verbose(2,"YES %s %d rr %3.1f rl %d ln %d %s iF %d maxE %d bestAli %d isp %d millibad %d match %d cover %3.1f rp %d polyA %d len %d start %d overlap ratio %d/%d %s\n",
            psl->qName,psl->tStart,((float)rep/(float)(psl->tEnd-psl->tStart) ),rep, 
            psl->tEnd-psl->tStart,psl->tName, pg->intronCount, 
            maxExons , pg->bestAliCount, pg->exonCover,
            calcMilliScore(psl),  psl->match + psl->misMatch + psl->repMatch , 
            minCoverPseudo * (float)psl->qSize, pg->tReps , pg->polyA, 
            pg->polyAlen, pg->polyAstart ,
            maxOverlap,psl->match+psl->misMatch+psl->repMatch , pg->overName
            );
       if ((pg->exonCover*2) - pg->conservedSpliceSites < 2)
           {
           verbose(4, "%s fake single exon ; exonCover %d - consSS %d < 2 \n",
                         psl->qName,  pg->exonCover , pg->conservedSpliceSites );

           dyStringAppend(reason,"singleExon;");
           pg->label = PSEUDO;
           outputLink(psl, pg, reason, bestParentPsl);
           }
       else if (bestParentPsl == NULL)
           {
           dyStringAppend(reason,"noBest;");
           pg->label = NOTPSEUDO;
           outputLink(psl, pg, reason, bestParentPsl);
           }
       else if (kg == NULL && mgc == NULL && gp == NULL)
           {
           dyStringAppend(reason,"mrna");
           pg->label = PSEUDO;
           outputLink(psl, pg, reason, bestParentPsl);
           }
       else
           {
           dyStringAppend(reason,"pseudogene");
           pg->label = PSEUDO;
           outputLink(psl, pg, reason, bestParentPsl);
           }
           keepChecking = FALSE;
       }
    }
else
    {
    if (bestParentPsl == NULL)
        dyStringAppend(reason,"noBest;");
    if (pg->exonCover < 1)
        dyStringAppend(reason,"exonCover;");
    if (pg->bestAliCount < 1)
        dyStringAppend(reason,"noAli;");
    if (trfRatio > maxTrf)
        dyStringAppend(reason,"trf;");
    if (pg->intronCount > 0)
        dyStringAppend(reason,"introns;");
    if (pg->exonCover == 1)
        dyStringAppend(reason,"singleExon;");
    if (maxExons <= 1)
        dyStringAppend(reason,"singleExonParent;");
    if (calcMilliScore(psl) < milliMinPseudo)
        dyStringAppend(reason,"milliBad;");
    if (psl->match + psl->misMatch + psl->repMatch < minCoverPseudo * (float)psl->qSize)
        dyStringAppend(reason,"coverage;");
   /* NOTE: This is doing the same thing in both if and else */
   if (bestParentPsl == NULL)
       {
       pg->label = NOTPSEUDO;
       outputLink(psl, pg, reason, bestParentPsl);
       }
    else
       {
       pg->label = NOTPSEUDO;
       outputLink(psl, pg, reason, bestParentPsl);
       }
    verbose(2,"NO. %s %s %d rr %3.1f rl %d ln %d %s iF %d maxE %d bestAli %d isp %d score %d match %d cover %3.1f rp %d\n",
        reason->string, psl->qName,psl->tStart,((float)rep/(float)(psl->tEnd-psl->tStart) ),rep, 
        psl->tEnd-psl->tStart,psl->tName, pg->intronCount, maxExons , pg->bestAliCount, pg->exonCover,
        calcMilliScore(psl),  psl->match + psl->misMatch + psl->repMatch , 
        minCoverPseudo * (float)psl->qSize, pg->tReps );
    }
dyStringFree(&iString);
dyStringFree(&reason);
}

void processBestMulti(char *acc, struct psl *pslList)
/* This function is passed a list of all hits of a single mrna */
/* Find psl's that are align best anywhere along their length. */

{
struct psl *bestParentPsl = NULL, *psl, *bestSEPsl = NULL;
int qSize = 0;
int *scoreTrack = NULL;
int maxExons = 0;
int milliScore;
int goodAliCount = 0;
int bestAliCount = 0;
int milliMin = 1000*minAli;
int bestStart = -1, bestEnd = -1;
static char bestStrand[3];
static char bestSEStrand[3];
int bestSEStart = 0, bestSEEnd = 0;
char *bestChrom = NULL, *bestSEChrom = NULL;
int bestScore = 0, bestSEScore = 0;

if (pslList == NULL)
    return;

/* calculate size of scoreArray by finding the longest aligment - some have polyA tail stripped off */
for (psl = pslList; psl != NULL; psl = psl->next)
    if (psl->qSize > qSize)
        qSize = psl->qSize;

AllocArray(scoreTrack, qSize+1);

/* search all alignments and store the best score for each base in scoreArray*/
for (psl = pslList; psl != NULL; psl = psl->next)
    {
    int blockIx;
    char strand = psl->strand[0];

    assert (psl!= NULL);
    milliScore = calcMilliScore(psl);
    verbose(2,"checking %s %s:%d-%d milliScore %d milliMin %d\n",psl->qName, psl->tName, psl->tStart+1, psl->tEnd, milliScore, milliMin);
    if (milliScore >= milliMin)
	{
	++goodAliCount;
	milliScore += sizeFactor(psl, rmskHash, trfHash);
        verbose(5,"@ %s %s:%d\n", psl->qName, psl->tName, psl->tStart+1);
	for (blockIx = 0; blockIx < psl->blockCount; ++blockIx)
	    {
	    int start = psl->qStarts[blockIx];
	    int size = psl->blockSizes[blockIx];
	    int end = start+size;
	    int i;
            /* if on - strand, reverse coords so always checking same sequence
               and therefore the same bases in the same order */
            if (strand == '-')
	        reverseIntRange(&start, &end, psl->qSize);
	    if (start < 0 || end > psl->qSize || psl->qSize > qSize)
		{
		warn("Error: qName %s tName %s qSize %d psl->qSize %d start (0 based) %d end %d",
		    psl->qName, psl->tName, qSize, psl->qSize, start, end);
		}
            verbose(5,"milliScore: %d qName: %s tName: %s:%d-%d qSize: %d psl->qSize %d i %d start %d end %d \n",
                    milliScore, psl->qName, psl->tName, psl->tStart+1, psl->tEnd, qSize, i, psl->qSize, start, end);
            /* picking highest scoring PSL for each base of the alignment */
	    for (i=start; i<end; ++i)
		{
                assert(i<=qSize);
		if (milliScore > scoreTrack[i])
                    {
                    verbose(6," %d base %d",milliScore, i);
                    scoreTrack[i] = milliScore;
                    }
		}
	    }
	}
    }
/* scoreTrack contains best alignment score per base of retro */
verbose(2,"---finding best---\n");
/* Print out any alignments that are within minTop% of top score for at least . */
bestScore = 0;
bestSEScore = 0;
for (psl = pslList; psl != NULL; psl = psl->next)
    {
    struct psl *pslMerge;
    int score = calcSizedScore(psl, rmskHash, trfHash);
    
    verbose(3,"milli %d > %d match %d > %d score: %d best %d qName %s tName %s:%d-%d \n",
            calcMilliScore(psl), milliMin, 
            (psl->match + psl->repMatch + psl->misMatch) , round(minCover * psl->qSize),
            score, bestScore, psl->qName, psl->tName, psl->tStart+1, psl->tEnd );
    if (
        calcMilliScore(psl) >= milliMin && closeToTop(psl, scoreTrack, score) 
        && (psl->match + psl->repMatch + psl->misMatch) >= round(minCover * psl->qSize))
	{
        ++bestAliCount;
        AllocVar(pslMerge);
        pslMergeBlocks(psl, pslMerge, 30);
        verbose(4,"merge blockCount %d \n", pslMerge->blockCount);
        assert (pslMerge->blockCount > 0);
        
        if (score  > bestScore && pslMerge->blockCount > 1)
            {
            bestParentPsl = psl;
            bestStart = psl->tStart;
            bestEnd = psl->tEnd;
            bestChrom = cloneString(psl->tName);
            bestScore = score;
            safef(bestStrand, sizeof(bestStrand), psl->strand );
            verbose(2,"BEST score: %d tName %s:%d \n",score,psl->tName,psl->tStart+1);
            }
        if (score  > bestSEScore )
            {
            bestSEPsl = psl;
            bestSEStart = psl->tStart;
            bestSEEnd = psl->tEnd;
            bestSEChrom = cloneString(psl->tName);
            bestSEScore = score  ;
            safef(bestSEStrand, sizeof(bestSEStrand), psl->strand );
            verbose(2,"BEST single Exon score: %d tName %s:%d \n",score,psl->tName,psl->tStart+1);
            }
        if (pslMerge->blockCount > maxExons )
            maxExons = pslMerge->blockCount;
        pslFree(&pslMerge);
	}
        else
        {
        verbose(3,"NOT BEST milli %d > %d match %d > %d score: %d best %d qName %s tName %s:%d-%d \n",
            calcMilliScore(psl), milliMin, 
            (psl->match + psl->repMatch + psl->misMatch) , round(minCover * psl->qSize),
            score, bestScore, psl->qName, psl->tName, psl->tStart+1, psl->tEnd );
        }
    }
if (bestScore== 0)
    { /* take best single exon hit, if no multi exon */
    bestParentPsl = bestSEPsl;
    bestStart = bestSEStart;
    bestEnd = bestSEEnd;
    bestChrom = bestSEChrom;
    bestScore = bestSEScore;
    safef(bestStrand , sizeof(bestStrand), bestSEStrand) ;
    }
if (bestChrom != NULL)
    verbose(2,"---DONE finding best--- %s:%d-%d\n",bestChrom, bestStart+1, bestEnd);
/* bestParentPsl contains "best" parent alignment for each retro */
/* output parent genes, retrogenes, and calculate feature vector */
/* this is the main loop */
for (psl = pslList; psl != NULL; psl = psl->next)
    {
    int score = calcSizedScore(psl, rmskHash, trfHash);

    verbose(2,"checking qName for extension %s \n",psl->qName, psl->tName, psl->tStart+1, psl->tEnd);
    if (strstr(psl->qName, "-") == NULL)
        {
        verbose(2,"skipping score step for blat alignment %s %s:%d-%d\n",psl->qName, psl->tName, psl->tStart+1, psl->tEnd);
        continue;
        }
    if (
        calcMilliScore(psl) >= milliMin && closeToTop(psl, scoreTrack, score) 
        && psl->match + psl->misMatch + psl->repMatch >= minCover * psl->qSize)
        {
        verbose(3,"\n##not checking %s:%d-%d %s best %s:%d-%d milliscore %d threshold %d aligning %d / qSize %d  closeToTop %d\n", 
             psl->tName, psl->tStart+1, psl->tEnd, psl->qName,  bestChrom, bestStart+1, bestEnd, 
             calcMilliScore(psl), milliMin, psl->match + psl->misMatch + psl->repMatch,psl->qSize ,closeToTop(psl, scoreTrack, score) );
        /* write out the PSL for these best alignments to a tab-separated file
           (arg9) */
        pslTabOut(psl, bestFile);
        }
    else 
        {
        /* calculate various features of pseudogene and 
           output feature records */
        pseudoFeaturesCalc(psl, bestParentPsl, maxExons, bestAliCount, 
            bestChrom, bestStart, bestEnd, bestStrand );
        }
    }
freeMem(scoreTrack);
}

struct hash *readPslQnameHash(char *pslFileName)
{
/* hash Qnames in a psl list */
struct hash *hash = newHash(0);
char *row[21] ;
struct lineFile *pf = lineFileOpen(pslFileName , TRUE);
while (lineFileNextRow(pf, row, ArraySize(row)))
    {
    struct psl *psl = pslLoad(row);
    hashAdd(hash, psl->qName, psl);
    }
lineFileClose(&pf);
return hash;
}

void addBlatAlignment( char *name, struct psl **pslList)
    /* append blat alignments to blastz mrna alignemnts to
       make sure parent has all exons that blastz sometimes misses */
{
char choppedName[256];
struct hashEl *el, *elist = NULL;
assert(name!=NULL);
/* Remove the suffix from the name (id) and result is in choppedName */
safef(choppedName, sizeof(choppedName), "%s",name);
chopSuffix(choppedName);
elist = hashLookup(qNameHash, choppedName);
if (elist == NULL)
    {
    safef(choppedName, sizeof(choppedName), "%s",name);
    elist = hashLookup(qNameHash, choppedName);
    }
if (elist == NULL)
    {
    verbose(3, "chopping %s at - \n", name);
    safef(choppedName, sizeof(choppedName), "%s",name);
    chopSuffixAt(choppedName,'-');
    elist = hashLookup(qNameHash, choppedName);
    }
if (elist == NULL)
    {
    verbose(1, "Warning: blat alignment for mrna %s or %s not found in all_mrna.\n", choppedName,name);
    }
else
    for (el = elist; el != NULL ; el = el->next)
        {
        struct psl *aPsl = NULL;
        if (el == NULL)
            break;
        aPsl = el->val;
        if (sameString(el->name, choppedName))
            {
            verbose(3, "adding blat alignment %s:%d-%d %s\n",aPsl->tName, aPsl->tStart+1, aPsl->tEnd, aPsl->qName);
            slAddHead(pslList, aPsl);
            }
        else
            {
            //verbose(4, "scanning blat alignment %s %s\n",choppedName, aPsl->qName);
            break;
            }
        }
    //slFreeList(&elist);
}

void pslPseudo(char *inName, char *bestAliName, char *pseudoFileName, char *linkFileName, char *axtFileName, char *orthoFileName)
/* find best alignments with and without introns.
 * Put pseudogenes  in pseudoFileName. 
 * store link between pseudogene and gene in LinkFileName */
{
struct lineFile *in = pslFileOpen(inName);
int lineSize;
char *line;
char *words[32];
char sep = '.';
int wordCount;
struct psl *pslList = NULL, *psl;
char lastName[256] = "nofile";
int aliCount = 0;
if (bestAliName != NULL)
    bestFile = mustOpen(bestAliName, "w");
else 
    errAbort("bestAliName cannot be blank\n");
if (pseudoFileName != NULL)
    pseudoFile = mustOpen(pseudoFileName, "w");
if (linkFileName != NULL)
    linkFile = mustOpen(linkFileName, "w");
if (axtFileName != NULL)
    axtFile = mustOpen(axtFileName, "w");
if (orthoFileName != NULL)
    orthoFile = mustOpen(orthoFileName, "w");
else 
    errAbort("orthoFileName cannot be blank\n");

verbose(1,"Processing %s to %s and %s\n", inName, bestAliName, pseudoFileName);
 if (!noHead)
     pslWriteHead(bestFile);
safef(lastName, sizeof(lastName),"nofile");
/* read input split PSL file line by line */
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
    if (sameString(lastName, "nofile"))
        safef(lastName, sizeof(lastName), "%s", psl->qName);
    if (stripVersion)
        chopSuffix(psl->qName);
    verbose(2,"scoring %s version %d lastName %s\n",psl->qName, stripVersion, lastName);
    /* if new qName encountered, not the same as the last one read */
    if (!samePrefix(lastName, psl->qName, sep))
	{
        if (!skipBlatMerge)
            addBlatAlignment(lastName, &pslList);
        else
            verbose(2,"Skipping merge of blat alignments ");
        slReverse(&pslList);
	processBestMulti(lastName, pslList);
	//pslFreeList(&pslList);
        pslList = 0x0;
	safef(lastName, sizeof(lastName), "%s", psl->qName);
	}
    slAddHead(&pslList, psl);
    }
/* why is this here */
addBlatAlignment(lastName, &pslList);
//slReverse(&pslList);
pslList = 0x0;
processBestMulti(lastName, pslList);
pslFreeList(&pslList);
lineFileClose(&in);
fclose(bestFile);
fclose(pseudoFile);
fclose(linkFile);
fclose(axtFile);
fclose(orthoFile);
verbose(1,"Processed %d alignments\n", aliCount);
}

int main(int argc, char *argv[])
/* Process command line. */
{
char *cdsFile = NULL;

optionInit(&argc, argv, optionSpecs);
if (argc != 20)
    usage(argc);
verboseSetLogFile("stdout");
verbosity = optionInt("verbose", 1);
verboseSetLevel(verbosity);
verbose(1,"version is %s\n",rcsid);
ss = axtScoreSchemeDefault();
/* smaller gap to handle spliced introns */
ss->gapExtend = 5;

// database name
db = cloneString(argv[1]);
minAli = optionFloat("minAli", minAli);
maxRep = optionFloat("maxRep", maxRep);
maxTrf = optionFloat("maxTrf", maxTrf);
minAliPseudo = optionFloat("minAliPseudo", minAliPseudo);
nearTop = optionFloat("nearTop", nearTop);
splicedOverlapRatio = optionFloat("splicedOverlapRatio", splicedOverlapRatio);
minCover = optionFloat("minCover", minCover);
minCoverPseudo = optionFloat("minCoverPseudo", minCoverPseudo);
minNearTopSize = optionInt("minNearTopSize", minNearTopSize);
minIntronSize = optionInt("minIntronSize", minIntronSize);
intronSlop = optionInt("intronSlop", intronSlop);
spliceDrift = optionInt("spliceDrift", spliceDrift);
maxBlockGap = optionInt("maxBlockGap" , maxBlockGap) ;
ignoreSize = optionExists("ignoreSize");
skipBlatMerge = optionExists("skipBlatMerge");
skipExciseRepeats = optionExists("skipExciseRepeats");
noIntrons = optionExists("noIntrons");
stripVersion = optionExists("stripVersion");
showall = optionExists("showall");
noHead = optionExists("nohead");
cdsFile = optionVal("cdsFile", NULL);
initWeights();

srand(time(NULL));
//sleep((float)rand()*10/RAND_MAX);
verbose(1,"Scanning %s\n", argv[13]);
//hashFileList(argv[13]);
verbose(1,"Loading mrna sequences from %s\n",argv[14]);
//mrnaList = faReadAllMixed(argv[14]);
/* Open and store information on genome sequence and mRNA sequences */
genomeSeqFile = twoBitOpen(argv[13]);
mrnaFile = twoBitOpen(argv[14]);
/* get list of ids from mrna 2bit file */
mrnaList = twoBitSeqNames(argv[14]);
verbose(1,"Loading genes from %s\n",argv[15]);
/* store genePreds from GENE2 in a list */
gpList1 = genePredLoadAll(argv[15]);
verbose(1,"Loading genes from %s\n",argv[16]);
/* store genePreds from GENE3 in a list */
gpList2 = genePredLoadAll(argv[16]);
verbose(1,"Loading genes from %s\n",argv[17]);
/* store genePreds from GENE1 (knownGene or equivalent) in a list */
kgList = genePredLoadAll(argv[17]);
verbose(1,"Loading net %s\n",argv[5]);
splitPath(argv[5], NULL, orthoNet1, NULL);
synHash = readNetToBinKeeper(argv[3], argv[5]);
verbose(1,"Loading net %s\n",argv[6]);
splitPath(argv[6], NULL, orthoNet2, NULL);
syn2Hash = readNetToBinKeeper(argv[3], argv[6]);
verbose(1,"Loading net %s\n",argv[18]);
splitPath(argv[18], NULL, orthoNet3, NULL);
syn3Hash = readNetToBinKeeper(argv[3], argv[18]);

verbose(1,"Loading Trf Bed %s\n",argv[7]);
trfHash = readBedCoordToBinKeeper(argv[3], argv[7], BEDCOUNT);

verbose(1,"Reading Repeats from %s\n",argv[4]);
rmskHash = readBedCoordToBinKeeper(argv[3], argv[4], BEDCOUNT);

verbose(1,"Reading mrnas from %s\n",argv[8]);
exprHash = readPslToBinKeeper(argv[3], argv[8]);
qNameHash = readPslQnameHash(argv[8]);
if (cdsFile != NULL)
    {
    verbose(1,"Reading cds from %s\n",cdsFile);
    loadCdsFile(cdsFile);
    }
verbose(1,"Scoring alignments from %s.\n",argv[2]);

pslPseudo(argv[2], argv[9], argv[10], argv[11], argv[12], argv[19]);

verbose(1,"freeing everything\n");
binKeeperPslHashFree(&exprHash);
binKeeperHashFree(&synHash);
binKeeperHashFree(&syn2Hash);
binKeeperHashFree(&syn3Hash);
binKeeperHashFree(&trfHash);
genePredFreeList(&gpList1);
genePredFreeList(&gpList2);
genePredFreeList(&kgList);
//freeDnaSeqList(&mrnaList);
twoBitClose(&genomeSeqFile);
twoBitClose(&mrnaFile);

if (abortAtEnd)
    {
    errAbort("mrna mismatches, pipeline aborting\n");
    return -4;
    }
else
    return 0;
}
