/* countFrameShift - using PSL alignments of mRNAs, count frame shifts of retrogene relative to parent. 
 */
#include "common.h"
#include "options.h"
#include "jksql.h"
#include "genbank.h"
#include "genePred.h"
#include "psl.h"
#include "hash.h"
#include "linefile.h"
#include "verbose.h"

static char const rcsid[] = "$Id: countFrameShift.c,v 1.3 2006/09/06 16:59:45 baertsch Exp $";

struct frameStats {
    struct frameStats *next;
    char *chrom;
    int start;
    int end;
    char strand[3];
    char *name;
    int score;
    int frame[3]; /* count of bases in each frame */
    int totalBases;
    int codingBases;
    int multiple3indel; /* count of multiple of 3 indels */
    int multiple3bases; /* count of bases in multiple of 3 indels */
    int clipStart;      /* number of bases retroGene start shorter than parent */
    int clipEnd;        /* number of bases retroGene stop shorter than parent */
};
/* command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"db", OPTION_STRING},
    {"cdsDb", OPTION_STRING},
    {"cdsFile", OPTION_STRING},
    {"requireUtr", OPTION_BOOLEAN},
    {"smallInsertSize", OPTION_INT},
    {"insertMergeSize", OPTION_INT},
    {"smallInsertSize", OPTION_INT},
    {"cdsMergeSize", OPTION_INT},
    {"cdsMergeMod3", OPTION_BOOLEAN},
    {"utrMergeSize", OPTION_INT},
    {"genePredExt", OPTION_BOOLEAN},
    {"allCds", OPTION_BOOLEAN},
    {"keepInvalid", OPTION_BOOLEAN},
    {"quiet", OPTION_BOOLEAN},
    {"verbose", OPTION_INT},
    {NULL, 0}
};

/* command line options */
static int gCdsMergeSize = -1;
static int gUtrMergeSize = -1;
static unsigned gPslOptions = genePredPslDefaults;
static boolean gRequireUtr = FALSE;
static boolean gKeepInvalid = FALSE;
static boolean gAllCds = FALSE;
static boolean gGenePredExt = FALSE;
static boolean gQuiet = FALSE;
struct genePred *gpList1 = NULL;

/* hash table of accession to CDS */
static struct hash *gCdsTable = NULL;

void usage()
/* Explain usage and exit. */
{
errAbort(
  "countFrameShift - count frame shifts in PSL alignments of mRNAs to retrogene\n"
  "usage:\n"
  "   countFrameShift [options] psl gene.gp outputFile\n"
  "\n"
  "Count frames shifts multi of 3 indels  in retrogenes based on alignments \n"
  "with parent mRNA. CDS annotation of mRNA comes from genbank. \n"
  "Accessions without valids CDS are optionally dropped. A best attempt is\n"
  "made to convert incomplete CDS annotations. CDS is clipped using overlapping\n"
  "gene annotation if it is smaller than the mapped retrogene annotation.\n"
  "\n"
  "The psl argument may either be a PSL file or a table in a databases,\n"
  "depending on options.  CDS maybe obtained from the database or file.\n"
  "Accession in PSL files are tried with and with out genbank versions.\n"
  "gene.gp can be any gene annotation tracck with cds.\n"
  "\n"
  "Options:\n"
  "  -db=db - get PSLs and CDS from this database, psl specifies the table.\n"
  "  -cdsDb=db - get CDS from this database, psl is a file.\n"
  "  -cdsFile=file - get CDS from this database, psl is a file.\n"
  "   File is table seperate with accession as the first column and\n"
  "   CDS the second\n"
  "  -allCds - consider PSL to be all CDS.\n"
  "  -keepInvalid - Keep sequences with invalid CDS.\n"
  "  -quiet - Don't print print info about dropped sequences.\n"
  "  -verbose=N - turn on debugging messages. 1 - default.\n"
  "\n", genePredStdInsertMergeSize, genePredStdInsertMergeSize);
}

void loadCdsFile(char *cdsFile)
/* read a CDS file into the global hash */
{
struct lineFile *lf = lineFileOpen(cdsFile, TRUE);
char *row[2];

gCdsTable = hashNew(20);
while (lineFileNextRowTab(lf, row, 2))
    hashAdd(gCdsTable, row[0], 
            lmCloneString(gCdsTable->lm, row[1]));

lineFileClose(&lf);
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

char *cdsQuery(struct sqlConnection *conn, char *acc, char *cdsBuf, int cdsBufSize)
/* query for a CDS, either in the hash table or database */
{
if (gCdsTable != NULL)
    return hashFindVal(gCdsTable, acc);
else
    {
    char query[512];
    safef(query, sizeof(query),
          "SELECT cds.name FROM cds,gbCdnaInfo WHERE (gbCdnaInfo.acc = '%s') AND (gbCdnaInfo.cds !=0) AND (gbCdnaInfo.cds = cds.id)",
          acc);
    return sqlQuickQuery(conn, query, cdsBuf, cdsBufSize);
    }
}

char *getCdsForAcc(struct sqlConnection *conn, char *acc, char *cdsBuf, int cdsBufSize)
/* look up a cds, trying with and without version */
{
char *cdsStr = cdsQuery(conn, acc, cdsBuf, cdsBufSize);
if (cdsStr == NULL)
    {
    int dotIdx = getGbVersionIdx(acc);
    if (dotIdx >= 0)
        {
        acc[dotIdx] = '\0';
        cdsStr = cdsQuery(conn, acc, cdsBuf, cdsBufSize);
        acc[dotIdx] = '.';
        }
    }
return cdsStr;
}

struct genbankCds getCds(struct sqlConnection *conn, struct psl *psl)
/* Lookup the CDS, either in the database or hash, or generate for query.  If
 * not found and looks like a it has a genbank version, try without the
 * version.  If allCds is true, generate a cds that covers the query.  Conn
 * maybe null if gCdsTable exists or gAllCds is true.  If CDS can't be
 * obtained, start and end are both set to -1.  If there is an error parsing
 * it, start and end are both set to 0. */
{
struct genbankCds cds;
ZeroVar(&cds);
if (gAllCds)
    {
    cds.start = psl->qStart;
    cds.end = psl->qEnd;
    if (psl->strand[0] == '-')
        reverseIntRange(&cds.start, &cds.end, psl->qSize);
    cds.startComplete = TRUE;
    cds.endComplete = TRUE;
    }
else
    {
    char cdsBuf[4096];
    char *cdsStr = getCdsForAcc(conn, psl->qName, cdsBuf, sizeof(cdsBuf));
    if (cdsStr == NULL)
        {
        if (!gQuiet)
            fprintf(stderr, "Warning: no CDS for %s\n", psl->qName);
        cds.start = cds.end = -1;
        }
    else
        {
        if (!genbankCdsParse(cdsStr, &cds))
            {
            if (!gQuiet)
                fprintf(stderr, "Warning: invalid CDS for %s: %s\n",
                        psl->qName, cdsStr);
            }
        else if ((cds.end-cds.start) > psl->qSize)
            {
            if (!gQuiet)
                fprintf(stderr, "Warning: CDS for %s (%u..%u) longer than qSize (%u)\n",
                        psl->qName, cds.start, cds.end, psl->qSize);
            cds.start = cds.end = -1;
            }
        }
    }
return cds;
}

static int getFrame(struct psl *psl, int start, int end,
                    struct genbankCds* cds)
/* get the starting frame for an exon of a mRNA.  start and end are
 * the bounds of the exon within the mRNA.  This may cover multiple psl
 * blocks due to merging of small gaps. */
{
/* use the 3' end is used if it's complete, as it is more often accurate when
 * genes are defined from mRNAs sequenced with reverse-transcriptase. */
int frame = -1;
/* map to mRNA coords in CDS since frame for an exon is in direction of
 * transcription. */
if (psl->strand[0] == '-')
    reverseIntRange(&start, &end, psl->qSize);
if (start < cds->start)
    start = cds->start;
if (end > cds->end)
    end = cds->end;

if (start < end)
    {
    /* compute from end if it is complete in mRNA */
    if (cds->endComplete)
        {
        int fr = (cds->end-start) % 3;
        frame = (fr == 2) ? 1 : ((fr == 1) ? 2 : 0);
        }
    else
        frame = (start-cds->start) % 3;
    }
else
    {
    verbose(3, "not in frame start %d>= end %d\n",start, end);
    }
return frame;
}

void calcFrame ( int *currentFrame, int frameShift)
/* calculate new frame after frameShift */
{
    verbose(3, "before %d + %d ",*currentFrame, frameShift);
*currentFrame += frameShift;
*currentFrame %= 3;
if (*currentFrame < 0)
    *currentFrame += 3;
if (*currentFrame < 0)
    errAbort("frameshift %d currentFrame %d\n",frameShift, *currentFrame);
verbose(3,"= %d\n",*currentFrame);
assert(*currentFrame >= 0);
assert(*currentFrame < 3);
}

int maxFrameRatio(struct frameStats *f)
/* score (0-1000) cases high if they have
   lots of frame shifts by scoring highly 
   if they have an even spread among frames 
   and score low, if all bases are in one frame. */
{
int score = 0;
if (f->codingBases == 0)
    return 0;

score = max(f->frame[0],f->frame[1]);
score = max(score,f->frame[2]);
score = 100-(score*100/f->codingBases);
assert(score <= 100);
//assert(score >=0);
return (score *10);
}
struct frameStats *calcFrameStats(struct psl *psl, struct genbankCds *cds)
    /* count frame shifts in psl */
{
int iExon = 0, startIdx = 0, stopIdx = 0, idxIncr = 0, iBlk = 0;
int prevTend = 0 , prevQend = 0, tBlkSize = 0;
int currentFrame = 0 ;
int qTot = 0, tTot = 0;
int geneOverlap = -1;
struct frameStats *frameStats = NULL;
struct genePred *gp = getOverlappingGene(&gpList1, "refGene", 
        psl->tName, psl->tStart, psl->tEnd , psl->qName, &geneOverlap);

if ((cds->start == -1 && cds->end == -1) || gp == NULL)
    return NULL;
verbose(3, "%s resulting gene(s) is %s count %d cds %d-%d\n",psl->qName, 
        gp->name, slCount(gp), gp->cdsStart, gp->cdsEnd);
verbose(3,"%s parent cds %d %d\n",
        psl->qName, cds->start,cds->end);
AllocVar(frameStats);

frameStats->chrom = psl->tName;
frameStats->start = psl->tStart;
frameStats->end = psl->tEnd;
frameStats->strand[0] = psl->strand[0];
frameStats->strand[1] = psl->strand[1];
frameStats->strand[2] = '\0';
frameStats->name = psl->qName;
frameStats->frame[0] = 0;
frameStats->frame[1] = 0;
frameStats->frame[2] = 0;
frameStats->totalBases = 0;
frameStats->codingBases = 0;

/* traverse psl in postive target order */
if (psl->strand[1] == '-')
    {
    startIdx = psl->blockCount-1;
    stopIdx = -1;
    idxIncr = -1;
    }
else
    {
    startIdx = 0;
    stopIdx = psl->blockCount;
    idxIncr = 1;
    }

if (psl->strand[0] == '-')
    {
    prevQend = psl->qEnd;
    }
prevTend = psl->tStarts[startIdx];
iExon = -1;  /* indicate none have been added */
for (iBlk = startIdx; iBlk != stopIdx; iBlk += idxIncr)
    {
    int tStart = psl->tStarts[iBlk];
    int tEnd = tStart + psl->blockSizes[iBlk];
    int qStart = psl->qStarts[iBlk];
    int qEnd = qStart + psl->blockSizes[iBlk];
    int tGap = tStart - prevTend;
    int qGap = qStart - prevQend;
    int fr;
    int frameShift = (tGap - qGap) % 3 ;
    int cdsStartOffset = 0;
    int cdsEndOffset = 0;
    int retroStartOffset = 0;
    int retroEndOffset = 0;
    if (psl->strand[1] == '-')
        reverseIntRange(&tStart, &tEnd, psl->tSize);
    if (psl->strand[0] == '-')
        reverseIntRange(&qStart, &qEnd, psl->qSize);
    qTot += (qEnd - qStart);
    tTot += (tEnd - tStart);
    verbose(3, "%s tTot %d s-e %d %d %s \n",psl->qName, tTot, tStart, tEnd, psl->strand);
    verbose(3, "%s tGap %d tStart %d prevTend %d\n",psl->qName, tGap, tStart, prevTend);
    if (psl->strand[0] == '-')
        {
        qGap = qEnd - prevQend - 1;
        frameShift = (tGap + qGap) ;
        verbose(3,"%s neg Gap %d end %d prev %d\n",
                psl->qName, qGap , qEnd , prevQend);
        }

    /* tBlkSize not used */
    if (frameShift == 0)
        tBlkSize += tEnd-tStart;
    else
        tBlkSize = tEnd-tStart;

    /* clip interval to the mapped parent coding start and end */
    if (cds->start > qStart && cds->start <= qEnd)
        cdsStartOffset = cds->start - qStart;
    if (cds->end < qEnd && cds->end >= qStart)
        cdsEndOffset = qEnd - cds->end;
    /* clip interval to the retro coding start and end */
    if (gp->cdsStart > tStart && gp->cdsEnd >= tStart)
        retroStartOffset = gp->cdsStart - tStart;
    if (gp->cdsEnd < tEnd && gp->cdsEnd >= tStart)
        retroEndOffset = tEnd - gp->cdsEnd;
    verbose(3, "%d %d compared with gene %d %d offsets %d %d\n",
            tStart, tEnd, gp->cdsStart, gp->cdsEnd,
            retroStartOffset, retroEndOffset);
    /* if the cds annotations do not match, take the smaller of the two */
    if (retroStartOffset > cdsStartOffset)
        {
        verbose(3, "%s *** retro cds override from %d to %d cds %d %d\n",
                psl->qName, cdsStartOffset, retroStartOffset, gp->cdsStart, gp->cdsEnd);
        cdsStartOffset = retroStartOffset;
        frameStats->clipStart += retroStartOffset;
        }
    if (retroEndOffset > cdsEndOffset)
        {
        verbose(3, "%s *** retro cds override from %d to %d cds %d %d\n",
                psl->qName, cdsEndOffset, retroEndOffset, gp->cdsStart, gp->cdsEnd);
        cdsEndOffset = retroEndOffset;
        frameStats->clipEnd += retroEndOffset;
        }
    /* this is used just to determine whether we are in coding or utr */
    fr = getFrame(psl, psl->qStarts[iBlk], psl->qStarts[iBlk]+psl->blockSizes[iBlk], cds);
    /* if coding not utr */
    if (fr >= 0 )
        {
        int addBases = (tEnd-tStart)-cdsStartOffset-cdsEndOffset;
        if (addBases < 0)
            {
            addBases = 0;
            frameShift = 0;
            }
        calcFrame(&currentFrame, frameShift);
        if ((tGap > 0 || qGap > 0) && frameShift ==0 && addBases > 0)
            {
            frameStats->multiple3indel++;
            frameStats->multiple3bases+= abs(tGap-qGap);
            verbose(3,"add mod3 cnt %d bases %d frshift %d\n",
                frameStats->multiple3indel, abs(tGap-qGap), frameShift);
            }
        frameStats->frame[currentFrame] += addBases;
        frameStats->codingBases += addBases;;
        verbose(3,"%s exon %d Tgap %d shift %d tBlksz %d prev %d t %d %d %s \n",
                psl->qName, iBlk,tGap, frameShift, tBlkSize, prevTend, tStart, tEnd,psl->tName );
        verbose(3,"%s exon %d Qgap %d frame %d  qsize %d prev %d q %d %d %c \n",
                psl->qName, iBlk,qGap, currentFrame, qEnd-qStart,prevQend, qStart, qEnd, psl->strand[0] );
        verbose(3,"%s FRAME %d ADD %d\n",
                psl->qName, currentFrame, addBases);
        }
    
    verbose(3,"%s coding? = %d ------------------------------------------<<<\n",psl->qName, fr);
    prevTend = tEnd;
    if (psl->strand[0] == '-')
        prevQend = qStart - 1;
    else
        prevQend = qEnd;
    }
frameStats->totalBases = tTot;
frameStats->score = maxFrameRatio(frameStats);
genePredFree(&gp); 
return (frameStats);
}

void printStats(struct frameStats *frameStats, FILE *outfile)
{
verbose(2,"%s frame 0/1/2 %d/%d/%d frameTot %d \n",
        frameStats->name, frameStats->frame[0], frameStats->frame[1], frameStats->frame[2],
        frameStats->frame[0]+frameStats->frame[1]+frameStats->frame[2]);
verbose(2,"%s frame %4.1f%% %4.1f%% %4.1f%% cdsTot %d total %d\n",
        frameStats->name, (float)frameStats->frame[0]/frameStats->codingBases*100, 
        (float)frameStats->frame[1]/frameStats->codingBases*100, 
        (float)frameStats->frame[2]/frameStats->codingBases*100, frameStats->codingBases, frameStats->totalBases);
verbose(2,"%s multi of three indels: %d bases %d\n",frameStats->name, frameStats->multiple3indel, frameStats->multiple3bases);
if( frameStats->frame[0] > 0 || 
    frameStats->frame[1] > 0 ||
    frameStats->frame[2] > 0) 
fprintf(outfile, "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%4.1f\t%4.1f\t%4.1f\t%d\t%d\n",
        frameStats->chrom, frameStats->start, frameStats->end,
        frameStats->name, frameStats->score, frameStats->strand,
        frameStats->frame[0], frameStats->frame[1], frameStats->frame[2], 
        frameStats->codingBases, frameStats->totalBases,
        (float)frameStats->frame[0]/frameStats->codingBases*100, 
        (float)frameStats->frame[1]/frameStats->codingBases*100, 
        (float)frameStats->frame[2]/frameStats->codingBases*100,
        frameStats->clipStart, frameStats->clipEnd);
}
void convertPslTableRow(char **row, FILE *outputFh)
/* A record from the PSL file  that includes CDS */
{
struct psl *psl = pslLoad(row+1);
struct  genbankCds cds;
struct frameStats *frameStats = NULL;
genbankCdsParse(row[0], &cds);
frameStats = calcFrameStats(psl, &cds);
if (frameStats != NULL)
    printStats(frameStats, outputFh);
pslFree(&psl);
}

void convertPslTable(struct sqlConnection *conn, char *pslTable, FILE *outputFh)
/* A row from the PSL query that includes CDS */
{
char query[512], **row;
struct sqlResult *sr;

/* generate join of cds with psls */
safef(query, sizeof(query),
      "SELECT cds.name,matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts "
      "FROM cds,%s,gbCdnaInfo WHERE (%s.qName = gbCdnaInfo.acc) AND (gbCdnaInfo.cds !=0) AND (gbCdnaInfo.cds = cds.id)",
      pslTable, pslTable);
sr = sqlGetResult(conn, query);
while ((row = sqlNextRow(sr)) != NULL)
    convertPslTableRow(row, outputFh);
sqlFreeResult(&sr);
}

void convertPslFileRow(struct sqlConnection *conn, char **row, FILE *outputFh)
/* A row from the PSL file, getting CDS */
{
struct psl *psl = pslLoad(row);
struct  genbankCds cds = getCds(conn, psl);
struct frameStats *frameStats = calcFrameStats(psl, &cds);
if (frameStats != NULL)
    printStats(frameStats, outputFh);
pslFree(&psl);
}

void convertPslFile(struct sqlConnection *conn, char *pslFile, FILE *outputFh)
{
struct lineFile *lf = pslFileOpen(pslFile);
char *row[PSL_NUM_COLS];

while (lineFileNextRowTab(lf, row, PSL_NUM_COLS))
    convertPslFileRow(conn, row, outputFh);
lineFileClose(&lf);
}

void countFrameShift(char *db, char *cdsDb, char *cdsFile, char *pslSpec,
                char *geneFile, char *outputFile)
/* countFrameShift - using PSL alignments of mRNAs, count frame shifts of retrogene relative to parent. */
{
struct sqlConnection *conn = NULL;
FILE* outputFh;

if (db != NULL)
    conn = sqlConnect(db);
else if (cdsDb != NULL)
    conn = sqlConnect(cdsDb);
outputFh = mustOpen(outputFile, "w");

if (cdsFile != NULL)
    loadCdsFile(cdsFile);
if (geneFile != NULL)
    gpList1 = genePredLoadAll(geneFile);

if (db == NULL)
    convertPslFile(conn, pslSpec, outputFh);
else
    convertPslTable(conn, pslSpec, outputFh);

if (ferror(outputFh))
    errAbort("error writing %s", outputFile);
carefulClose(&outputFh);
sqlDisconnect(&conn);
}

int main(int argc, char *argv[])
/* Process command line. */
{
char *db, *cdsDb, *cdsFile, *pslSpec, *geneFile, *outputFile;
int optCnt;
int verbosity = 1;

optionInit(&argc, argv, optionSpecs);
if (argc != 4)
    usage();
pslSpec = argv[1];
geneFile = argv[2];
outputFile = argv[3];
db = optionVal("db", NULL);
cdsDb = optionVal("cdsDb", NULL);
cdsFile = optionVal("cdsFile", NULL);
gRequireUtr = optionExists("requireUtr");
if (optionExists("cdsMergeMod3") && !optionExists("cdsMergeSize"))
    errAbort("must specify -cdsMergeSize with -cdsMergeMod3");
if (optionExists("cdsMergeSize") || optionExists("utrMergeSize"))
    {
    gCdsMergeSize = optionInt("cdsMergeSize", -1);
    gUtrMergeSize = optionInt("utrMergeSize", -1);
    if (optionExists("cdsMergeMod3"))
        gPslOptions |= genePredPslCdsMod3;
    if (optionExists("smallInsertSize") || optionExists("insertMergeSize"))
        errAbort("can't specify -smallInsertSize or -insertMergeSize with -cdsMergeSize or -utrMergeSize");
    }
else
    {
    int insertMergeSize = genePredStdInsertMergeSize;
    if (optionExists("smallInsertSize"))
        insertMergeSize = optionInt("smallInsertSize", genePredStdInsertMergeSize);
    insertMergeSize = optionInt("insertMergeSize", genePredStdInsertMergeSize);
    gCdsMergeSize = gUtrMergeSize = insertMergeSize;
    }
gGenePredExt = optionExists("genePredExt");
gKeepInvalid = optionExists("keepInvalid");
gAllCds = optionExists("allCds");
gQuiet = optionExists("quiet");
verbosity = optionInt("verbose", verbosity);
verboseSetLevel(verbosity);

if (gAllCds && ((cdsDb != NULL) || (cdsFile != NULL)))
    errAbort("can't specify -allCds with -cdsDb or -cdsFile");
if (gAllCds && gRequireUtr)
    errAbort("can't specify -allCds with -requireUtr");
/* this is a bit of work to implement */
if (gAllCds && (db != NULL))
    errAbort("can't specify -allCds with -db");

optCnt = 0;
if (db != NULL)
    optCnt++;
if (cdsDb == NULL)
    optCnt++;
if (cdsFile != NULL)
    optCnt++;
if (gAllCds)
    optCnt++;

if (optCnt == 1)
    errAbort("must specify one and only one of -db, -cdsDb, -cdsFile, or -allCds");

countFrameShift(db, cdsDb, cdsFile, pslSpec, geneFile, outputFile);
return 0;
}

