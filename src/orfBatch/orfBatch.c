/* orfBatch - Predict best orf on a list of regions using borf. */
#include "common.h"
#include "hdb.h"
#include "linefile.h"
#include "bed.h"
#include "borf.h"
#include "ucscRetroInfo.h"
#include "genePred.h"


static char const rcsid[] = "$Id: orfBatch.c,v 1.2 2011/11/05 21:56:02 baertsch Exp $";
float minScore = 50.0;

void usage()
/* Explain usage and exit. */
{
errAbort(
  "orfBatch - Predict best orf on a list of regions using borf\n"
  "usage:\n"
  "   orfBatch database predict.bed output.bed output.gp\n"
  );
}


void bedMergeBlocks(struct bed *bed, struct bed *outBed,
                       int insertMergeSize)
/* merge together blocks separated by small inserts. */
{
int iBlk, iExon = -1;
int startIdx, stopIdx, idxIncr;

assert(outBed!=NULL);
outBed->chromStarts = needMem(bed->blockCount*sizeof(unsigned));
outBed->blockSizes = needMem(bed->blockCount*sizeof(unsigned));

//if (psl->strand[1] == '-')
//    {
//    startIdx = psl->blockCount-1;
//    stopIdx = -1;
//    idxIncr = -1;
//    }
//else
//    {
    startIdx = 0;
    stopIdx = bed->blockCount;
    idxIncr = 1;
//    }

for (iBlk = startIdx; iBlk != stopIdx; iBlk += idxIncr)
    {
    unsigned tStart = bed->chromStarts[iBlk];
    unsigned size = bed->blockSizes[iBlk];
//    if (psl->strand[1] == '-')
//        reverseIntRange(&tStart, &tEnd, psl->tSize);
    if ((iExon < 0) || ((tStart - (outBed->chromStarts[iExon]+outBed->blockSizes[iExon])) > insertMergeSize))
        {
        iExon++;
        outBed->chromStarts[iExon] = tStart;
        outBed->blockSizes[iExon] = size;
	}
    else
        outBed->blockSizes[iExon] += size;
    }
outBed->blockCount = iExon+1;
strcpy(outBed->strand, bed->strand);
outBed->name = cloneString(bed->name);
outBed->thickStart = bed->thickStart;
outBed->thickEnd = bed->thickEnd;
outBed->chrom = cloneString(bed->chrom);
outBed->chromStart = bed->chromStart;
outBed->chromEnd = bed->chromEnd;
}
struct genePred *genePredFromBorfBed(struct bed *bed, struct borf *borf)
/* fill in genePred from bed and borf annotation */
{
char score[128];
int i;
struct genePred *gp = bedToGenePred(bed);

gp->optFields |= genePredName2Fld;
sprintf(score,"%4.1f",borf->score);
gp->name2 = cloneString(score);
gp->optFields |= genePredCdsStatFld;
if (borf != NULL )
    {
    gp->cdsStartStat = cdsComplete;
    gp->cdsEndStat = cdsComplete;
    }
else
    {
    gp->cdsStartStat = cdsUnknown;
    gp->cdsEndStat = cdsUnknown;
    }
if (gp->cdsStart < gp->txStart)
    {
    gp->cdsStart = gp->txStart;
    gp->cdsStartStat = cdsIncomplete;
    }
if (gp->cdsEnd > gp->txEnd)
    {
    gp->cdsEnd = gp->txEnd;
    gp->cdsEndStat = cdsIncomplete;
    }
if (borf != NULL)
    {
    AllocArray(gp->exonFrames, gp->exonCount);
    gp->optFields |= genePredExonFramesFld ;
    for (i=0; i<gp->exonCount; ++i)
        {
        gp->exonFrames[i]=0; //abs(atoi(borf->frame))-1;
        }
    }
return gp;
}

void doIt(char *database, char *fileName, char *outBed, char *geneFile)
/* read all bed records and generate genePred and update thickStart/End in bed */
{
struct bed *bed ;
struct bed *mBed ;
struct borf *borf = NULL;
FILE *outF = mustOpen(outBed,"w");
FILE *predF = mustOpen(geneFile,"w");
struct ucscRetroInfo *pg = NULL, *pgList = NULL;
char name[256];

//hSetDb(database);
pgList = ucscRetroInfoLoadAll(fileName);
    //bed = bedLoadN(row, wordCount);
for (pg = pgList ; pg != NULL ; pg = pg->next)
    {
    struct genePred *gp = NULL;
    //int lastBlock = pg->blockCount-1;
    AllocVar(bed);
    bed->chrom = pg->chrom;
    bed->chromStart = pg->chromStart;
    bed->chromEnd = pg->chromEnd;
    bed->score = pg->score;
    bed->blockCount = 0;
    if (pg->strand[0] == '+')
        bed->chromEnd = pg->chromEnd;
    else
        bed->chromEnd = pg->chromEnd-1;
    /*
    bed->chromStarts = needMem(bed->blockCount*sizeof(unsigned));
    bed->blockSizes = needMem(bed->blockCount*sizeof(unsigned));
    bed->chromStarts[0] = pg->chromStart;
    bed->blockSizes[0] = pg->chromEnd-pg->chromStart;
    bed->blockCount = pg->blockCount;
    bed->blockSizes[i] = pg->blockSizes[i];
    for (i = 0 ; i<(bed->blockCount) ; i++)
        {
        bed->chromStarts[i] = pg->chromStarts[i];
        bed->blockSizes[i] = pg->blockSizes[i];
        }
    if (pg->strand[0] == '+')
        {
        if (bed->chromEnd != (bed->chromStart)+(bed->chromStarts[lastBlock])+(bed->blockSizes[lastBlock]))
            bed->chromEnd = (bed->chromStart)+(bed->chromStarts[lastBlock])+(bed->blockSizes[lastBlock]);
        }
    else
        {
        if (bed->chromEnd != (bed->chromStart)+(bed->chromStarts[lastBlock])+(bed->blockSizes[lastBlock])-1)
            bed->chromEnd = (bed->chromStart)+(bed->chromStarts[lastBlock])+(bed->blockSizes[lastBlock]-1);
        }
        */
    safef(bed->strand ,sizeof(bed->strand),"%s",pg->strand);
    sprintf(name,"%s",pg->name);
    //sprintf(name,"%s.%s.%d",pg->name, pg->chrom, pg->chromStart);
    bed->name = name;
    if (bed->chromStart > bed->chromEnd)
        errAbort("start after end %d %d > %d %d of %s %s", pg->chromStart, bed->chromStart, pg->chromEnd, bed->chromEnd, pg->chrom, pg->name);
    if (bed->blockCount > 1)
        {
        AllocVar(mBed);
        bedMergeBlocks(bed, mBed, 10);
        }
    else
        mBed = bed;
    borf = borfFromGenomeBed(database, mBed);
    borfOutput(borf , stdout, '\t', '\n');
    /* set thickStart/thickEnd based on cds */
    if (pg->strand[0] == '+')
        {
        mBed->thickStart = pg->thickStart= borf->cdsStart+pg->chromStart;
        mBed->thickEnd = pg->thickEnd= borf->cdsEnd+pg->chromStart+3;
        }
    else
        {
        mBed->thickStart = pg->thickStart= borf->size-borf->cdsEnd+pg->chromStart-3;
        mBed->thickEnd = pg->thickEnd= (borf->size-borf->cdsStart+pg->chromStart);
        }
    printf("blockcount before %d after %d\n",bed->blockCount, mBed->blockCount);
    borf->feature = cloneString(pg->chrom);
    borf->cdsStart = pg->chromStart;
    borf->cdsEnd = pg->chromEnd;
    gp = genePredFromBorfBed(mBed, borf);
    genePredOutput(gp, predF, '\t', '\n');
    pg->posConf = round(borf->score);
    ucscRetroInfoOutput(pg, outF, '\t', '\n');
    genePredFree(&gp);
    }
carefulClose(&outF);
}

int main(int argc, char *argv[])
/* Process command line. */
{
if (argc < 5)
    usage();
doIt(argv[1], argv[2], argv[3], argv[4]);
return 0;
}
