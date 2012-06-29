/* genePredFilter - filter genePred files */
#include "common.h"
#include "options.h"
#include "chromBins.h"
#include "binRange.h"
#include "genePred.h"
#include "genePredReader.h"
#include "hash.h"
#include "localmem.h"
#include "linefile.h"
#include "verbose.h"

static char const rcsid[] = "";

int cdsExons = 1000;

/* Command line option specifications */
static struct optionSpec optionSpecs[] = {
    {"cdsExons", OPTION_INT},
    {NULL, 0}
};

static boolean cdsOverlaps(struct genePred *gp1, struct genePred *gp2)
/* determine if there is any overlap in the CDS of two genes */
{
int iExon1, start1, end1, iExon2, start2, end2;

for (iExon1 = 0; iExon1 < gp1->exonCount; iExon1++)
    {
    if (genePredCdsExon(gp1, iExon1, &start1, &end1))
        {
        for (iExon2 = 0; iExon2 < gp2->exonCount; iExon2++)
            {
            if (genePredCdsExon(gp2, iExon2, &start2, &end2))
                if (positiveRangeIntersection(start1, end1, start2, end2))
                    return TRUE;
            }
        }
    }
return FALSE;
}

struct genePred *findOverGenes(struct binKeeper *chrGenes,
                               struct genePred *seedGp)
/* find all genes who's CDS overlaps the specified gene and remove them
 * from chrGenes.  seedGp must already have been removed from chrGenes. */
{
struct genePred *overCds = NULL;
struct binElement *overs = binKeeperFind(chrGenes, seedGp->cdsStart, seedGp->cdsEnd);
struct binElement *gpEl;

while ((gpEl = slPopHead(&overs)) != NULL)
    {
    struct genePred *gp = gpEl->val; 
    if (cdsOverlaps(seedGp, gp))
        {
        binKeeperRemove(chrGenes, gp->cdsStart, gp->cdsEnd, gp);
        slAddHead(&overCds, gp);
        }
    freez(&gpEl);
    }
return overCds;
}

struct genePred *getOverGenes(struct binKeeper *chrGenes, struct genePred *nextGp)
/* get list of genes who's CDS overlaps the specified gene.  The returned
 * list will included the supplied gene and all genes will be removed from
 * chrGenes */
{
struct genePred *overGenes = nextGp, *gp;
boolean anyFound = TRUE;
binKeeperRemove(chrGenes, nextGp->cdsStart, nextGp->cdsEnd, nextGp);

/* continue to grow the list until no more overlapping are found.
 * this requires rechecking the expanded list each pass */
do
    {
    anyFound = FALSE;
    for (gp = overGenes; gp != NULL; gp = gp->next)
        {
        struct genePred *moreOver = findOverGenes(chrGenes, gp);
        if (moreOver != NULL)
            {
            anyFound = TRUE;
            overGenes = slCat(overGenes, moreOver);
            }
        }
    }
while (anyFound);

return overGenes;
}


/* return count of coding exons */
int genePredCountCodingExons(struct genePred *gp)
{
int i;
int count = 0;
for (i=0; i<(gp->exonCount); i++)
    {
    if ( (gp->cdsStart <= gp->exonEnds[i]) &&  
         (gp->cdsEnd >= gp->exonStarts[i]) )
         count++;
    }
return count;
}

static void genePredFilter(char *inGpFile, char *outGpFile)
/* create single-coverage genePred files */
{
struct genePredReader *gpr = genePredReaderFile(inGpFile, NULL);
struct genePred *gp;
FILE *outFh = mustOpen(outGpFile, "w");

while ((gp = genePredReaderNext(gpr)) != NULL)
    if (genePredCountCodingExons(gp) >= cdsExons)
        genePredTabOut(gp, outFh);

carefulClose(&outFh);
}

static void usage(char *msg)
/* Explain usage and exit. */
{
errAbort("%s\n"
    "\n"
    "genePredFilter - filter genePred files\n"
    "\n"
    "genePredFilter [options] inGenePred outGenePred\n"
    "\n"
    "Options:\n"
    "  -cdsExons=N - require at least N coding exons.\n"
    "\n", msg);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, optionSpecs);
cdsExons = optionInt("cdsExons" , cdsExons) ;
if (argc != 3)
    usage("wrong # args");
genePredFilter(argv[1], argv[2]);

return 0;
}
