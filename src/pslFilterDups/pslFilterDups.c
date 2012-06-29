/* pslFilterDups - used for rna alignments , keep only best hits to a genomic region */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "memalloc.h"
#include "psl.h"
#include "cheapcgi.h"
#include "portable.h"

static char const rcsid[] = "$Id: pslFilterDups.c,v 1.1 2011/11/04 23:01:44 baertsch Exp $";

int reward = 1;
int cost = 1;
int gapOpenCost = 4;
double gapSizeLogMod = 1.0;
int minScore = 15;
int minMatch = 30;
int minUniqueMatch = 20;
int maxBadPpt = 700;
int minAli = 600;
int minAliT = 0;

void usage()
/* Print usage instructions and exit. */
{
errAbort(
    "pslFilterDups - used for rna alignments , keep only best hits to a genomic region \n"
    "    pslFilterDups in.psl out.psl \n"
    "options\n"
    "     -dir  Input files are directories rather than single files\n"
    "    -reward=N (default %d) Bonus to score for match\n"
    "    -cost=N (default %d) Penalty to score for mismatch\n"
    "    -gapOpenCost=N (default %d) Penalty for gap opening\n"
    "    -gapSizeLogMod=N (default %4.2f) Penalty for gap sizes\n"
    "    -minScore=N (default %d) Minimum score to pass filter\n"
    "    -minMatch=N (default %d) Min match (including repeats to pass)\n"
    "    -minUniqueMatch (default %d) Min non-repeats to pass)\n"
    "    -maxBadPpt (default %d) Maximum divergence in parts per thousand\n"
    "    -minAli (default %d) Minimum ratio query in alignment in ppt\n"
    "    -minAliT (default %d) Like minAli for target\n",
    reward, cost, gapOpenCost, gapSizeLogMod, minScore, minMatch,
    minUniqueMatch, maxBadPpt, minAli, minAliT);
}

boolean filterOk(struct psl *psl)
/* Return TRUE if psl passes filter. */
{
int score;
int totalAli = psl->match + psl->repMatch + psl->misMatch;
double totAli1000 = totalAli * 1000;

if (psl->match + psl->repMatch < minMatch)
    return FALSE;
if (psl->match < minUniqueMatch)
    return FALSE;
if (totAli1000 / psl->qSize < minAli)
    return FALSE;
if (totAli1000 / psl->tSize < minAliT)
    return FALSE;
if (pslCalcMilliBad(psl, FALSE) > maxBadPpt)
    return FALSE;
score = (psl->match + psl->repMatch)*reward - psl->misMatch*cost 
        - (psl->qNumInsert + psl->tNumInsert + 1) * gapOpenCost
	- log(psl->qBaseInsert + psl->tBaseInsert + 1) * gapSizeLogMod;
if (score < minScore)
    return FALSE;
return TRUE;
}

void removeOverlap(struct psl **pslList)
/* remove alignments that overlap the highest score alignment */
{
struct psl *top = NULL, *psl, *tmp;

if (pslList == NULL) return;
if (*pslList == NULL) return;
top = *pslList;
if (top->next == NULL) return;
for (psl = top->next; psl != NULL ; psl = psl->next)
    {
    if (rangeIntersection(top->tStart, top->tEnd, psl->tStart, psl->tEnd) > 0)
        {
        tmp = psl->next; /* save next element in list before delete*/
        slRemoveEl(pslList, psl);
        pslFree(&psl);
        psl = tmp; /* restore pointer to next element */
        if (psl == NULL) break;
        }
    }
    tmp = top->next;
    removeOverlap(&(tmp->next));
}
void pslFilterDups(char *inName, char *outName)
/* Filter inName into outName. */
{
struct lineFile *in = pslFileOpen(inName);
FILE *out = mustOpen(outName, "w");
struct psl *psl, *p, *pslList = NULL;
int passCount = 0;
int totalCount = 0;
char prev[128] = "first";

printf("Filtering %s to %s\n", inName, outName);
//pslWriteHead(out);
while ((psl = pslNext(in)) != NULL)
    {
    ++totalCount;
    if (sameString(prev,"first"))
        {
        pslList = psl;
        }
    else if (sameString(prev,psl->qName))
        {
        slAddHead(&pslList, psl);
        continue;
        }
    else
        {
        assert(pslList != NULL);
        slReverse(&pslList);
        removeOverlap(&pslList);
        for (p = pslList; p != NULL ; p = p->next)
            {
            ++passCount;
            pslTabOut(p,out);
            }
        pslFree(&pslList);
        pslList = psl;
        }
     
    /*
    if (filterOk(psl))
	{
	++passCount;
	pslTabOut(psl, out);
	}
    pslFree(&psl);
    */
    safef(prev, sizeof(prev), psl->qName);
    }
slReverse(&pslList);
removeOverlap(&pslList);
for (p = pslList; p != NULL ; p = p->next)
    pslTabOut(p,out);
pslFree(&pslList);
pslList = psl;
lineFileClose(&in);
fclose(out);
printf("%d of %d passed filter\n", passCount, totalCount);
}

void pslFilterDir(char *inDir, char *outDir)
/* Filter all .psl and .pslx files in directory. */
{
struct slName *inList = listDir(inDir, "*.psl*");
struct slName *inEl;
int inCount = slCount(inList), inIx = 0;
char inName[512], outName[512];
for (inEl = inList; inEl != NULL; inEl = inEl->next)
    {
    printf("%d of %d ", ++inIx, inCount);
    fflush(stdout);
    sprintf(inName, "%s/%s", inDir, inEl->name);
    sprintf(outName, "%s/%s", outDir, inEl->name);
    pslFilterDups(inName, outName);
    }
}

int main(int argc, char *argv[])
{
cgiSpoof(&argc, argv);
if (argc != 3)
    usage();
reward = cgiUsualInt("reward", reward);
cost = cgiUsualInt("cost", cost);
gapOpenCost = cgiUsualInt("gapOpenCost", gapOpenCost);
gapSizeLogMod = cgiUsualDouble("gapSizeLogMod", gapSizeLogMod);
minScore = cgiUsualInt("minScore", minScore);
minMatch = cgiUsualInt("minMatch", minMatch);
minUniqueMatch = cgiUsualInt("minUniqueMatch", minUniqueMatch);
maxBadPpt = cgiUsualInt("maxBadPpt", maxBadPpt);
minAli = cgiUsualInt("minAli", minAli);
if (cgiBoolean("dir") || cgiBoolean("dirs"))
    {
    pslFilterDir(argv[1], argv[2]);
    }
else
    {
    pslFilterDups(argv[1], argv[2]);
    }
return 0;
}
