/* hgCountTissueMicroarray - Create a copy of microarray database that contains the median value of replicas. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "obscure.h"
#include "portable.h"
#include "jksql.h"
#include "hdb.h"
#include "hgRelate.h"
#include "expRecord.h"
#include "expData.h"
#include "microarray.h"

static char const rcsid[] = "$Id: hgCountTissueMicroarray.c,v 1.1 2009/11/27 00:00:53 baertsch Exp $";

char *tabDir = ".";
boolean doLoad;
int limit = 0;
float cutoff = 2.0;
int minExps = 2;

void usage()
/* Explain usage and exit. */
{
errAbort(
  "hgCountTissueMicroarray - count number of tissues expressed per probe higher than some threshold \n"
  "usage:\n"
  "   hgCountTissueMicroarray db table expTable specFile.ra output\n"
  "The specFile.ra contains a line for each expression column you want in\n"
  "the new table.  The format of the line is\n"
  " 'name' 'group' ix1 ix2 ... ixN\n"
  "where 'name' is what will appear in the column label in the Family browser\n"
  "and 'group' is what whill appear in the row label in dense mode in the\n"
  "genome browser, and ix1 ... ixN are the (zero based) indices of columns you\n"
  "are taking the median of in the original data set\n"
  "options:\n"
  "   -cutoff=N.N - Only count tissues that have expression greater than N.N (default %f4.2)\n"
  "   -tab=dir - Output tab-separated files to directory.\n"
  "   -noLoad  - If true don't load database and don't clean up tab files\n"
  "   -limit=N - Only do limit rows of table, for testing\n"
  "   -minExps=N - Minimum number of experiments to use, default %d\n"
  , cutoff, minExps
  );
}

static struct optionSpec options[] = {
   {"tab", OPTION_STRING},
   {"cutoff", OPTION_FLOAT},
   {"noLoad", OPTION_BOOLEAN},
   {"limit", OPTION_INT},
   {"minExps", OPTION_INT},
   {NULL, 0},
};

void makeNewExpTable(char *oldTable, struct maMedSpec *medList, char *newTable)
/* Create new expTable in hgFixed that is very similar
 * to oldExpTable, but with rows defined by medList. */
{
struct maMedSpec *med;
struct expRecord *oldExp, newExp;
struct sqlConnection *conn = sqlConnect("hgFixed");
FILE *f = hgCreateTabFile(tabDir, newTable);
char query[256], **row;
struct sqlResult *sr;
int curId = 0;

for (med = medList; med != NULL; med = med->next)
    {
    /* Load expression record from old table of first
     * thing in median. */
    safef(query, sizeof(query),
    	"select * from %s where id = %d", oldTable, med->ids[0]);
    sr = sqlGetResult(conn, query);
    if ((row = sqlNextRow(sr)) == NULL)
        errAbort("Can't find id %d in %s\n", med->ids[0], oldTable);
    oldExp = expRecordLoad(row);
    sqlFreeResult(&sr);
    if (oldExp->numExtras < 3)
        errAbort("Can only deal with old tables with 3 extras or more");


    /* Create new expression record, mostly just a shallow copy of old. */
    newExp = *oldExp;
    newExp.id = curId;
    ++curId;
    newExp.name = newExp.description = med->name;
    newExp.extras[2] = med->group;

    /* Save new one, free old one. */
    expRecordTabOut(&newExp, f);
    expRecordFree(&oldExp);
    }

if (doLoad)
    {
    expRecordCreateTable(conn, newTable);
    hgLoadTabFile(conn, tabDir, newTable, &f);
    hgRemoveTabFile(tabDir, newTable);
    }
sqlDisconnect(&conn);
}

struct expData *countTissuesFromSpec(struct expData *expList, struct maMedSpec *specList, 
        float cutoff, int minExps)
/* count number of tissues that have log ratio >= cutoff 
   return list of probes and counts */
{
struct expData *exp , *retList = NULL;

for (exp = expList ; exp != NULL ; exp = exp->next)
    {
    struct expData *countExp;
    int sizeOne;
    char *row[1] ;
    struct maMedSpec *spec;
    int j;
    row[0] = "0.0";
    AllocVar(countExp);
    countExp->expCount = 1;
    countExp->name = cloneString(exp->name);
    sqlFloatDynamicArray(row[0], &countExp->expScores, &sizeOne);
    assert(sizeOne == countExp->expCount);
    for (spec = specList ; spec != NULL ; spec = spec->next)
        for (j = 0; j < spec->count; j++)
            {
            int id = spec->ids[j];
            assert (id <= exp->expCount);
            /* look for positive values */
            if (exp->expScores[id] >= cutoff)
                countExp->expScores[0]++;
            /* look for negative values */
            if (exp->expScores[id] <= -cutoff)
                countExp->expScores[0]++;
            }
    slAddHead(&retList, countExp);
    }
return retList;
}

void makeNewDataTable(char *database, char *oldTable, struct maMedSpec *medList, char *newTable)
/* Create new table in database based on medians of data
 * in old table as defined by medList. */
{
struct sqlConnection *conn = sqlConnect(database);
FILE *f = hgCreateTabFile(tabDir, newTable);
struct expData *expList, *tissueCountList, *exp;

expList = expDataLoadTableLimit(conn, oldTable, limit);
tissueCountList = countTissuesFromSpec(expList, medList, cutoff, minExps);
for (exp = tissueCountList; exp != NULL; exp = exp->next)
    expDataTabOut(exp, f);
/*
if (doLoad)
    {
    expDataCreateTable(conn, newTable);
    hgLoadTabFile(conn, tabDir, newTable, &f);
    hgRemoveTabFile(tabDir, newTable);
    }
    */
expDataFreeList(&expList);
//expDataFreeList(&medianExpList);
sqlDisconnect(&conn);
}


void hgCountTissueMicroarray(char *database, char *expTable, 
	char *specFile, char *newTable, char *newExpTable)
/* hgCountTissueMicroarray - count number of tissues expressed per probe higher than some threshold */
{
struct maMedSpec *medList = maMedSpecReadAll(specFile);
//makeNewExpTable(expTable, medList, newExpTable);
makeNewDataTable(database, expTable, medList, newTable);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 5)
    usage();
doLoad = !optionExists("noLoad");
minExps = optionInt("minExps", minExps);
cutoff = optionFloat("cutoff", cutoff);
if (optionExists("tab"))
    {
    tabDir = optionVal("tab", tabDir);
    makeDir(tabDir);
    }
limit = optionInt("limit", limit);

hgCountTissueMicroarray(argv[1], argv[2], argv[3], argv[4], argv[5]);
return 0;
}
