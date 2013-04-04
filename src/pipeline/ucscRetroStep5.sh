#!/bin/bash 
# post process cluster jobs that runs retroFinder
# and then filters out zinc fingers and immunoglobin retros if PFAM is defined in the DEF file.
# the retro results, alignments and cds coding of the parent mapped to the retro is loaded into the database.
#
# script analyzeExpress.sh is called at the end to check for expression level of retros.
#
DEF=$1
source $DEF
#set -beEu -o pipefail
CWD=`pwd`
## called from ~baertsch/baertsch/scripts/ucscRetroStep4.sh
echo "-----------------------------------"
echo "Starting ucscRetroStep5.sh $DEF on $HOST"
echo "-----------------------------------"
date
cd $OUTDIR
pwd
echo cat chr*_NoOverlap.bed to pseudoGeneLinkNoOverlap.bed
cat chr*_NoOverlap.bed | $SCRIPT/removeTandemDups > pseudoGeneLinkNoOverlap.bed
wc -l chr*_NoOverlap.bed
tawk '$5>=350{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkNoOverlapFilter.bed                       
tawk '$5>=425{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLink425.bed
tawk '$5>=510{print $0}' pseudoGeneLinkNoOverlap.bed > retroMrnaInfo.raw.bed
tawk '$5>=650{print $0}' pseudoGeneLinkNoOverlap.bed > retroMrnaInfo650.bed
#cut -f 1-12 retroMrnaInfo650.bed > retroMrnaInfo.12.bed
wc -l pseudoGeneLinkNoOverlap.bed pseudoGeneLinkNoOverlapFilter.bed pseudoGeneLink425.bed retroMrnaInfo.raw.bed retroMrnaInfo650.bed 

if [[ -n $PFAM ]] 
then
   #remove immunoglobin
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where (a1.alias like 'IGH%' or a1.alias like 'IGK%' or a1.alias like 'IGL%') and a1.kgID = a2.kgID"  > kgImmuno.lst
    hgsql $DB -N -B -e "     select gbCdnaInfo.acc
          from gbCdnaInfo, geneName, description 
             where gbCdnaInfo.geneName=geneName.id and gbCdnaInfo.description = description.id and (geneName.name like 'IGH%' or geneName.name like 'IKL%' or geneName.name like 'IGL%')" >> kgImmuno.lst
    hgsql $DB -N -B -e "select name from knownGene where proteinID like 'IGH%' or proteinID like 'IKL%' or proteinID like 'IGL%'" >> kgImmuno.lst
#remove znf and NBPF
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'ZNF%' and a1.kgID = a2.kgID" |sort -u > kgZnf.lst
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'NBPF%' and a1.kgID = a2.kgID" |sort -u >> kgZnf.lst

#    cat kgZnf.lst refZnf.lst kgImmuno.lst > bothZnf.lst
cat kgZnf.lst kgImmuno.lst > bothZnf.lst

# grap genes with pfam domains (zinc finger, immunoglobin, NBPF, and olfactory receptor
    hgsql $DB -N -B -e "select  k.name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, \
                exonStarts, exonEnds from $GENEPFAM k, $PFAM p \
                where k.name = p.$PFAMIDFIELD and p.$PFAMDOMAIN in (\
            'PF00096', 'PF01352', 'PF06758', 'PF00047', 'PF07654', 'PF00642'  );" > zincKg.gp
    echo zcat $GENEPFAM.tab.gz to  $GENEPFAM.tab
    zcat $GENEPFAM.tab.gz > $GENEPFAM.tab
    overlapSelect zincKg.gp $GENEPFAM.tab zincKg2.gp -inFmt=genePred
    cut -f 1 zincKg2.gp |sort | uniq > zincKg.lst
    echo "Zinc fingers excluded"
    wc -l zincKg2.gp zincKg.lst

    echo "$SCRIPT/selectById 1 zincKg.lst 1 $GENEPFAM.tab to kgZnf.gp"
    $SCRIPT/selectById 1 zincKg.lst 1 $GENEPFAM.tab > kgZnf.gp
else
echo "xxx" > zincKg.lst
echo "skipping zinc finger and immunoglobin filtering"
fi


$SCRIPT/selectById 1 zincKg.lst 44 retroMrnaInfo.raw.bed > retroMrnaInfoZnf.bed
$SCRIPT/selectById -not 1 zincKg.lst 44 retroMrnaInfo.raw.bed > retroMrnaInfoLessZnf.bed
echo "before and after zinc finger filtering"
rm -f $TABLE.bed
cp -pf retroMrnaInfoLessZnf.bed $TABLE.bed
sort -k4,4 $TABLE.bed > $TABLE.sort.bed
sort -k4,4 ortho.txt | join -1 4 -2 1 -o 2.1 2.2 2.3 $TABLE.sort.bed -|awk '{$1=$1;OFS="\t";print $0}' > ortho.filter.txt
wc -l retroMrnaInfo.raw.bed retroMrnaInfoLessZnf.bed retroMrnaInfoZnf.bed $TABLE.bed
cut -f 1-12 $TABLE.bed > retroMrnaInfo.12.bed
textHistogram -col=5 $TABLE.bed -binSize=50 -maxBinCount=50
echo creating $ALIGN.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' $TABLE.bed > pseudoGeneLinkSelect.tab
pslSelect -qtStart=pseudoGeneLinkSelect.tab pseudo.psl $ALIGN.psl
wc -l $ALIGN.psl pseudoGeneLinkSelect.tab
hgLoadBed $DB -verbose=9 -renameSqlTable -allowNegativeScores -noBin ucscRetroInfoXX -sqlTable=/cluster/home/$USER/kent/src/hg/lib/ucscRetroInfo.sql $TABLE.bed
mkdir -p $RETRODIR
rm -f $RETRODIR/$TABLE.bed
cp -p $TABLE.bed $RETRODIR
hgsql $DB -e "drop table if exists $TABLE;"
hgsql $DB -e "alter table ucscRetroInfoXX rename $TABLE;"
hgLoadPsl $DB $ALIGN.psl
rm -f $RETRODIR/$ALIGN.psl
cp -p $ALIGN.psl $RETRODIR
hgLoadSqlTab $DB ${ORTHOTABLE} ~/kent/src/hg/lib/ucscRetroOrtho.sql ortho.filter.txt

zcat cds.tab.gz |tawk '{print $1"."$2,$3}' | sort | uniq > ucscRetroCds${VERSION}.tab 
hgLoadSqlTab $DB ucscRetroCds${VERSION} ~/kent/src/hg/lib/ucscRetroCds.sql ucscRetroCds${VERSION}.tab
rm -f $RETRODIR/ucscRetroCds${VERSION}.tab
cp -p ucscRetroCds${VERSION}.tab $RETRODIR
cp -p DEF $RETRODIR

#writing TrackDb.ra entry to temp file
$SCRIPT/makeTrackDb.sh $DEF > $OUTDIR/trackDb.retro
echo "writing template trackDb.ra entry to $OUTDIR/trackDb.retro"

echo "database loaded, update trackDb.ra entry"
echo "run $SCRIPT/analyseExpress.sh $DEF to update expression level of each retro"
