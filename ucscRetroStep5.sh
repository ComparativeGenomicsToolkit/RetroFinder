#!/bin/bash 
# post process cluster jobs that runs retroFinder
# and then filters out zinc fingers and immunoglobin retros if PFAM is defined in the DEF file.
# the retro results, alignments and cds coding of the parent mapped to the retro is loaded into the database.
#
# script analyzeExpress.sh is called at the end to check for expression level of retros.
#
source $1
#set -beEu -o pipefail
CWD=`pwd`
## called from ~baertsch/baertsch/scripts/ucscRetroStep4.sh
echo "-----------------------------------"
echo "Starting ucscRetroStep5.sh on $HOST"
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

    cat kgZnf.lst refZnf.lst kgImmuno.lst > bothZnf.lst

# grap genes with pfam domains (zinc finger, immunoglobin, NBPF, and olfactory receptor
    hgsql $DB -N -B -e "select  k.name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENEPFAM k, $PFAM p \
                where k.name = p.$PFAMIDFIELD and p.$PFAMDOMAIN in (\
            'PF00096', 'PF01352', 'PF06758', 'PF00047', 'PF07654'  );" > zincKg.gp
    cut -f 1 zincKg.gp |sort | uniq > zincKg.lst
    echo "Zinc fingers excluded"
    wc -l zincKg.gp zincKg.lst

    echo zcat $GENEPFAM.tab.gz to  $GENEPFAM.tab
    zcat $GENEPFAM.tab.gz > $GENEPFAM.tab
    echo "$SCRIPT/selectById 1 zincKg.lst 4 $GENEPFAM.tab to kgZnf.gp"
    $SCRIPT/selectById 1 zincKg.lst 1 $GENEPFAM.tab > kgZnf.gp
else
echo "xxx" > zincKg.lst
echo "skipping zinc finger and immunoglobin filtering"
fi


grep -v -F -f zincKg.lst retroMrnaInfo.raw.bed > retroMrnaInfoLessZnf.bed
echo "before and after zinc finger filtering"
cp retroMrnaInfoLessZnf.bed $TABLE.bed
wc -l retroMrnaInfo.raw.bed retroMrnaInfoLessZnf.bed $TABLE.bed
cut -f 1-12 $TABLE.bed > retroMrnaInfo.12.bed
textHistogram -col=5 $TABLE.bed -binSize=50 -maxBinCount=50
echo creating $ALIGN.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' $TABLE.bed > pseudoGeneLinkSelect.tab
pslSelect -qtStart=pseudoGeneLinkSelect.tab pseudo.psl $ALIGN.psl
wc -l $ALIGN.psl pseudoGeneLinkSelect.tab
hgLoadBed $DB -verbose=9 -allowNegativeScores -noBin retroMrnaInfoXX -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/retroMrnaInfo.sql $TABLE.bed
hgsql $DB -e "drop table if exists $TABLE;"
# kaku is no longer being used and some values are not loaded correctly
# as represented in exponential notation so they get replaced by inf.
hgsql $DB -e "alter table retroMrnaInfoXX set kaku = 0;"
hgsql $DB -e "alter table retroMrnaInfoXX rename $TABLE;"
hgLoadPsl $DB $ALIGN.psl

hgLoadSqlTab $DB ucscRetroCds ~/kent/src/hg/lib/ucscRetroCds.sql ucscRetroCds.tab

####################
# analyze expression 
####################
#make exon list
cd $OUTDIR
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 where exonCount > 1" > $GENE1.multiExon.genePred
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2 where exonCount > 1" > $GENE2.multiExon.genePred
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3 where exonCount > 1" > $GENE3.multiExon.genePred
echo "genePredFilter $GENE1.multiExon.genePred $GENE1.multiCDSExon.genePred -cdsExons=2"
/cluster/home/baertsch/bin/x86_64/genePredFilter $GENE1.multiExon.genePred $GENE1.multiCDSExon.genePred -cdsExons=2
echo "genePredFilter $GENE2.multiExon.genePred $GENE2.multiCDSExon.genePred -cdsExons=2"
/cluster/home/baertsch/bin/x86_64/genePredFilter $GENE2.multiExon.genePred $GENE2.multiCDSExon.genePred -cdsExons=2
echo "genePredFilter $GENE3.multiExon.genePred $GENE3.multiCDSExon.genePred -cdsExons=2"
/cluster/home/baertsch/bin/x86_64/genePredFilter $GENE3.multiExon.genePred $GENE3.multiCDSExon.genePred -cdsExons=2
wc -l $GENE1.multiCDSExon.genePred $GENE2.multiCDSExon.genePred $GENE1.multiCDSExon.genePred
echo "cat $GENE1.multiCDSExon.genePred $GENE2.multiCDSExon.genePred $GENE3.multiCDSExon.genePred to  all.multiCds.gp"
cat $GENE1.multiCDSExon.genePred $GENE2.multiCDSExon.genePred $GENE3.multiCDSExon.genePred > all.multiCds.gp

echo "ldHgGene $DB rbRefGeneMulti $GENE2.multiCDSExon.genePred -predTab"
ldHgGene $DB rbRefGeneMulti $GENE2.multiCDSExon.genePred -predTab
featureBits $DB $GENE2:cds -bed=$GENE2.cds.bed
featureBits $DB rbRefGeneMulti;
#52619505 bases of 2881515245 (1.826%) in intersection
featureBits $DB $GENE2;     
#54738314 bases of 2881515245 (1.900%) in intersection
featureBits $DB rbRefGeneMulti $GENE2.cds.bed -bed=$GENE2.multiCds.bed 
#29845142 bases of 2881515245 (1.036%) in intersection

pwd
ls -l retroMrnaInfoLessZnf.bed 
echo "calculate age of retros"

for bed in retroMrnaInfoLessZnf retroMrnaInfo650 ; do $SPLITBYAGE ${bed}.bed ${bed}.ancient.bed ${bed}.recent.bed; done
tawk '{print NF}' retroMrnaInfoLessZnf.bed |uniq

echo "overlap with ests - use cluster jobs to speed this step"
mkdir -p estSplit
pslSplitOnTarget estFiltered.psl.gz estSplit
cd estSplit
rm -f jobList
for i in `ls *.psl` ; do echo overlapSelect $i ../retroMrnaInfoLessZnf.bed pseudoEst.${i%%.psl}.bed ; done >> jobList
cd ..
ssh -T $CLUSTER "cd $OUTDIR/estSplit ; /parasol/bin/para make jobList"

cat estSplit/pseudoEst*.bed > pseudoEstAll.bed
wc -l retroMrnaInfoLessZnf.bed pseudoEstAll.bed

echo "mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid "
mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid > all_mrna.log 2> all_mrna.err
mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 -utrMergeSize=10 all_mrna_utr.gp -cdsDb=$DB -keepInvalid > all_mrna_utr.log 2> all_mrna_utr.err

awk '$8>1{print }' all_mrna.gp > all_mrna_multiExon.gp
awk '$8>1{print }' all_mrna_utr.gp > all_mrna_multiExonUTR.gp

mkdir -p $EXPDIR
cd $EXPDIR
echo pwd
pwd
$SCRIPT/analyseExpress.sh ../$1
cd ..
