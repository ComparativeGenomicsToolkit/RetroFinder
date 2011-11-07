#!/bin/bash
set -bevEu -o pipefail
source $1
echo '-------- script analyseAge.sh -------------------'
#########################
# analyze age of retros 
#########################

echo "calculate age of retros"
cd $OUTDIR
cd $EXPDIR

#split retros by age

for bed in retroMrnaInfo650 retroMrnaInfoLessZnf ; do $SPLITBYAGE $OUTDIR/${bed}.bed $OUTDIR/ortho.filter.txt ${bed}.ancient.bed ${bed}.recent.bed $ANCIENT1 $ANCIENT2 $SCRIPT ; wc -l ${bed}.ancient.bed ${bed}.recent.bed ; done

for i in Coding Noncoding Nearcoding Antisense ; do $SCRIPT/selectById -tsv 1 ortho.ancient.id 4 retroKg${i}.bed |sort -k5,5nr> retroKg${i}.ancient.bed ; $SCRIPT/selectById -not -tsv 1 ortho.ancient.id 4 retroKg${i}.bed |sort -k5,5nr> retroKg${i}.recent.bed ; done

for i in ${GENE1}Cds ${GENE2}Cds ${GENE1}Cds50 ${GENE2}Cds50 ; do $SCRIPT/selectById -tsv 1 ortho.ancient.id 4 pseudo${i}.bed |sort -k5,5nr> pseudo${i}.ancient.bed ; $SCRIPT/selectById -not -tsv 1 ortho.ancient.id 4 pseudo${i}.bed |sort -k5,5nr> pseudo${i}.recent.bed ; done

if [[ $GENE1 == "knownGene" ]] 
then
for i in ${GENE1}ProtCds ${GENE2}ProtCds; do $SCRIPT/selectById -tsv 1 ortho.ancient.id 4 pseudo${i}.bed |sort -k5,5nr> pseudo${i}.ancient.bed ; $SCRIPT/selectById -not -tsv 1 ortho.ancient.id 4 pseudo${i}.bed |sort -k5,5nr> pseudo${i}.recent.bed ; done
fi
for i in ucscRetroProtein ucscRetroTranscript ; do $SCRIPT/selectById -tsv 1 ortho.ancient.id 4 ${i}.bed |sort -k5,5nr> ${i}.ancient.bed ; $SCRIPT/selectById -not -tsv 1 ortho.ancient.id 4 ${i}.bed |sort -k5,5nr> ${i}.recent.bed ; done

for bed in pseudo${GENE2} pseudoEstMrna.filter pseudoEstAll pseudoExpressed pseudoEst5Mrna pseudoEst5 pseudoEst5AndMrna pseudoEst100AA pseudo5Est100AA pseudoEst5AndMrna100AA pseudo${GENE1}Cds pseudo${GENE1}Cds50 pseudo${GENE2}Cds pseudo${GENE2}Cds50 pseudoEst5MrnaNotKg pseudoEst10Mrna pseudoEst10MrnaNotKg pseudoEst5AndMrnaNotKg; do $SCRIPT/selectById -tsv 1 ortho.ancient.id 4 ${bed}.bed |sort -k5,5nr> ${bed}.ancient.bed ; $SCRIPT/selectById -not -tsv 1 ortho.ancient.id 4 ${bed}.bed |sort -k5,5nr> ${bed}.recent.bed ; done

#for i in Coding Noncoding Nearcoding Antisense ; do $SPLITBYAGE retroKg${i}.bed retroKg${i}.ancient.bed retroKg${i}.recent.bed; done
#TODO: convert these to selectById
#for i in ${GENE1}Cds ${GENE2}Cds ${GENE1}Cds50 ${GENE2}Cds50 ; do $SPLITBYAGE pseudo${i}.bed pseudo${i}.ancient.bed pseudo${i}.recent.bed ; done


if [ -f altSplice.bed ] 
then
 $SPLITBYAGE retroSplice.bed retroSplice.ancient.bed retroSplice.recent.bed 
fi

$SCRIPT/doAge ../DEF pseudoEst5Mrna
$SCRIPT/doAge ../DEF retroMrnaInfoLessZnf 
$SCRIPT/doAge ../DEF retroKgCoding 


echo '-------- END script analyseAge.sh -------------------'
echo '------ run script ucscRetroStep6.sh to make Html pages -------------------'
