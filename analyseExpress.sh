#!/bin/bash
set -beEu -o pipefail
source $1
echo '-------- script analyseExpress.sh -------------------'
#expression according to knownGene
hgsql $PDB -N -B -e "select * from hgncXref " > hgncXref.tab
if [[ $GENE1 == "knownGene" ]] 
then
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'coding'" > kgCoding.gp
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'noncoding'" > kgNoncoding.gp
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'nearcoding'" > kgNearcoding.gp
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'antisense'" > kgAntisense.gp
for i in Coding Noncoding Nearcoding Antisense ; do overlapSelect kg${i}.gp ../$TABLE.bed retroKg${i}.bed; $SPLITBYAGE retroKg${i}.bed retroKg${i}.ancient.bed retroKg${i}.recent.bed; done
zcat ../knownGene.tab.gz |sort > ../knownGene.sort.gp
hgsql $DB -N -B -e "select distinct kgID from kgSpAlias, uniProt.proteinEvidence where spId = acc and proteinEvidenceType <=1" |sort> kgProtEvidence.id
hgsql $DB -N -B -e "select distinct kgID from kgSpAlias, uniProt.proteinEvidence where spId = acc and proteinEvidenceType <=2" |sort> kgTransEvidence.id
$SCRIPT/selectById -tsv 1 kgProtEvidence.id 1 ../knownGene.sort.gp > knownGeneProt.gp
$SCRIPT/selectById -tsv 1 kgTransEvidence.id 1 ../knownGene.sort.gp > knownGeneTrans.gp

overlapSelect knownGeneProt.gp ../${GENE1}.cds.bed knownGeneProtCds.bed -selectFmt=genePred
overlapSelect knownGeneTrans.gp ../${GENE1}.cds.bed knownGeneTransCds.bed -selectFmt=genePred
overlapSelect knownGeneProtCds.bed ../$TABLE.bed ucscRetroProtein.bed
overlapSelect knownGeneTransCds.bed ../$TABLE.bed ucscRetroTranscript.bed
fi

echo "Exon Shuffling"
# exon shuffling
cat ../$TABLE.bed ../retroMrnaInfoZnf.bed > ucscRetroInfoZnf.bed

overlapSelect ../${GENE1}.multiCds.bed ../$TABLE.bed pseudo${GENE1}Cds.bed
overlapSelect ../${GENE1}.multiCds.bed ../$TABLE.bed -statsOutput pseudo${GENE1}Cds.out
overlapSelect ../${GENE1}.multiCds.bed ../$TABLE.bed pseudo${GENE1}Cds50.bed -overlapThreshold=0.50
overlapSelect ../${GENE2}.multiCds.bed ../$TABLE.bed pseudo${GENE2}Cds.bed
overlapSelect ../${GENE2}.multiCds.bed ../$TABLE.bed -statsOutput pseudo${GENE2}Cds.out
overlapSelect ../${GENE2}.multiCds.bed ../$TABLE.bed pseudo${GENE2}Cds50.bed -overlapThreshold=0.50
overlapSelect -selectFmt=genePred ../${GENE2}.tab.gz pseudo${GENE2}Cds.bed shuffleEns.bed
overlapSelect -selectFmt=genePred ../${GENE2}.multiCDSExon.genePred pseudo${GENE2}Cds.bed shuffleEnsMulti.bed
overlapSelect pseudo${GENE1}Cds.bed ../pseudoEstAll.bed pseudoEstAllNotShuffle.bed -nonOverlapping
for i in ${GENE1}Cds ${GENE2}Cds ${GENE1}Cds50 ${GENE2}Cds50 ; do $SPLITBYAGE pseudo${i}.bed pseudo${i}.ancient.bed pseudo${i}.recent.bed ; done

if [[ $GENE1 == "knownGene" ]] 
then
overlapSelect knownGeneProt.gp pseudo${GENE1}Cds.bed pseudo${GENE1}ProtCds.bed
overlapSelect knownGeneProt.gp pseudo${GENE2}Cds.bed pseudo${GENE2}ProtCds.bed
for i in ${GENE1}ProtCds ${GENE2}ProtCds; do $SPLITBYAGE pseudo${i}.bed pseudo${i}.ancient.bed pseudo${i}.recent.bed ; done
fi

ls -l pseudo${GENE1}*.bed
wc -l pseudo${GENE1}*.bed
cut -f 4 pseudo${GENE1}Cds.bed > shuffle.list
awk '{print "update '$TABLE' set type = \"expressed shuffle\" where name = \""$1"\" ;"}' shuffle.list> updateShuffle.sql
hgsql $DB < updateShuffle.sql

# expression analysis based on mrna and est overlap
# grab results of cluster job
cat ../estSplit/pseudoEst.*.bed > ../pseudoEstAll.bed
cat ../estSplit/est.*.id > est.id
cat ../estSplit/stat.*.out > stat.out
cat ../estSplit/statagg.*.out |sort |grep -v "#"> statagg.out
awk '{print $1}' stat.out |uniq -c |awk '{print $2,$1}' |sort> estCount.out
join estCount.out statagg.out > estCoverage.out
overlapSelect ../all_mrnaFiltered.psl.gz ../$TABLE.bed -statsOutput stdout |sort > mrna.out
zcat ../splicedEst.psl.gz |cut -f 10 |sort |uniq> splicedEst.id
overlapSelect ../splicedEst.psl.gz ../$TABLE.bed -statsOutput stdout | sort > splicedstat.out
awk '{print $1}' splicedstat.out |uniq -c |awk '{print $2,$1}' |sort> splicedEstCount.out
overlapSelect ../splicedEst.psl.gz ../$TABLE.bed -statsOutput -aggregate stdout | sort > splicedagg.out
join splicedEstCount.out splicedagg.out > splicedEstCoverage.out

awk '$3>0.50{print $1}' mrna.out |sort > mrna.id
awk '$2>=10 && $3>0.50{print $0}' estCoverage.out > est10.out
awk '$2>=5 && $3>0.50{print $0}' estCoverage.out > est5.out
awk '$2>=1 && $3>0.50{print $0}' estCoverage.out > est.out
awk '$2>=5 && $3>0.50{print $0}' splicedEstCoverage.out > splicedEst5.out
awk '$2>=10 && $3>0.50{print $0}' splicedEstCoverage.out > splicedEst10.out
awk '{print $1}' est.out |sort > est1.id
awk '{print $1}' est5.out |sort > est5.id
awk '{print $1}' est10.out |sort > est10.id
awk '{print $1}' splicedEst5.out |sort > splicedEst5.id
awk '{print $1}' splicedEst10.out |sort > splicedEst10.id

join est10.id mrna.id > est10Mrna.id
join est5.id mrna.id > est5Mrna.id
grep -F -f est5.id ../$TABLE.bed >  pseudoEst5.bed
cat mrna.id est10.id |sort |uniq > mrnaEst10.id
cat mrna.id est5.id |sort |uniq > mrnaEst5.id
cat mrna.id est1.id |sort |uniq > mrnaEst1.id
$SCRIPT/selectById -tsv 1 est5Mrna.id 4 ../$TABLE.bed > pseudoEst5AndMrna.bed
$SCRIPT/selectById -tsv 1 est10Mrna.id 4 ../$TABLE.bed > pseudoEst10AndMrna.bed
$SCRIPT/selectById -tsv 1 mrnaEst1.id 4 ../$TABLE.bed > pseudoEstMrna.filter.bed
$SCRIPT/selectById -tsv 1 mrnaEst5.id 4 ../$TABLE.bed > pseudoEst5Mrna.bed
$SCRIPT/selectById -tsv 1 mrnaEst10.id 4 ../$TABLE.bed >  pseudoEst10Mrna.bed
$SCRIPT/selectById -tsv 1 est5.id 4 ../$TABLE.bed > pseudoEst5.bed
$SCRIPT/selectById -tsv 1 est10.id 4 ../$TABLE.bed > pseudoEst10.bed
overlapSelect pseudoEst10Mrna.bed pseudo${GENE2}Cds.bed  est10MrnaRefSeq.bed

overlapSelect retroKgCoding.bed pseudoEst5Mrna.bed pseudoEst5MrnaNotKg.bed -nonOverlapping
overlapSelect retroKgCoding.bed pseudoEst10Mrna.bed pseudoEst10MrnaNotKg.bed -nonOverlapping
overlapSelect retroKgCoding.bed pseudoEst5AndMrna.bed pseudoEst5AndMrnaNotKg.bed -nonOverlapping

#sort -k2,2nr estCoverage.out |head -50|awk '{print $1}' > top50
#echo $SCRIPT/selectById -tsv 1 top50 4 ../$TABLE.bed   top50.bed
#$SCRIPT/selectById -tsv 1 top50 4 ../$TABLE.bed  > top50.bed

for i in ucscRetroProtein ucscRetroTranscript ; do $SPLITBYAGE ${i}.bed ${i}.ancient.bed ${i}.recent.bed ; done
echo "expressed retro stats"
wc -l *.id |sort -n
overlapSelect ../${GENE2}.tab.gz ../$TABLE.bed pseudo${GENE2}.bed -selectFmt=genePred

#regenerate gene predictions from expressed retros

#set with 1 or more EST (tight) not shuffle events
#look for orf in window 200 bp around retro
tawk '$2>200{$2=$2-200; $3=$3+200;print $0}$2<=200{$3=$3+200;print $0}' pseudoEstAllNotShuffle.bed > pseudoEstAll.bed
echo "orfBatch $DB pseudoEstAll.bed pseudoEstAll.out pseudoEstAll.gp to borfEst.out"
orfBatch $DB pseudoEstAll.bed pseudoEstAll.out pseudoEstAll.gp > borfEst.out
echo gene-check  -nib-dir $NIB pseudoEstAll.gp checkEst.rdb
~markd/compbio/genefinding/GeneTools/bin/x86_64/opt/gene-check  -genome-seqs $NIB pseudoEstAll.gp checkEst.rdb
awk '$7=="ok"&& $25=="noStart" || $25==""{print $1}' checkEst.rdb > goodOrf.list

$SCRIPT/selectById -tsv 1 goodOrf.list 4 pseudoEstAll.out > pseudoEstAll.good.out

echo "histogram of best.orf scores"
awk '{print $52}' pseudoEstAll.good.out |textHistogram stdin -maxBinCount=40 -binSize=10
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed weak\", posConf = \""$52"\" where name = \""$4"\" ;"}' pseudoEstAll.good.out > updateExp.sql
echo "hgsql $DB updateExp.sql"
hgsql $DB < updateExp.sql

#set with 5 Est and 1 spliced mRNA
overlapSelect pseudo${GENE1}Cds.bed pseudoEst5AndMrna.bed pseudoEst5AndMrnaNotShuffle.bed -nonOverlapping
tawk '$2>200{$2=$2-200; $3=$3+200;print $0}$2<=200{$3=$3+200;print $0}' pseudoEst5AndMrnaNotShuffle.bed > pseudoEst5AndMrna.window.bed
orfBatch $DB pseudoEst5AndMrna.window.bed pseudoEst5AndMrna.out pseudoEst5AndMrna.gp > borfEst5AndMrna.out
echo gene-check  -nib-dir $NIB pseudoEst5AndMrna.gp checkEst5AndMrna.rdb
~markd/compbio/genefinding/GeneTools/bin/x86_64/opt/gene-check  -genome-seqs $NIB pseudoEst5AndMrna.gp checkEst5AndMrna.rdb
tawk '$7=="ok" && $25=="noStart" || $25==""{print $1}' checkEst5AndMrna.rdb > goodOrf5AndMrna.list
$SCRIPT/selectById -tsv 1 goodOrf5AndMrna.list 4 pseudoEst5AndMrna.out > pseudoEst5AndMrna.good.out
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed strong\", posConf = \""$52"\" where name = \""$4"\" ;"}' pseudoEst5AndMrna.good.out > updateExp5.sql
echo "hgsql $DB  updateExp5.sql"
hgsql $DB < updateExp5.sql

#bad orfs
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' checkEst.rdb > badOrf.list
awk '{print "update '$TABLE' set type = \"expressed weak noOrf\" where name = \""$1"\" ;"}' badOrf.list> updateBad.sql
hgsql $DB < updateBad.sql
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' checkEst5AndMrna.rdb > badOrf5.list
awk '{print "update '$TABLE' set type = \"expressed strong noOrf\" where name = \""$1"\" ;"}' badOrf5.list> updateBad5.sql
hgsql $DB < updateBad5.sql
wc -l update*.sql

$SCRIPT/selectById -tsv 1 goodOrf.list 1 pseudoEstAll.gp > pseudoExpressed.gp
echo "ldHgGene $DB pseudoExpressed pseudoExpressed.gp -genePredExt -predTab"
ldHgGene $DB ucscRetroExpressed pseudoExpressed.gp -genePredExt -predTab
genePredToBed pseudoExpressed.gp > pseudoTmp.bed
$SCRIPT/selectById -tsv 1 goodOrf.list 4 pseudoEstAll.bed > pseudoExpressed.bed
echo "length histogram of coding region"
awk '{print $7-$6}' pseudoExpressed.gp|textHistogram stdin -maxBinCount=100 -binSize=50

tawk '($8-$7)>300{print }' pseudoTmp.bed>  pseudoTmp2.bed
$SCRIPT/selectById -tsv 4 pseudoTmp2.bed 4 ../$TABLE.bed > pseudoEst100AA.bed
echo tawk '$8-$7>300{print }' pseudoExpressed.bed redirect pseudoEst100AA.bed
$SCRIPT/selectById -tsv 4 pseudoEst5AndMrna.good.out 4 ../$TABLE.bed > pseudo5Est100AA.bed
wc -l pseudoExpressed.bed pseudoEst100AA.bed goodOrf.list pseudoEst5AndMrna.bed pseudoEst5AndMrna.good.out
overlapSelect pseudoEst5AndMrna.bed pseudo5Est100AA.bed pseudoEst5AndMrna100AA.bed

#split retros by age

for bed in pseudo${GENE2} pseudoEstMrna.filter pseudoEstAll pseudoExpressed pseudoEst5Mrna pseudoEst5 pseudoEst5AndMrna pseudoEst100AA pseudo5Est100AA pseudoEst5AndMrna100AA pseudo${GENE1}Cds pseudo${GENE1}Cds50 pseudo${GENE2}Cds pseudo${GENE2}Cds50 pseudoEst5MrnaNotKg pseudoEst10MrnaNotKg pseudoEst5AndMrnaNotKg; do $SPLITBYAGE ${bed}.bed ${bed}.ancient.bed ${bed}.recent.bed; done

mkdir -p age
cd age
$SCRIPT/doAge ../../DEF pseudoEst5Mrna
$SCRIPT/doAge ../../DEF retroMrnaInfoLessZnf 
$SCRIPT/doAge ../../DEF retroKgCoding 

#retros involved in alt-splicing

if [[ -n $ALTSPLICE ]] 
then
echo "hgsql $DB -N -B -e select * from $ALTSPLICE " 
hgsql $DB -N -B -e "select * from $ALTSPLICE " |cut -f2-19 | agxToBed stdin altSplice.bed
fi
if [ -f altSplice.bed ] 
then
 echo "overlapSelect altSplice.bed ../$TABLE.bed retroSplice.bed"
 overlapSelect altSplice.bed ../$TABLE.bed retroSplice.bed
 $SPLITBYAGE retroSplice.bed retroSplice.ancient.bed retroSplice.recent.bed 
fi

echo '-------- END script analyseExpress.sh -------------------'
echo '------ run script ucscRetroStep6.sh to make Html pages -------------------'
