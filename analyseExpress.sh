#!/bin/bash
set -beEu -o pipefail
source $1
echo '-------- script analyseExpress.sh -------------------'
hgsql $PDB -N -B -e "select * from hgncXref " > hgncXref.tab
overlapSelect ../$GENE2.multiCds.bed ../$TABLE.bed pseudoRefGeneCds.bed
overlapSelect ../$GENE2.multiCds.bed ../retroMrnaInfo.12.bed -statsOutput pseudoRefGeneCds.out
overlapSelect ../$GENE2.multiCds.bed ../retroMrnaInfo.12.bed pseudoRefGeneCds50.bed -overlapThreshold=0.50
overlapSelect -selectFmt=genePred ../$GENE2.tab.gz pseudoRefGeneCds.bed shuffleEns.bed
overlapSelect -selectFmt=genePred ../$GENE2.multiCDSExon.genePred pseudoRefGeneCds.bed shuffleEnsMulti.bed
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
awk '$2>=5 && $3>0.50{print $0}' splicedEstCoverage.out > splicedEst5.out
awk '$2>=10 && $3>0.50{print $0}' splicedEstCoverage.out > splicedEst10.out
awk '{print $1}' est5.out |sort > est5.id
awk '{print $1}' est10.out |sort > est10.id
awk '{print $1}' splicedEst5.out |sort > splicedEst5.id
awk '{print $1}' splicedEst10.out |sort > splicedEst10.id

join est10.id mrna.id > est10Mrna.id
join est5.id mrna.id > est5Mrna.id
grep -F -f est5.id ../$TABLE.bed >  pseudoEst5.bed
cat mrna.id est10.id |sort |uniq > mrnaEst10.id
cat mrna.id est5.id |sort |uniq > mrnaEst5.id
$SCRIPT/selectById -tsv 1 est5Mrna.id 4 ../$TABLE.bed > pseudoEst5AndMrna.bed
$SCRIPT/selectById -tsv 1 est10Mrna.id 4 ../$TABLE.bed > pseudoEst10AndMrna.bed
$SCRIPT/selectById -tsv 1 mrnaEst5.id 4 ../$TABLE.bed > pseudoEstMrna.filter.bed
$SCRIPT/selectById -tsv 1 mrnaEst5.id 4 ../$TABLE.bed > pseudoEst5Mrna.bed
$SCRIPT/selectById -tsv 1 mrnaEst10.id 4 ../$TABLE.bed >  pseudoEst10Mrna.bed
$SCRIPT/selectById -tsv 1 est5.id 4 ../$TABLE.bed > pseudoEst5.bed
$SCRIPT/selectById -tsv 1 est10.id 4 ../$TABLE.bed > pseudoEst10.bed
overlapSelect pseudoEst10Mrna.bed pseudoRefGeneCds.bed  est10MrnaRefSeq.bed
echo "expressed retro stats"
wc -l *.id |sort -n
overlapSelect ../refGene.tab.gz ../$TABLE.bed pseudoRefGene.bed -selectFmt=genePred

#regenerate gene predictions from expressed retros
wc -l ../pseudoEstAll.bed
tawk '$5 > 600 && $2 > 25{$2=$2-25; $3=$3+25;print $0}$5 > 600 && $2 <= 25{ $3=$3+25;print $0}' pseudoEstOrMrna.bed > pseudoEstOrMrna600.bed
echo "orfBatch $DB pseudoEstOrMrna600.bed pseudoEstOrMrna600.out pseudoEstOrMrna600.gp to borf.out"
/cluster/home/baertsch/bin/i386/orfBatch $DB pseudoEstOrMrna600.bed pseudoEstOrMrna600.out pseudoEstOrMrna600.gp >borf.out
echo "genePredSingleCover pseudoEstOrMrna600.gp pseudoEstOrMrna600.single.gp"
genePredSingleCover pseudoEstOrMrna600.gp pseudoEstOrMrna600.single.gp
awk '{print "mrna."$1}' pseudoEstOrMrna600.single.gp |sort > pseudoEstOrMrna600.ids
echo "predicted orfs from retros"
wc -l pseudoEstOrMrna600.single.gp


#tawk '{print $3-$2}' ../pseudoGeneLinkNoOverlapFilter.bed > lengthAll.txt
#tawk '{print $4"."$1"."$2,$0}' ../pseudoGeneLinkNoOverlapFilter.bed > pseudoJoin.bed
#tawk '$14 == "expressed"{print $3-$2}' ../pseudoGeneLinkNoOverlapFilter.bed > lengthExp.txt
#tawk '$14 != "expressed"{print $3-$2}' ../pseudoGeneLinkNoOverlapFilter.bed > lengthNonExp.txt
#tawk '$14 == "expressed"{print $1,$2,$3,$4"."$1"."$2,$5,$6}' ../pseudoGeneLinkNoOverlapFilter.bed | sort -k5,5nr >pseudoExpressed.bed
#tawk '$14 != "expressed"{print $1,$2,$3,$4"."$1"."$2,$5,$6}' ../pseudoGeneLinkNoOverlapFilter.bed | sort -k5,5nr >pseudoNotExpressed.bed
#overlapSelect /san/sanvol1/scratch/pseudo/hg18/sortedKnownGene.tab.gz pseudoExpressed.bed -selectFmt=genePred -inFmt=bed exp.bed -nonOverlapping -selectCds
#tawk '{$2=$2-2000;$3=$3+2000;print $0}' exp.bed > expWindow.bed
 
#tawk 'BEGIN{print "name","Parent","score","child","parent","Spliced","knownGene"}' > exp.label
#tawk '$14 == "expressed"{print $4"."$1"."$2,$41,$5,$1,$16,$24"/"$20,$47}' ../pseudoGeneLinkNoOverlapFilter.bed | sort -k5,5nr >>exp.label
#~markd/bin/bedToHtmlDir -no-sort -dir-frame-per 25 -title "expressed retroGenes not overlapping refSeq" -label-tsv exp.label -browser-url http://hgwdev-baertsch.cse.ucsc.edu hg18 exp.bed ~/.html/expressed
#tawk '$14=="expressed"{print $0}' ../pseudoGeneLinkNoOverlapFilter.bed > pseudoGeneLinkExpressed.bed

#orfBatch hg18 pseudoGeneLinkExpressed.bed expOut.bed pseudoExpressed.genePred > borf.tab
#gene-check  -nib-dir $NIB -incl-ok pseudoExpressed.genePred check.rdb 
#sort check.rdb > check.rdb.sort
#awk '{print "update pseudoGeneLink set thickStart = "$7", thickEnd = "$8" , posConf = "$55" where chrom = \""$1"\" and chromStart = \""$2"\" and chromEnd = \""$3"\";"}' expOut.bed >updateExp.sql
##hgsql hg18 -N -B < deleteExp.sql
#tawk 'BEGIN{print "name","mRNA","Parent","score","child","parent","Spliced","knownGene","borf score","strand","orf st","orf","problem"}' > expOrf.label
#tawk '{print $4"."$1"."$2,$4,$41,$5,$1,$16,$24"/"$20,$47,$55,$6}' expOut.bed | sort | \
#  join - check.rdb.sort -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.5,2.6,2.13 | \
#  awk '{OFS="\t";$1=$1;print $0}' | sort -k4,4nr >>expOrf.label
#~markd/bin/bedToHtmlDir -dir-below -no-sort -dir-frame-per 25 -title "expressed retroGenes " -label-tsv expOrf.label -browser-url http://hgwdev-baertsch.cse.ucsc.edu hg18 expWindow.bed ~/.html/expressedOrf
#overlapSelect notMacaca.bed pseudoExpressed.bed pseudoNotMacaca.bed -overlapThreshold=0.80
#~markd/bin/bedToHtmlDir -no-sort -title "Expressed retroGene not in macaca" -browser-url http://genome-test.cse.ucsc.edu hg18 pseudoNotMacaca.bed ~/.html/notMacaca

#HUMAN SPECIFIC
#overlapSelect chimpDels.bed pseudoExpressed.bed retroHumExp.bed
#overlapSelect chimpDels.bed pseudoNotExpressed.bed retroHumNotExp.bed
#tawk '{print $1,$2,$3}' retroHumExp.bed > retroHumExp.id
#grep -F -f retroHumExp.id pseudoExpressedLong.bed > retroHumLongExp.bed  
#overlapSelect indelible.txt pseudoExpressed.bed retroHumanKrish.bed -selectFmt=bed
#tawk '{print $1,$2,$3}' retroHumanKrish.bed > retroHumanKrish.id
#grep -F -f retroHumanKrish.id pseudoExpressedLong.bed > retroHumLongKrish.bed  

#MACACA SPECIFIC
#for i in `ls /cluster/bluearc/macaca/axtRecipBest3/*.axt` ;do b=`basename $i` ; echo $b; axtToBed $i ${b%%.axt}_macaca.bed ; done   

#HIGH CONF PRIMATE RETRO
#cat ../pseudoGeneLinkNoOverlapFilter.bed  | \
#    tawk '$35 > 700 && $5 > 425 && $14 == "expressed" && $24 > 1{print $0}'  | sort -k5,5n > retroPrimateExpressed.bed
#orfBatch hg18 retroPrimateExpressed.bed expPrimate.bed retroPrimateExpressed.genePred > borfPrimate.tab
#echo gene-check  -nib-dir $NIB retroPrimateExpressed.genePred checkPrimate.rdb 
#gene-check  -nib-dir $NIB retroPrimateExpressed.genePred checkPrimate.rdb 
#sort checkPrimate.rdb > checkPrimate.rdb.sort
#tawk 'BEGIN{print "name","mRNA","Parent","score","child","parent","Spliced","knownGene","borf score","strand","orf st","orf","problem"}' > retroPrimateExpressed.label
#tawk '{print $4"."$1"."$2,$4,$41,$5,$1,$16,$24"/"$20,$47,$55,$6}' retroPrimateExpressed.bed | sort | \
#  join - checkPrimate.rdb.sort -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.5,2.6,2.13 | \
#  awk '{OFS="\t";$1=$1;print $0}' | sort -k4,4nr >>retroPrimateExpressed.label 
#

##### new set

# look for orf in bigger window around retro

#echo "ldHgGene $DB pseudoExpressed pseudoExpressed.gp -genePredExt -predTab"
#ldHgGene $DB pseudoExpressed pseudoExpressed.gp -genePredExt -predTab

#set with 1 or more EST (tight)
tawk '$2>50{$2=$2-50; $3=$3+50;print $0}$2<=50{$3=$3+50;print $0}' ../pseudoEstAll.bed > pseudoEstAll.bed
echo "orfBatch $DB pseudoEstAll.bed pseudoEstAll.out pseudoEstAll.gp to borfEst.out"
orfBatch $DB pseudoEstAll.bed pseudoEstAll.out pseudoEstAll.gp > borfEst.out
echo gene-check  -nib-dir $NIB pseudoEstAll.gp checkEst.rdb
~markd/compbio/genefinding/GeneTools/bin/x86_64/opt/gene-check  -genome-seqs $NIB pseudoEstAll.gp checkEst.rdb
awk '$7=="ok"&& $25=="noStart" || $25==""{print $1}' checkEst.rdb > goodOrf.list

#fgrep -f goodOrf.list pseudoEstAll.out > pseudoEstAll.good.out
$SCRIPT/selectById -tsv 1 goodOrf.list 4 pseudoEstAll.out > pseudoEstAll.good.out

echo "histogram of best.orf scores"
awk '{print $55}' pseudoEstAll.good.out |textHistogram stdin -maxBinCount=40 -binSize=10
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed weak\", posConf = \""$55"\" where name = \""$4"\" ;"}' pseudoEstAll.good.out > updateExp.sql
head updateExp.sql
echo "hgsql $DB updateExp.sql"
hgsql $DB < updateExp.sql
#set with 5Est and 1spliced mRNA
tawk '$2>50{$2=$2-50; $3=$3+50;print $0}$2<=50{$3=$3+50;print $0}' pseudoEst5AndMrna.bed > pseudoEst5AndMrna.window.bed
orfBatch $DB pseudoEst5AndMrna.window.bed pseudoEst5AndMrna.out pseudoEst5AndMrna.gp > borfEst5AndMrna.out
echo gene-check  -nib-dir $NIB pseudoEst5AndMrna.gp checkEst5AndMrna.rdb
~markd/compbio/genefinding/GeneTools/bin/x86_64/opt/gene-check  -genome-seqs $NIB pseudoEst5AndMrna.gp checkEst5AndMrna.rdb
tawk '$7=="ok" && $25=="noStart" || $25==""{print $1}' checkEst5AndMrna.rdb > goodOrf5AndMrna.list
#tawk '$6=="ok"&& $7=="ok"{print $1}' checkEst5AndMrna.rdb > goodOrf5AndMrna.list
#fgrep -f goodOrf5AndMrna.list pseudoEst5AndMrna.out > pseudoEst5AndMrna.good.out
$SCRIPT/selectById -tsv 1 goodOrf5AndMrna.list 4 pseudoEst5AndMrna.out > pseudoEst5AndMrna.good.out
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed strong\", posConf = \""$55"\" where name = \""$4"\" ;"}' pseudoEst5AndMrna.good.out > updateExp5.sql
head updateExp5.sql
echo "hgsql $DB  updateExp5.sql"
hgsql $DB < updateExp5.sql
#bad orfs
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' checkEst.rdb > badOrf.list
awk '{print "update '$TABLE' set type = \"expressed weak noOrf\" where name = \""$1"\" ;"}' badOrf.list> updateBad.sql
hgsql $DB < updateBad.sql
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' checkEst5AndMrna.rdb > badOrf5.list
awk '{print "update '$TABLE' set type = \"expressed strong noOrf\" where name = \""$1"\" ;"}' badOrf5.list> updateBad5.sql
hgsql $DB < updateBad5.sql

#fgrep -f goodOrf.list pseudoEstAll.gp > pseudoExpressed.gp
$SCRIPT/selectById -tsv 1 goodOrf.list 1 pseudoEstAll.gp > pseudoExpressed.gp
echo "ldHgGene $DB pseudoExpressed pseudoExpressed.gp -genePredExt -predTab"
ldHgGene $DB ucscRetroExpressed pseudoExpressed.gp -genePredExt -predTab
genePredToBed pseudoExpressed.gp > pseudoTmp.bed
$SCRIPT/selectById -tsv 1 goodOrf.list 4 ../$TABLE.bed > pseudoExpressed.bed
echo "length histogram of coding region"
awk '{print $7-$6}' pseudoExpressed.gp|textHistogram stdin -maxBinCount=100 -binSize=50

tawk '($8-$7)>300{print }' pseudoTmp.bed>  pseudoTmp2.bed
$SCRIPT/selectById -tsv 4 pseudoTmp2.bed 4 ../$TABLE.bed > pseudoEst100AA.bed
echo tawk '$8-$7>300{print }' pseudoExpressed.bed redirect pseudoEst100AA.bed
$SCRIPT/selectById -tsv 4 pseudoEst5AndMrna.good.out 4 ../$TABLE.bed > pseudo5Est100AA.bed
wc -l pseudoExpressed.bed pseudoEst100AA.bed goodOrf.list pseudoEst5AndMrna.bed pseudoEst5AndMrna.good.out
overlapSelect pseudoEst5AndMrna.bed pseudo5Est100AA.bed pseudoEst5AndMrna100AA.bed

#split retros by age

for bed in pseudoRefGene pseudoEstMrna.filter pseudoEstAll pseudoExpressed pseudoEst5Mrna pseudoEst5 pseudoEst5AndMrna pseudoEst100AA pseudo5Est100AA pseudoEst5AndMrna100AA pseudoRefGeneCds pseudoRefGeneCds50 shuffleEns shuffleEnsMulti ; do $SPLITBYAGE ${bed}.bed ${bed}.ancient.bed ${bed}.recent.bed; done

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
echo '-------- run script ucscRetroStep6.sh to make Html pages -------------------'
