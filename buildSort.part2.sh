# run from kkstore02
DB=$2
BASE=$3
RESULT=$1
pushd $1
cd ..
pwd
echo cat chr*pseudoNoOverlap.bed to $BASE/pseudoGeneLinkNoOverlap.bed
cat chr*pseudoNoOverlap.bed > $BASE/pseudoGeneLinkNoOverlap.bed
popd
tawk '$5>=350{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkNoOverlapFilter.bed                       
tawk '$5>=425{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLink425.bed
tawk '$5>=510{print $0}' pseudoGeneLinkNoOverlap.bed > retroMrnaInfo.raw.bed

echo creating retroMrnaAli.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkSelect.tab
pslSelect -qtStart=pseudoGeneLinkSelect.tab pseudo.psl retroMrnaAli.psl
wc -l retroMrnaAli.psl pseudoGeneLinkSelect.tab
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
fgrep -w -f bothZnf.lst sortedKnownGene.tab > kgZnf.gp
grep -v -F -f bothZnf.lst retroMrnaInfo.raw.bed > retroMrnaInfoLessZnf.bed
cp retroMrnaInfoLessZnf.bed retroMrnaInfo.bed
wc -l retroMrnaInfo.raw.bed pseudoGeneLink425.bed retroMrnaInfo.bed
#echo Loading Bed
#hgLoadBed $DB pseudoGeneLink pseudoGeneLinkNoOverlap.bed -hasBin -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/pseudoGeneLink.sql
##hgLoadBed $DB pseudoGeneLink pseudoGeneLink425.bed -hasBin -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/pseudoGeneLink.sql
##echo Loading Psl
##hgLoadPsl $DB pseudoMrna2.psl
hgLoadBed $DB -noBin retroMrnaInfo -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/retroMrnaInfo.sql retroMrnaInfo.bed
hgLoadPsl $DB retroMrnaAli.psl

#make exon list
hgsql $DB -N -B < refGeneMultiExon.sql |cut -f2-12 > $RESULT/refGeneMultiExon.genePred
hgsql $DB -N -B < kgMultiExon.sql > $RESULT/kgMultiExon.genePred
hgsql $DB -N -B < mgcMultiExon.sql > $RESULT/mgcMultiExon.genePred
genePredFilter refGeneMultiExon.genePred refGeneMultiCDSExon.genePred -cdsExons=2
genePredFilter kgMultiExon.genePred kgMultiCDSExon.genePred -cdsExons=2
genePredFilter mgcMultiExon.genePred mgcMultiCDSExon.genePred -cdsExons=2
cat kgMultiCDSExon.genePred mgcMultiCDSExon.genePred refGeneMultiCDSExon.genePred > kgMgcRefGeneMultiCds.gp

ldHgGene $DB rbRefGeneMulti refGeneMultiCDSExon.genePred -predTab
featureBits $DB refGene:cds -bed=refGeneCds.bed
featureBits $DB rbRefGeneMulti;
#52619505 bases of 2881515245 (1.826%) in intersection
featureBits $DB refGene;     
#54738314 bases of 2881515245 (1.900%) in intersection
featureBits $DB rbRefGeneMulti refGeneCds.bed -bed=refGeneMultiCds.bed 
#29845142 bases of 2881515245 (1.036%) in intersection

#make extract extra columns for html pages
makeRetroExtraAttr.sh hg18

#check for expression
#zcat all_mrna.psl.gz |sort -k10,10> all_mrna.qName.psl
#for i in `chromsNoY` chrY ; do echo $i ; hgsql $DB -N -B -e "select * from ${i}_est" |cut -f2-22 >> est.psl ; done
#nohup nice pslCDnaGenomeMatch all_mrna.qName.psl S1.len mrnaNoversion.2bit /cluster/data/$DB/nib all_mrnaFiltered.psl -score=mrnaMatch.tab -bedOut=mrnaMis.bed  
#zcat est.psl |sort -k10,10> est.qName.psl
#pslSplit nohead /san/sanvol1/scratch/est/hg18 est.qName.psl -chunkSize=12000
#mkdir run.est
#
#cd run.est
#cat <! gsub
#LOOP
#pslCDnaGenomeMatch $(path1) chrom.sizes /san/sanvol1/scratch/est/hg18/est.2bit /scratch/hg/hg18/nib $(dir1)/../output/$(file1).filter.psl -score=$(dir1)/../log/$(file1).mrnaMatch.tab -bedOut=$(dir1)/../log/$(file1).mrnaMis.bed -minDiff=5
#ENDLOOP
#ssh pk -T "cd /cluster/data/hg18/bed/pseudo/run.est ; para make spec"
#cd ..
#cat /san/sanvol1/scratch/est/hg18/output/* | sort -k14,14 -k16,16n >~/hg18/pseudo/estFiltered.psl
#awk -f splitTarget.awk < estFiltered.psl 
#for i in `ls chr*.psl` ; do echo $i ; mv $i ${i%%.psl}_estFiltered.psl ; done
#for i in `ls chr*.psl` ; do hgLoadPsl $DB $i ; done
overlapSelect estFiltered.psl.gz retroMrnaInfoLessZnf.bed pseudoEstAll.bed
grep -F -f bothZnf.lst retroMrnaInfo.bed > retroZnf.bed

mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid > all_mrna.log 2> all_mrna.err
mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 -utrMergeSize=10 all_mrna_utr.gp -cdsDb=$DB -keepInvalid > all_mrna_utr.log 2> all_mrna_utr.err

awk '$8>1{print }' all_mrna.gp > all_mrna_multiExon.gp
awk '$8>1{print }' all_mrna_utr.gp > all_mrna_multiExonUTR.gp

cd stats
doStat
cd ..
echo pwd
pwd
analyseExpress.sh $DB
#nice -4 zcat pseudoMrnaLink*.txt.gz | awk '$28 == 1 && $23 < 50 {print $0}' | 
#awk 'NF==47{$37=$37" -1 0";print $0}NF==48{$37=$37" 0"; print $0}' pseudoGeneLinkSort.bed > pseudoGeneLinkSort.bed
#awk 'NF==47{$37="0 0";print $0}NF==48{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkNoOverlap.bed
# awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
echo pwd
pwd
echo 'BUILD frame - catting axts'
pushd $RESULT/axt
zcat *.axt.gz > all.axt
popd
mkdir -p frames/axt
mkdir -p frames/mafRaw
echo cat $RESULT/axt/all.axt pipe axtSplitByTarget stdin frames/axt
pwd
cat $RESULT/axt/all.axt | axtSplitByTarget stdin frames/axt
cd frames/axt
awk '$2~/chr*/{print $5}' chr*.axt | sort -u > qlist
awk '{i=index($1,"-");print substr($1,1,i-1),$1}' qlist |sort > qlist.ext
pwd
echo sort ../../final.len > mrna.sort.len
wc -l ../../final.len
sort ../../final.len > mrna.sort.len
join mrna.sort.len qlist.ext|awk '{OFS="\t";print $3,$2}' > mrnaExt.len 

for i in `ls *.axt` ; do echo $i ; axtToMaf -tPrefix=${DB}. -qPrefix=mrna. $i S1.len mrnaExt.len ../mafRaw/${i%%.axt}.maf ; done

