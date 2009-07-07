## called from ~baertsch/baertsch/scripts/buildSort.sh
RESULT=$1
DB=$2
BASE=$3
SCRIPT=$4
VERSION=4
TABLE=retroMrnaInfo$VERSION
echo "---------------------------"
echo "Starting buildSort.part2.sh"
echo "---------------------------"
date
cd $BASE
pwd
echo cat chr*_NoOverlap.bed to pseudoGeneLinkNoOverlap.bed
cat chr*_NoOverlap.bed > pseudoGeneLinkNoOverlap.bed
wc -l chr*_NoOverlap.bed
tawk '$5>=350{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkNoOverlapFilter.bed                       
tawk '$5>=425{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLink425.bed
tawk '$5>=510{print $0}' pseudoGeneLinkNoOverlap.bed > retroMrnaInfo.raw.bed
tawk '$5>=650{print $0}' pseudoGeneLinkNoOverlap.bed > retroMrnaInfo650.bed
cut -f 1-12 retroMrnaInfo650.bed > retroMrnaInfo.12.bed
wc -l pseudoGeneLinkNoOverlap.bed pseudoGeneLinkNoOverlapFilter.bed pseudoGeneLink425.bed retroMrnaInfo.raw.bed retroMrnaInfo650.bed retroMrnaInfo.12.bed

echo creating retroMrnaAli.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkSelect.tab
pslSelect -qtStart=pseudoGeneLinkSelect.tab pseudo.psl retroMrnaAli$VERSION.psl
wc -l retroMrnaAli$VERSION.psl pseudoGeneLinkSelect.tab
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
hgsql $DB -N -B -e "select  k.name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from knownGene k, knownToPfam p \
            where k.name = p.name and p.value in (
        'PF00096',
        'PF01352',
        'PF06758',
        'PF00047',
        'PF07654'  
        );" > zincKg.gp
cut -f 1 zincKg.gp |sort | uniq > zincKg.lst

zcat sortedKnownGene.tab.gz > sortedKnownGene.tab
#fgrep -w -f bothZnf.lst sortedKnownGene.tab > kgZnf.gp
#grep -v -F -f bothZnf.lst retroMrnaInfo.raw.bed > retroMrnaInfoLessZnf.bed
fgrep -w -f zincKg.lst sortedKnownGene.tab > kgZnf.gp
grep -v -F -f zincKg.lst retroMrnaInfo.raw.bed > retroMrnaInfoLessZnf.bed
cp retroMrnaInfoLessZnf.bed retroMrnaInfo$VERSION.bed
textHistogram -col=5 retroMrnaInfo4.bed -binSize=50 -maxBinCount=50
wc -l retroMrnaInfo.raw.bed pseudoGeneLink425.bed retroMrnaInfo$VERSION.bed
#echo Loading Bed
#hgLoadBed $DB pseudoGeneLink pseudoGeneLinkNoOverlap.bed -hasBin -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/pseudoGeneLink.sql
##hgLoadBed $DB pseudoGeneLink pseudoGeneLink425.bed -hasBin -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/pseudoGeneLink.sql
##echo Loading Psl
##hgLoadPsl $DB pseudoMrna2.psl
hgLoadBed $DB -noBin retroMrnaInfoXX -sqlTable=/cluster/home/baertsch/kent/src/hg/lib/retroMrnaInfo.sql retroMrnaInfo$VERSION.bed
hgsql $DB -e "drop table $TABLE;"
hgsql $DB -e "alter table retroMrnaInfoXX rename $TABLE;"
hgLoadPsl $DB retroMrnaAli$VERSION.psl

#make exon list
cd $BASE
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from refGene where exonCount > 1" > refGeneMultiExon.genePred
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from knownGene where exonCount > 1" > kgMultiExon.genePred
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from mgcGenes where exonCount > 1" > mgcMultiExon.genePred
echo "genePredFilter refGeneMultiExon.genePred refGeneMultiCDSExon.genePred -cdsExons=2"
genePredFilter refGeneMultiExon.genePred refGeneMultiCDSExon.genePred -cdsExons=2
echo "genePredFilter kgMultiExon.genePred kgMultiCDSExon.genePred -cdsExons=2"
genePredFilter kgMultiExon.genePred kgMultiCDSExon.genePred -cdsExons=2
echo "genePredFilter mgcMultiExon.genePred mgcMultiCDSExon.genePred -cdsExons=2"
genePredFilter mgcMultiExon.genePred mgcMultiCDSExon.genePred -cdsExons=2
wc -l kgMultiCDSExon.genePred refGeneMultiCDSExon.genePred mgcMultiCDSExon.genePred
echo "cat kgMultiCDSExon.genePred mgcMultiCDSExon.genePred refGeneMultiCDSExon.genePred to  kgMgcRefGeneMultiCds.gp"
cat kgMultiCDSExon.genePred mgcMultiCDSExon.genePred refGeneMultiCDSExon.genePred > kgMgcRefGeneMultiCds.gp

echo "ldHgGene $DB rbRefGeneMulti refGeneMultiCDSExon.genePred -predTab"
ldHgGene $DB rbRefGeneMulti refGeneMultiCDSExon.genePred -predTab
featureBits $DB refGene:cds -bed=refGeneCds.bed
featureBits $DB rbRefGeneMulti;
#52619505 bases of 2881515245 (1.826%) in intersection
featureBits $DB refGene;     
#54738314 bases of 2881515245 (1.900%) in intersection
featureBits $DB rbRefGeneMulti refGeneCds.bed -bed=refGeneMultiCds.bed 
#29845142 bases of 2881515245 (1.036%) in intersection

#make extract extra columns for html pages
$SCRIPT/makeRetroExtraAttr.sh $DB $TABLE

#check for expression
#zcat all_mrna.psl.gz |sort -k10,10> all_mrna.qName.psl
#for i in `chromsNoY` chrY ; do echo $i ; hgsql $DB -N -B -e "select * from ${i}_est" |cut -f2-22 >> est.psl ; done
#nohup nice pslCDnaGenomeMatch all_mrna.qName.psl S1.len mrnaNoversion.2bit /cluster/data/$DB/nib all_mrnaFiltered.psl -score=mrnaMatch.tab -bedOut=mrnaMis.bed  
#zcat est.psl |sort -k10,10> est.qName.psl
#pslSplit nohead /san/sanvol1/scratch/est/$DB est.qName.psl -chunkSize=12000
#mkdir run.est
#
#cd run.est
#cat <! gsub
#LOOP
#pslCDnaGenomeMatch $(path1) chrom.sizes /san/sanvol1/scratch/est/$DB/est.2bit /scratch/hg/$DB/nib $(dir1)/../output/$(file1).filter.psl -score=$(dir1)/../log/$(file1).mrnaMatch.tab -bedOut=$(dir1)/../log/$(file1).mrnaMis.bed -minDiff=5
#ENDLOOP
#ssh pk -T "cd /cluster/data/$DB/bed/pseudo/run.est ; para make spec"
#cd ..
#cat /san/sanvol1/scratch/est/$DB/output/* | sort -k14,14 -k16,16n >~/$DB/pseudo/estFiltered.psl
#awk -f splitTarget.awk < estFiltered.psl 
#for i in `ls chr*.psl` ; do echo $i ; mv $i ${i%%.psl}_estFiltered.psl ; done
#for i in `ls chr*.psl` ; do hgLoadPsl $DB $i ; done
echo "overlapSelect estFiltered.psl.gz retroMrnaInfoLessZnf.bed pseudoEstAll.bed"
overlapSelect estFiltered.psl.gz retroMrnaInfoLessZnf.bed pseudoEstAll.bed
wc -l retroMrnaInfoLessZnf.bed pseudoEstAll.bed
echo "grep -F -f bothZnf.lst retroMrnaInfo$VERSION.bed to retroZnf.bed"
grep -F -f bothZnf.lst retroMrnaInfo$VERSION.bed > retroZnf.bed
#hgsql $DB -N -B < bothZnf.sql > zincFinger.psl

cp all_mrna.psl.gz all_mrnaFiltered.psl.gz
echo "mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid "
mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid > all_mrna.log 2> all_mrna.err
mrnaToGene all_mrnaFiltered.psl.gz -cdsMergeSize=10 -utrMergeSize=10 all_mrna_utr.gp -cdsDb=$DB -keepInvalid > all_mrna_utr.log 2> all_mrna_utr.err

awk '$8>1{print }' all_mrna.gp > all_mrna_multiExon.gp
awk '$8>1{print }' all_mrna_utr.gp > all_mrna_multiExonUTR.gp

echo pwd
pwd
$SCRIPT/analyseExpress.sh $DB $TABLE $VERSION
#nice -4 zcat pseudoMrnaLink*.txt.gz | awk '$28 == 1 && $23 < 50 {print $0}' | 
#awk 'NF==47{$37=$37" -1 0";print $0}NF==48{$37=$37" 0"; print $0}' pseudoGeneLinkSort.bed > pseudoGeneLinkSort.bed
#awk 'NF==47{$37="0 0";print $0}NF==48{print $0}' pseudoGeneLinkNoOverlap.bed > pseudoGeneLinkNoOverlap.bed
# awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html

mkdir -p stats
cd stats
$SCRIPT/doStat $DB $TABLE
cd ..
cd frames

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
echo sort ../../raw.len > mrna.sort.len
wc -l ../../raw.len
sort ../../raw.len > mrna.sort.len
join mrna.sort.len qlist.ext|awk '{OFS="\t";print $3,$2}' > mrnaExt.len 

for i in `ls *.axt` ; do echo $i ; axtToMaf -tPrefix=${DB}. -qPrefix=mrna. $i ../../S1.len mrnaExt.len ../mafRaw/${i%%.axt}.maf ; done

cd ..
mkdir -p maf
mkdir -p mafFrames
cd mafRaw
for i in `ls *.maf` ; do echo $i ; ~/bin/x86_64/mafFilter -componentFilter=../pseudoExpressed.ids $i > ../maf/$i; done
cd ../mafRaw
#echo "genePredToMafFrames $DB $DB pseudoExpressed.single.gp mafFrames/$DB.mafFrames maf/*.maf"
echo "genePredToMafFrames $DB maf/*.maf mafFrames/$DB.mafFrames $DB pseudoExpressed.single.gp "
for i in `ls maf` ; do genePredToMafFrames $DB maf/$i mafFrames/$DB.$i.mafFrames $DB pseudoExpressed.single.gp ; done
cat mafFrames/$DB.*.maf.mafFrames > mafFrames/$DB.mafFrames

#get gene preds for parent genes

hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,1)-1 as cdsStart, substring(c.name,4,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,2,2) = '..' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,4,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version " >mrnaGene.gp
hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,2)-1 as cdsStart, substring(c.name,5,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,3,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,5,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version " >>mrnaGene.gp
hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,3)-1 as cdsStart, substring(c.name,6,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,4,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,6,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version " >>mrnaGene.gp

#get gene preds for refSeq genes
hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,1)-1 as cdsStart, substring(c.name,4,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,2,2) = '..' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,4,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version and name like 'NM%'" >refGene.gp
hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,2)-1 as cdsStart, substring(c.name,5,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,3,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,5,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version and name like 'NM%'" >>refGene.gp
hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,3)-1 as cdsStart, substring(c.name,6,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, retroMrnaInfo p where cds = c.id and substring(c.name,4,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,6,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version and name like 'NM%'" >>refGene.gp
#calculate frames for parent genes 
for i in `ls maf` ; do genePredToMafFrames $DB maf/$i mafFrames/mrna.$i.mafFrames mrna mrnaGene.gp ; echo $i ;done
cat mafFrames/mrna.*.maf.mafFrames > mafFrames/mrna.mafFrames
rm -f mafFrames/mrna.*.maf.mafFrames
echo "hgLoadMafFrames $DB blastzRetroFrames mafFrames/$DB.mafFrames mafFrames/mrna.mafFrames"
#load both sets of frames
hgLoadMafFrames $DB blastzRetroFrames mafFrames/$DB.mafFrames mafFrames/mrna.mafFrames
