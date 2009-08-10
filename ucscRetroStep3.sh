#!/bin/bash 
set -beEu -o pipefail
source $1
echo "starting ucscRetroStep3.sh $1 for $DB"
mkdir -p $OUTDIR
cd $OUTDIR
mkdir -p result
mkdir -p result/axt
mkdir -p log
mkdir -p out

ls $NIB/ > $OUTDIR/S1.lst
cp $GENOME/$DB/chrom.sizes .

hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from all_mrna" > all_mrna.psl 
hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from refSeqAli" >> all_mrna.psl
rm -f all_mrna.psl.gz
gzip all_mrna.psl
#cat $RMSK/*.out |awk '{OFS="\t";print $5,$6,$7}' | grep -v position|grep -v sequence | tawk 'length($0)>2{print $0}' > rmsk.bed
rm -f rmsk.bed
if [ $RMSK == "rmsk" ]; then hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from rmsk" >> rmsk.bed ;
else 
for i in `cut -f 1 chrom.sizes` ;do hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from ${i}_rmsk" >> rmsk.bed ; done  ; fi
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET1 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET1.txt
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET2 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET2.txt 
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET3 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET3.txt 
hgsql $DB -N -B -e "select chrom, chromStart, chromEnd from simpleRepeat" > simpleRepeat.bed 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 order by chrom, txStart, txEnd" | sort -k2,2 -k4,4n > $GENE1.tab
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2" > $GENE2.tab 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3" > $GENE3.tab 
hgsql $DB -N -B -e "select acc, version, name, type from refSeqAli a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" > cds.tab
hgsql $DB -N -B -e "select acc, version, name, type from all_mrna a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" >> cds.tab
rm -f *.tab.gz *.txt.gz *.bed.gz
gzip *.tab &
gzip *.txt &
gzip *.bed 
#cp -p *.gz /san/sanvol1/scratch/retro/$DB
#cp -p S1.lst /san/sanvol1/scratch/retro/$DB

#filter mrna and est with pslCDnaGenomeMatch
mkdir -p $TMPEST/mrna

rm -f mrna.2bit
ln $MRNABASE/mrna.2bit . -s
#gbGetSeqs genbank mrna mrna.fa -native -db=$DB -get=seq -gbRoot=/cluster/data/genbank &
#faToTwoBit mrna.fa /san/sanvol1/scratch/est/mrna/mrna.2bit
#/cluster/data/genbank/bin/x86_64/gbGetSeqs -db=hg -inclVersion -native -gbRoot=/cluster/data/genbank \
#   genbank mrna stdout | tr acgt ACGT > mrna.fa 
#/cluster/data/genbank/bin/x86_64/gbGetSeqs -db=hg -inclVersion -native -gbRoot=/cluster/data/genbank \
#   refSeq mrna stdout | tr acgt ACGT > refseq.fa 
#gbGetSeqs refSeq  mRNA refSeq.fa -native -db=$DB -get=seq -gbRoot=/cluster/data/genbank &
#gbGetSeqs refSeq  mRNA refSeqVers.fa -inclVersion -native -db=$DB -get=seq -gbRoot=/cluster/data/genbank

#hgsql $DB -N -B -e "select concat(acc,'.',version) as acc, r.qStart, r.qEnd from gbCdnaInfo i, cds c, refSeqAli r  where i.cds = c.id  and acc like 'NM%' and acc = r.qName" > refSeqCds.bed
#faToTwoBit refSeq.fa refSeq.2bit
#hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from refSeqAli" > refSeqAli.psl
#nohup nice pslCDnaGenomeMatch refSeqAli.new.psl /cluster/data/$DB/chrom.sizes refSeq.2bit  /cluster/data/$DB/nib/ refSeqAli.filter.psl -score=refSeq.tab -bedOut=refSeqMis.bed -passthru &
#gzip refSeqAli.filter.psl
#
#hgsql $DB -N -B -e "select * from all_mrna" | cut -f2-22 > all_mrnaNew.psl
#nohup nice pslCDnaGenomeMatch all_mrnaNew.psl /cluster/data/$DB/chrom.sizes mrna.2bit  /cluster/data/$DB/nib/ all_mrnaFiltered.psl -score=mrna.tab -bedOut=rmrnaMis.bed -passthru &
#gzip all_mrnaFiltered.psl
#zcat all_mrnaFiltered.psl.gz refSeqAli.filter.psl.gz > all_mrnaRefSeq.psl
#sort -k10 all_mrnaFiltered.psl.gz > all_mrnaNew.qName.psl
#gzip all_mrnaNew.qName.psl
#pslSplit nohead  $TMPEST/mrna/input/ all_mrnaNew.qName.psl.gz -chunkSize=2000
#cd run.mrna
#ls $TMPEST/mrna/input/*.psl > list
#gensub2 list single gsub spec
#para make
#cat $TMPEST/mrna/output/*.psl > all_mrnaFiltered.psl

####
# run retro pipeine on cluster
####
mkdir -p $RETRODIR
cd $OUTDIR
mkdir -p run.0
cd run.0
wc -l ../split/*.psl |grep -v total | grep -v acc.lst| sort -nr | awk '{print $2}' |sed -e 's,../split/tmp,,'|sed -e 's/.psl//' > list
cp $OUTDIR/$1 .
echo "#LOOP" > gsub
echo "$SCRIPT/doBuildpk.sh \$(path1) $OUTDIR/$1 {check out exists $RESULT/pseudoGeneLink\$(path1).bed} " >> gsub
echo "#ENDLOOP" >> gsub

gensub2 list single gsub jobList
echo "Job Count"
wc -l spec
ssh -T $CLUSTER "cd $OUTDIR/run.0 ; para make jobList"
#1190 jobs in batch
#2098 jobs (including everybody's) in Parasol queue.
#Checking finished jobs
#..................
#Completed: 1190 of 1190 jobs
#CPU time in finished jobs:      37794s     629.91m    10.50h    0.44d  0.001 y
#IO & Wait Time:                  3482s      58.03m     0.97h    0.04d  0.000 y
#Average job time:                  35s       0.58m     0.01h    0.00d
#Longest running job:                0s       0.00m     0.00h    0.00d
#Longest finished job:             239s       3.98m     0.07h    0.00d
#Submission to last job:          1345s      22.42m     0.37h    0.02d

#post process
#cd /cluster/data/$DB/bed/pseudo
echo "check parasol status and then run ucscRetroStep4.sh DEF"
cd ..

######################
# get ests and filter
######################
hgsql $DB -N -B -e "show tables like '%est'" > est.lst

rm -f est.psl est.psl.gz
for i in `cat est.lst` ; do hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from ${i}" >> est.psl ; done
gzip est.psl
rm -f splicedEst.psl splicedEst.psl.gz
hgsql $DB -N -B -e "show tables like '%intronEst'" > splicedEst.lst
for i in `cat splicedEst.lst` ; do hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from ${i}" >> splicedEst.psl ; done
gzip splicedEst.psl
ln est.psl.gz estFiltered.psl.gz -s

exit

sort -k10 est.psl > est.qName.psl
gbGetSeqs genbank est est.fa -native -db=$DB -get=seq -gbRoot=/cluster/data/genbank 
faToTwoBit -stripVersion -ignoreDups est.fa est.2bit
cp -p est.fa /san/sanvol1/scratch/est/$DB
pslSplit nohead  /san/sanvol1/scratch/est/$DB/input/ est.qName.psl -chunkSize=2000
cd run.est
rm -f est.psl est.qName.psl
ls /san/sanvol1/scratch/est/$DB/input/*.psl > list
gensub2 list single gsub spec
para make
cat /san/sanvol1/scratch/est/$DB/output/*.psl > all_estFiltered.psl
awk -f splitTarget.awk < estFiltered.psl 
for i in `ls chr*estF*psl` ; do echo $i ; hgLoadPsl $DB $i -noTNameIx -fastLoad; done



nohup nice pslCDnaGenomeMatch estSort.psl S1.len est.2bit /cluster/data/$DB/nib/ estFinal.psl -score=estOut.tab -bedOut=estMis.bed -passthru &

#filter  retros using expression
overlapSelect estFiltered.psl retroMrnaInfo600.fix.bed pseudoEstAll.bed

mrnaToGene all_mrnaFiltered.psl -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid > all_mrna.log 2> all_mrna.err
mrnaToGene all_mrnaFiltered.psl -cdsMergeSize=10 -utrMergeSize=10 all_mrna_utr.gp -cdsDb=$DB -keepInvalid > all_mrna_utr.log 2> all_mrna_utr.err

awk '$8>1{print }' all_mrna.gp > all_mrna_multiExon.gp
awk '$8>1{print }' all_mrna_utr.gp > all_mrna_multiExonUTR.gp

doConvertSvm.sh pseudoGeneLinkSortFilter.fix.bed 16
cat svm.1/output |awk '{print $1+(rand()*0.1)-(rand()*0.1)} '> svm.test16.fix.dat 
paste svm.test16.fix.dat pseudoScore.txt > svm.dat

awk '{print $1+(rand()*0.1)-(rand()*0.1)} ' svm.1/output | paste - pseudoScore.txt |awk $1<100 && $1 > -100{print $0}'> svm.dat


#overlap with vega
overlapSelect vegaGene.gp pseudoGeneLink425.fix.bed pseudoVegaGene.gp

#exon shuffle
hgsql $DB -N -B -e "select * from retroMrnaInfo where conservedSpliceSites= 0" |cut -f1-56 > retroMrnaInfo.shuffle.bed 
#remove znf
hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'ZNF%' and a1.kgID = a2.kgID" |sort -u > kgZnf.lst
cat kgZnf.lst refZnf.lst > bothZnf.lst
grep -v -F -f bothZnf.lst retroMrnaInfo.shuffle.bed > retroMrnaInfoLessZnf.bed
tawk '{$55=$14;$14=0;print $0}' retroMrnaInfoLessZnf.bed > retroMrnaInfo.fix.bed 
tawk '$5> 650{print $0}' retroMrnaInfo.fix.bed > retroMrnaInfo650.fix.bed
tawk '$5> 650{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' retroMrnaInfo.bed > retroMrnaInfo.12.bed
hgsql $DB -N -B < refGeneMultiExon.sql |cut -f2-12 > refGeneMultiExon.genePred
ldHgGene $DB rbRefGeneMulti refGeneMultiExon.genePred -predTab 

featureBits $DB refGene:cds -bed=refGeneCds.bed
featureBits $DB rbRefGeneMulti;
52619505 bases of 2881515245 (1.826%) in intersection
featureBits $DB refGene;     
54738314 bases of 2881515245 (1.900%) in intersection
featureBits $DB rbRefGeneMulti refGeneCds.bed -bed=refGeneMultiCds.bed 
29845142 bases of 2881515245 (1.036%) in intersection
overlapSelect refGeneMultiCds.bed retroMrnaInfo650.fix.bed pseudoRefGeneCds.bed
overlapSelect refGeneMultiCds.bed retroMrnaInfo.12.bed -statsOutput pseudoRefGeneCds.out
overlapSelect refGeneMultiCds.bed retroMrnaInfo.12.bed pseudoRefGeneCds50.bed -overlapThreshold=0.50

~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB pseudoRefGeneCds.bed ~/.html/shuffle/all 
overlapSelect pnasretroExp.from.hg17.bed pseudoRefGeneCds.bed notMatchKass.bed -nonOverlapping
~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB notMatchKass.bed ~/.html/shuffle/notK
wc -l notMatchKass.bed


awk '{OFS=";";$1=$1;print $4"."$1"."$2,$0}' retroMrnaInfo650.fix.bed|sort > retroMrnaInfo650.txt
awk '{OFS=";";$1=$1;print $4"."$1"."$2,$0}' pseudoShuffle.bed |sort> pseudoShuffle.txt 
join -t ";" pseudoShuffle.txt retroMrnaInfo650.txt|awk -F";" '{OFS="\t";$1=$1;print $0}' |cut -f14-70 > pseudoShuffle.long.bed
tawk '{print $4,$1,$2,$3,$5,$6,$41,$16,$17,$18,$19,$55,$32,$54,$31,$15,$33,$35,$24,$20,$10,$26,$39,$22,$23,$56,$40,$50}' head.fix > pseudoShuffle.long.txt
tawk '{print $4,$1,$2,$3,$5,$6,$41,$16,$17,$18,$19,$55,$32,$54,$31,$15,$33,$35,$24,$20,$10,$26,$39,$22,$23,$56,$40,$50}' pseudoShuffle.long.bed >> pseudoShuffle.long.txt
#
#ES cells
cut -f1-10 knownGene.tab |genePredSingleCover stdin knownGene.sc.gp
awk '{print $7"\t"$2}' ESExpressedGenes.final |sed -e 's/:/    /'|sed -e 's/-/ /' > ESExpressedGenes.bed
awk '{print $7"\t"$2}' ESExpressedlog2.final |sed -e 's/:/     /'|sed -e 's/-/ /' |bedSort stdin ESExpressedlog2.bed
tawk '{print $16,$17,$18,$4,$5,$6,$7,$8,$9,$10,$11,$12}' pseudoGeneLink600.bed > pseudoParent.bed
bedSort pseudoParent.bed pseudoParent.sort.bed
bedOverlap pseudoParent.sort.bed stdout | cut -f2-13 > pseudoParent.sc.bed
  3819 pseudoParent.sc.bed
overlapSelect ESExpressedGenes.bed pseudoParent.sc.bed ESEparent.bed
overlapSelect knownGene.sc.gp ESExpressedlog2.bed ESExpressedlog2.sc.bed
  2124 ESExpressedlog2.sc.bed
overlapSelect ESExpressedlog2.sc.bed pseudoParent.sc.bed ESEparentlog.bed
  2124 ESExpressedlog2.sc.bed
  3819 pseudoParent.sc.bed
   953 ESEparentlog.bed

##pnas compare since rhesus
 hgsql $DB -N -B -e "select chrom,chromStart, chromEnd, name, score from pnasretro" > pnasretro.bed
 overlapSelect pnasretro.bed pseudoExpressedRhesus.bed pseudoExpRhesusNotPnas.bed -nonOverlapping
 overlapSelect pnasretro.bed pseudoExpressedRhesus.bed pnasExpRhesusMatch.bed
 wc -l pnasExpRhesusMatch.bed pseudoExpRhesusNotPnas.bed pseudoExpressedRhesus.bed
    113 pnasExpRhesusMatch.bed
     75 pseudoExpRhesusNotPnas.bed
    188 pseudoExpressedRhesus.bed
    376 total
 cat pseudoExpRhesusNotPnas.bed pnasExpRhesusMatch.bed > pseudoExpressedRhesus.sort.bed
 ~markd/bin/bedToHtmlDir -browser-url http://hgwdev-baertsch.cse.ucsc.edu -context-bases 300 $DB pseudoExpressedRhesus.sort.bed ~/.html/retroRhesus/

##krish indel
# wget http://www.soe.ucsc.edu/~krish/indel/hg18_cat.bed
# mv hg18_cat.bed indelible.raw.bed
# awk '$3-$2>100{print $0}' indelible.raw.bed > indelible.bed
# grep hg18_insert indelible.bed > indelibleHuman.bed
# grep -v panTro2_delete indelible.bed > indelibleHumanIndel.bed
# overlapSelect indelibleHuman.bed pseudoGeneLink600.fix.bed pseudoHuman.bed


hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,1)-1 as cdsStart, substring(c.name,4,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, pseudoGeneLink3 p where cds = c.id and substring(c.name,2,2) = '..' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,4,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version" >mrnaGene.gp

hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,2)-1 as cdsStart, substring(c.name,5,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, pseudoGeneLink3 p where cds = c.id and substring(c.name,3,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,5,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version" >>mrnaGene.gp

hgsql $DB -N -B -e "select p.name, p.name as chrom, '+', 0 as txStart ,gEnd-gStart as txEnd,substring(c.name,1,3)-1 as cdsStart, substring(c.name,6,10) as cdsEnd,1 as exonCount, 0 as exonStarts, gEnd-gStart as exonEnds from cds c, gbCdnaInfo g, pseudoGeneLink3 p where cds = c.id and substring(c.name,4,2) = '..' and substring(c.name,1,1) != '<' and acc = substring(p.name,1,instr(p.name,'.')-1) and substring(c.name,6,1) != '>' and substring(p.name,instr(p.name,'.')+1,1) = version" >>mrnaGene.gp

hgsql $DB -N -B  -e "select acc, name from all_mrna a , gbCdnaInfo g , cds c where qName = acc and cds = c.id " > cds.txt
mrnaToGene all_mrna.psl -cdsFile=cds.txt  mrnaPred.gp -quiet

## blastzRetro
pushd /san/sanvol1/scratch/pseudo/results/axt
zcat *.axt.gz > all.axt
popd
cat /san/sanvol1/scratch/pseudo/$DB/result/axt/all.axt | axtSplitByTarget stdin frames/maf

hgsql $DB -N -B -e "select * from kgTxInfo t, knownGene k where k.name = t.name and category = 'coding'" > kgCoding.gp
hgsql $DB -N -B -e "select * from kgTxInfo t, knownGene k where k.name = t.name and category <> 'coding'" > kgNotCoding.gp

