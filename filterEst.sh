#!/bin/bash 
#
# extract est sequences and alignments, then filter to remove est with multiple hits that are ambigouous
#
set -beEu -o pipefail
source $1
echo "starting estFilter.sh $1 for $DB"
#for i in `chromsNoY` chrY chrM; do echo $i ; hgsql $DB -N -B -e "select * from ${i}_est" |cut -f2-22 >> est.psl ; done
rm -f est.psl
hgsql $DB -N -B -e "select acc from gbWarn" > gbWarn.id
hgsql $DB -N -B -e "select * from all_est" |cut -f2-22 |grep -v -F -f gbWarn.id>> est.psl 
tawk '{print $10}' est.psl|sort |uniq > est.list &
cat est.psl |sort -k10,10> est.qName.psl
pslStats -queryStats est.qName.psl est.stats 
awk 'NF > 1 && $3 > 1{print $1}' est.stats  > est.doublehit.list
awk 'NF > 1 && $3 == 1{print $1}' est.stats  > est.singlehit.list
pslSelect -queries=est.singlehit.list est.qName.psl est.qName.single.psl &
/cluster/data/genbank/bin/x86_64/gbGetSeqs -gbRoot=/cluster/data/genbank -db=$GBDB -native -inclVersion genbank est est.fa  -accFile=est.doublehit.list
faToTwoBit -stripVersion est.fa est.2bit
rm -f est.fa
pslSelect -queries=est.doublehit.list est.qName.psl est.qName.filter.psl
rm -rf est
mkdir -p $OUTDIR/est
pslSplit nohead $OUTDIR/est est.qName.filter.psl -chunkSize=12000


mkdir -p $OUTDIR/estOutput
mkdir -p $OUTDIR/estLog
mkdir -p run.est

cd run.est
ls ../est/*psl > list
echo "#LOOP" > template
echo "pslCDnaGenomeMatch \$(path1) S1.len $OUTDIR/est.2bit $NIB $OUTDIR/estOutput/\$(file1).filter.psl -score=$OUTDIR/estLog/\$(file1).mrnaMatch.tab -bedOut=$OUTDIR/estLog/\$(file1).mrnaMis.bed -minDiff=5" >> template
echo "#ENDLOOP" >> template
gensub2 list single template spec
ssh $CLUSTER -T "cd $OUTDIR/run.est ; para make spec"
cd ..
cat est.qName.single.psl $OUTDIR/estOutput/* | sort -k14,14 -k16,16n > $OUTDIR/estFiltered.psl
gzip estFiltered.psl
#mkdir -p $OUTDIR/estSplit
#pslSplitOnTarget est.psl estSplit
