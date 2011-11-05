#!/bin/bash 
#
# extract est sequences and alignments, then filter to remove est with multiple hits that are ambigouous
#
set -beEu -o pipefail
source $1
hgsql $DB -N -B -e "select acc from gbWarn" > gbWarn.id
if [[ $SPLIT_SPLICED_EST == 1 ]] ; then
    for chr in `cut -f 1 S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$SPLICED_EST" |cut -f2-22 |grep -v -F -f gbWarn.id>> splicedEst.psl ; done
else
    hgsql $DB -N -B -e "select * from $SPLICED_EST" |cut -f2-22 |grep -v -F -f gbWarn.id> splicedEst.psl 
fi
rm -f splicedEst.psl.gz
gzip splicedEst.psl
echo "starting estFilter.sh $1 for $DB"
rm -f est.psl
if [[ $SPLIT_EST == 1 ]] ; then
    for chr in `cut -f 1 S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$EST" |cut -f2-22 |grep -v -F -f gbWarn.id>> est.psl ; done
else
    hgsql $DB -N -B -e "select * from $EST" |cut -f2-22 |grep -v -F -f gbWarn.id> est.psl 
fi
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
echo "$BINDIR/pslCDnaGenomeMatch \$(path1) S1.len $OUTDIR/est.2bit $NIB $OUTDIR/estOutput/\$(file1).filter.psl -score=$OUTDIR/estLog/\$(file1).mrnaMatch.tab -bedOut=$OUTDIR/estLog/\$(file1).mrnaMis.bed -minDiff=4 -notAlignPenalty=3" >> template
echo "#ENDLOOP" >> template
gensub2 list single template spec
ssh $CLUSTER -T "cd $OUTDIR/run.est ; /parasol/bin/para make spec"
cd ..
find est.qName.single.psl $OUTDIR/estOutput -name '*.psl' |xargs cat | sort -k14,14 -k16,16n > $OUTDIR/estFiltered.psl
gzip estFiltered.psl
#mkdir -p $OUTDIR/estSplit
#pslSplitOnTarget est.psl estSplit
