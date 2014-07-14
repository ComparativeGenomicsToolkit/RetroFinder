#!/bin/bash 
#
# extract est sequences and alignments, then filter to remove est with multiple hits that are ambigouous
#
set -bevEu -o pipefail
source $1
#cd $OUTDIR
hgsql $DB -N -B -e "select acc from gbWarn" > $OUTDIR/gbWarn.id
if [[ $SPLIT_SPLICED_EST == 1 ]] ; then
    for chr in `cut -f 1 $OUTDIR/S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$SPLICED_EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/splicedEst.psl ; done
else
    hgsql $DB -N -B -e "select * from $SPLICED_EST" |cut -f2-22 |grep -v -F -f gbWarn.id> $OUTDIR/splicedEst.psl 
fi
rm -f $OUTDIR/splicedEst.psl.gz
gzip $OUTDIR/splicedEst.psl
echo "Starting estFilter.sh $1 for $DB"
rm -f $OUTDIR/est.psl
if [[ $SPLIT_EST == 1 ]] ; then
    for chr in `cut -f 1 S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/est.psl ; done
else
    hgsql $DB -N -B -e "select * from $EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id> $OUTDIR/est.psl 
fi
hgsql $DB -N -B -e "select * from all_est" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/est.psl 
tawk '{print $10}' $OUTDIR/est.psl|sort |uniq > $OUTDIR/est.list &
cat $OUTDIR/est.psl |sort -k10,10> $OUTDIR/est.qName.psl
pslStats -queryStats $OUTDIR/est.qName.psl $OUTDIR/est.stats 
awk 'NF > 1 && $3 > 1{print $1}' $OUTDIR/est.stats > $OUTDIR/est.doublehit.list
awk 'NF > 1 && $3 == 1{print $1}' $OUTDIR/est.stats > $OUTDIR/est.singlehit.list
pslSelect -queries=$OUTDIR/est.singlehit.list $OUTDIR/est.qName.psl $OUTDIR/est.qName.single.psl &
/cluster/data/genbank/bin/x86_64/gbGetSeqs -gbRoot=/cluster/data/genbank -db=$GBDB -native -inclVersion genbank est $OUTDIR/est.fa  -accFile=$OUTDIR/est.doublehit.list
faToTwoBit -stripVersion $OUTDIR/est.fa $OUTDIR/est.2bit
rm -f $OUTDIR/est.fa
pslSelect -queries=$OUTDIR/est.doublehit.list $OUTDIR/est.qName.psl $OUTDIR/est.qName.filter.psl
rm -rf $OUTDIR/est
mkdir -p $OUTDIR/est
${BINDIR}/pslSplit nohead $OUTDIR/est $OUTDIR/est.qName.filter.psl -chunkSize=12000

mkdir -p $OUTDIR/estOutput
mkdir -p $OUTDIR/estLog
mkdir -p run.est

# cd run.est
ls $OUTDIR/est/*psl > $OUTDIR/run.est/list
echo "#LOOP" > $OUTDIR/run.est/template
echo "pslCDnaGenomeMatch \$(path1) S1.len $OUTDIR/est.2bit $NIB $OUTDIR/estOutput/\$(file1).filter.psl -score=$OUTDIR/estLog/\$(file1).mrnaMatch.tab -bedOut=$OUTDIR/estLog/\$(file1).mrnaMis.bed -minDiff=4 -notAlignPenalty=3" >> $OUTDIR/run.est/template
echo "#ENDLOOP" >> $OUTDIR/run.est/template
gensub2 $OUTDIR/run.est/list single $OUTDIR/run.est/template $OUTDIR/run.est/spec
ssh $CLUSTER -T "cd $OUTDIR/run.est ; /parasol/bin/para make spec"

# Now run these commands in $OUTDIR
find $OUTDIR/est.qName.single.psl $OUTDIR/estOutput -name '*.psl' |xargs cat | sort -k14,14 -k16,16n > $OUTDIR/estFiltered.psl
gzip $OUTDIR/estFiltered.psl
#mkdir -p $OUTDIR/estSplit
#pslSplitOnTarget est.psl estSplit
