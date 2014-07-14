#!/bin/bash 
#
# extract mrna sequences and alignments, then filter to remove mrna with multiple hits that are ambigouous
#
set -beEu -o pipefail
source $1
echo "starting filterMrna.sh $1 for $DB"
tawk '{print $10}' $MRNABASE/gbMrnaOnly.psl |sort |uniq > $OUTDIR/all_mrna.list
cat $MRNABASE/gbMrnaOnly.psl |sort -k10,10> $OUTDIR/all_mrna.qName.psl
pslStats -queryStats $OUTDIR/all_mrna.qName.psl $OUTDIR/all_mrna.stats 
awk 'NF > 1 && $3 > 1{print $1}' $OUTDIR/all_mrna.stats > $OUTDIR/all_mrna.doublehit.list
#awk 'NF > 1 && $3 == 1{print $1}' all_mrna.stats  > all_mrna.singlehit.list
pslSelect -queries=$OUTDIR/all_mrna.doublehit.list $OUTDIR/all_mrna.qName.psl stdout | grep -v "_hap" |grep -v random |sort -k10,10> $OUTDIR/all_mrna_double.sort.psl
cut -f 10 $OUTDIR/all_mrna_double.sort.psl |sort |uniq > $OUTDIR/all_mrna.doublehit.list
$SCRIPT/selectById -not 1 $OUTDIR/all_mrna.doublehit.list 10 $MRNABASE/gbMrnaOnly.psl > $OUTDIR/all_mrna_kept.psl
faToTwoBit -stripVersion $MRNABASE/mrna.fa $OUTDIR/mrnaNoversion.2bit

pslCDnaGenomeMatch $OUTDIR/all_mrna_double.sort.psl $OUTDIR/S1.len $OUTDIR/mrnaNoversion.2bit $NIB $OUTDIR/all_mrna_double.filter.psl -score=$OUTDIR/all_mrna_double.score -bedOut=$OUTDIR/all_mrna_double.bed  -verbose=3 -minDiff=4 -notAlignPenalty=3 > $OUTDIR/all.log 
cat $OUTDIR/all_mrna_kept.psl $OUTDIR/all_mrna_double.filter.psl > $OUTDIR/all_mrnaFiltered.psl
gzip $OUTDIR/all_mrnaFiltered.psl

