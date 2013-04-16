#!/bin/bash 
#
# extract mrna sequences and alignments, then filter to remove mrna with multiple hits that are ambigouous
#
set -beEu -o pipefail
source $1
echo "starting filterMrna.sh $1 for $DB"
tawk '{print $10}' $MRNABASE/gbMrnaOnly.psl |sort |uniq > all_mrna.list
cat $MRNABASE/gbMrnaOnly.psl |sort -k10,10> all_mrna.qName.psl
pslStats -queryStats all_mrna.qName.psl all_mrna.stats 
awk 'NF > 1 && $3 > 1{print $1}' all_mrna.stats  > all_mrna.doublehit.list
#awk 'NF > 1 && $3 == 1{print $1}' all_mrna.stats  > all_mrna.singlehit.list
pslSelect -queries=all_mrna.doublehit.list all_mrna.qName.psl stdout | grep -v "_hap" |grep -v random |sort -k10,10> all_mrna_double.sort.psl
cut -f 10 all_mrna_double.sort.psl |sort |uniq > all_mrna.doublehit.list
$SCRIPT/selectById -not 1 all_mrna.doublehit.list 10 $MRNABASE/gbMrnaOnly.psl > all_mrna_kept.psl
faToTwoBit -stripVersion $MRNABASE/mrna.fa mrnaNoversion.2bit

pslCDnaGenomeMatch all_mrna_double.sort.psl S1.len mrnaNoversion.2bit /hive/data/genomes/$DB/nib all_mrna_double.filter.psl -score=all_mrna_double.score -bedOut=all_mrna_double.bed  -verbose=3 -minDiff=4 -notAlignPenalty=3 > all.log 
cat all_mrna_kept.psl all_mrna_double.filter.psl > all_mrnaFiltered.psl
gzip all_mrnaFiltered.psl

