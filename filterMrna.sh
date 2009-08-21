#!/bin/bash 
#
# extract mrna sequences and alignments, then filter to remove mrna with multiple hits that are ambigouous
#
set -beEu -o pipefail
source $1
echo "starting filterMrna.sh $1 for $DB"
hgsql $DB -N -B -e "select * from all_mrna" |cut -f2-22 > all_mrna.psl 
tawk '{print $10}' all_mrna.psl|sort |uniq > all_mrna.list
cat all_mrna.psl |sort -k10,10> all_mrna.qName.psl
pslStats -queryStats all_mrna.qName.psl all_mrna.stats 
awk 'NF > 1 && $3 > 1{print $1}' all_mrna.stats  > all_mrna.doublehit.list
#awk 'NF > 1 && $3 == 1{print $1}' all_mrna.stats  > all_mrna.singlehit.list
pslSelect -queries=all_mrna.doublehit.list all_mrna.qName.psl stdout | grep -v "_hap" |grep -v random |sort -k10,10> all_mrna_double.sort.psl
cut -f 10 all_mrna.qName.double.psl > all_mrna.doublehit.list
grep -F -v -f all_mrna.doublehit.list all_mrna.psl > all_mrna_kept.psl
/cluster/data/genbank/bin/x86_64/gbGetSeqs -gbRoot=/cluster/data/genbank -db=$GBDB -native -inclVersion genbank all_mrna all_mrna.fa  -accFile=all_mrna.doublehit.list
faToTwoBit -stripVersion all_mrna.fa all_mrna.2bit
rm -f all_mrna.fa

pslCDnaGenomeMatch all_mrna_double.sort.psl S1.len all_mrna.2bit /cluster/data/hg18/nib all_mrna_double.filter.psl -score=all_mrna_double.score -bedOut=all_mrna_double.bed  -verbose=3 -minDiff=5 > all.log 
cat all_mrna_kept.psl all_mrna_double.filter.psl > all_mrnaFiltered.psl
gzip all_mrnaFiltered.psl
