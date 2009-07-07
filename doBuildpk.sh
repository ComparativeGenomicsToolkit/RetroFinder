#!/bin/bash
set -beEu -o pipefail
SEQ=$1
DB=hg18
NIB=/scratch/data/$DB/nib/
SAN=/hive/users/baertsch
PAN=$SAN/retro/$DB/result
CACHE=$SAN/retro/$DB
ICACHE=$SAN/retro/$DB
LOG=$SAN/retro/$DB/log
OUT=$SAN/retro/$DB/out
cd /tmp
rm -f $LOG/pseudo$1.log
rm -f $LOG/err$1.log
ulimit -d 2200000
ulimit -v 2200000
echo "/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$ICACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $ICACHE/chrom.sizes $ICACHE/rmsk.bed.gz $ICACHE/mouseNet.txt.gz $ICACHE/dogNet.txt.gz $ICACHE/simpleRepeat.bed.gz $ICACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $ICACHE/S1.lst $CACHE/mrna.2bit $ICACHE/refGene.tab.gz $ICACHE/mgcGene.tab.gz $ICACHE/sortedKnownGene.tab.gz $ICACHE/rhesusNet.txt.gz"
/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$ICACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $ICACHE/chrom.sizes $ICACHE/rmsk.bed.gz $ICACHE/mouseNet.txt.gz $ICACHE/dogNet.txt.gz $ICACHE/simpleRepeat.bed.gz $ICACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $ICACHE/S1.lst $CACHE/mrna.2bit $ICACHE/refGene.tab.gz $ICACHE/mgcGene.tab.gz $ICACHE/sortedKnownGene.tab.gz $ICACHE/rhesusNet.txt.gz > /tmp/pseudo$1.log 2> /tmp/err$1.log
if [ $? == 0 ]; then
ulimit -a >> /tmp/err$1.log
 awk '$34 == 1 || $34==-2 {print $0}' /tmp/pseudoMrnaLink$1.txt >/tmp/pseudoGeneLink$1.bed
 mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 mv /tmp/err$1.log $LOG/err$1.log
 rm -f /tmp/pseudoMrnaLink$1.txt.gz
 gzip -f /tmp/pseudoMrnaLink$1.txt
# rm -f /tmp/pseudoMrnaLink$1.txt
 mv /tmp/pseudoGeneLink$1.bed $PAN
 rm -f /tmp/pseudo$1.axt.gz
 gzip /tmp/pseudo$1.axt
#cp -fp /tmp/pseudoMrnaLink$1.txt.gz $PAN
 mv /tmp/pseudo$1.axt.gz $PAN/axt


# awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
else
 mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 mv /tmp/err$1.log $LOG/err$1.log
exit 3
fi
