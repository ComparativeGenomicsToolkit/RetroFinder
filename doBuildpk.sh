#!/bin/bash
set -beEu -o pipefail
SEQ=$1
source $2
CACHE=$OUTDIR
#SAN=/hive/users/baertsch
#PAN=$SAN/retro/$DB/result
#LOG=$SAN/retro/$DB/log
#OUT=$SAN/retro/$DB/out
cd /tmp
rm -f $LOG/pseudo$1.log
rm -f $LOG/err$1.log
ulimit -d 2200000
ulimit -v 2200000
echo "/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz"
/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz > /tmp/pseudo$1.log 2> /tmp/err$1.log
if [ $? == 0 ]; then
ulimit -a >> /tmp/err$1.log
 /bin/awk '$34 == 1 || $34==-2 {print $0}' /tmp/pseudoMrnaLink$1.txt >/tmp/pseudoGeneLink$1.bed
 /bin/mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 /bin/mv /tmp/err$1.log $LOG/err$1.log
 /bin/rm -f /tmp/pseudoMrnaLink$1.txt.gz
 /bin/gzip -f /tmp/pseudoMrnaLink$1.txt
# /bin/rm -f /tmp/pseudoMrnaLink$1.txt
 /bin/mv /tmp/pseudoGeneLink$1.bed $RESULT
 /bin/rm -f /tmp/pseudo$1.axt.gz
 /bin/gzip /tmp/pseudo$1.axt
#/bin/cp -fp /tmp/pseudoMrnaLink$1.txt.gz $RESULT
 /bin/mv /tmp/pseudo$1.axt.gz $RESULT/axt


# /bin/awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    /bin/awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
else
 /bin/mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 /bin/mv /tmp/err$1.log $LOG/err$1.log
exit 3
fi
