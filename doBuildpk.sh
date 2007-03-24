#!/bin/bash
SEQ=$1
DB=hg18
NIB=/scratch/hg/gs.19/build36/bothMixedNib/
SAN=/san/sanvol1/scratch
#SAN=/cluster/panasas/home/store/
PAN=$SAN/pseudo/$DB/result
CACHE=$SAN/pseudo/$DB
ICACHE=$SAN/pseudo/$DB
#LOG=/cluster/panasas/home/store/$DB/pseudo/log
LOG=$SAN/pseudo/$DB/log
#LOG=/cluster/bluearc/$DB/pseudo/log
OUT=/cluster/bluearc/$DB/pseudo/out
#OUT=$SAN/pseudo/$DB/out
#OUT=/cluster/panasas/home/store/$DB/pseudo/out
cd /tmp
rm -f $LOG/pseudo$1.log
rm -f $LOG/err$1.log
ulimit -d 2200000
ulimit -v 2200000
echo "/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$ICACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $ICACHE/chrom.sizes $ICACHE/rmsk.bed.gz $ICACHE/mouseNet.txt.gz $ICACHE/dogNet.txt.gz $ICACHE/simpleRepeat.bed.gz $ICACHE/all_mrnaRefSeq.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $ICACHE/S1.lst $CACHE/mrna.2bit $ICACHE/refGene.tab.gz $ICACHE/mgcGene.tab.gz $ICACHE/sortedKnownGene.tab.gz $ICACHE/rheMac2Net.txt.gz"
/cluster/home/baertsch/bin/x86_64/pslPseudo -verbose=5 -minAli=0.98 -nearTop=0.005 -cdsFile=$ICACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $ICACHE/chrom.sizes $ICACHE/rmsk.bed.gz $ICACHE/mouseNet.txt.gz $ICACHE/dogNet.txt.gz $ICACHE/simpleRepeat.bed.gz $ICACHE/all_mrnaRefSeq.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $ICACHE/S1.lst $CACHE/mrna.2bit $ICACHE/refGene.tab.gz $ICACHE/mgcGene.tab.gz $ICACHE/sortedKnownGene.tab.gz $ICACHE/rheMac2Net.txt.gz > /tmp/pseudo$1.log 2> /tmp/err$1.log
if [ $? == 0 ]; then
ulimit -a >> /tmp/err$1.log
 awk '$34 == 1 || $34==-2 {print $0}' /tmp/pseudoMrnaLink$1.txt >/tmp/pseudoGeneLink$1.bed
 mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 mv /tmp/err$1.log $LOG/err$1.log
 gzip -f /tmp/pseudoMrnaLink$1.txt
# rm -f /tmp/pseudoMrnaLink$1.txt
 mv /tmp/pseudoGeneLink$1.bed $PAN
 gzip /tmp/pseudo$1.axt
#cp -fp /tmp/pseudoMrnaLink$1.txt.gz $PAN
 mv /tmp/pseudo$1.axt.gz $PAN/axt


#hgLoadPsl $DB pseudoMrna.psl
#cat << '_EOF_' > loadLink.sql
#delete from pseudoGeneLink where assembly = $DB;
#_EOF_
#hgsql $DB < loadLink.sql
#hgLoadBed hg16 pseudoGeneLink pseudoGeneLink.bed -oldTable
# awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
else
 mv /tmp/pseudo$1.log $LOG/pseudo$1.log
 mv /tmp/err$1.log $LOG/err$1.log
exit 3
fi
