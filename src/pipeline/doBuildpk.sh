#!/bin/bash
set -bevEu -o pipefail
SEQ=$1
source $2
CACHE=$OUTDIR
# TMP=/scratch/tmp/retro.$$.$USER
# $TMPDIR is an environment variable that points to /data/tmp
TMP=$TMPDIR/retro.$$.$USER
mkdir -p $TMP
cd $TMP
rm -f $LOG/pseudo$1.log
#ulimit -d 2800000
#ulimit -v 2800000
echo "${BINDIR}/pslPseudo $RETRO_OPTIONS -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt $TMP/pseudo$1.axt $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log"
${BINDIR}/pslPseudo $RETRO_OPTIONS -verbose=5 -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt /dev/null $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log
if [ $? == 0 ]; then
ulimit -a >> $TMP/pseudo$1.log
 /bin/awk '$32 == 1 || $32==-2 {print $0}' $TMP/pseudoMrnaLink$1.txt > $TMP/pseudoGeneLink$1.bed
 /bin/gzip $TMP/pseudo$1.log
 /bin/mv $TMP/pseudo$1.log.gz $LOG/pseudo$1.log.gz
 /bin/rm -f $TMP/pseudoMrnaLink$1.txt.gz
 /bin/gzip -f $TMP/pseudoMrnaLink$1.txt
# /bin/rm -f $TMP/pseudoMrnaLink$1.txt
 /bin/mv $TMP/pseudoGeneLink$1.bed $RESULT
# /bin/rm -f $TMP/pseudo$1.axt.gz
# /bin/gzip $TMP/pseudo$1.axt
 /bin/cp -fp $TMP/pseudoMrnaLink$1.txt.gz $RESULT
# /bin/mv $TMP/pseudo$1.axt.gz $RESULT/axt
# /bin/rm -rf $TMP


# /bin/awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    /bin/awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
else
 /bin/gzip $TMP/pseudo$1.log
 /bin/mv $TMP/pseudo$1.log.gz $LOG/pseudo$1.log.gz
exit 3
fi
