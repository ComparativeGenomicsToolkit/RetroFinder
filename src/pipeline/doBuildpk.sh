#!/bin/bash
## Summary: Runs pslPseudo which uses the alignments and other annotations to 
# score the candidate retrogenes. Output is pseudo*.axt and pseudoMrnaLink*.txt
# which are output to a temp directory. Pseudogenes and expressed retrogenes 
# are selected from pseudoMrnaLink*.txt and copied to pseudoMrnaLink*.bed. 
# These output files are moved to the result directory or result/axt directory 
# for the axt files.

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
#echo "${BINDIR}/pslPseudo $RETRO_OPTIONS -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt $TMP/pseudo$1.axt $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log"
echo "/hive/users/hartera/GencodeWG/retroFinder/trunk/bin2/pslPseudo $RETRO_OPTIONS -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt $TMP/pseudo$1.axt $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log"
#${BINDIR}/pslPseudo $RETRO_OPTIONS -verbose=5 -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt /dev/null $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log
# Run pslPseudo to score alignments using annotations and compile information
# for ucscRetroInfoXX table.
/hive/users/hartera/GencodeWG/retroFinder/trunk/bin2/pslPseudo $RETRO_OPTIONS -verbose=5 -cdsFile=$CACHE/cds.tab.gz $DB $CACHE/split/tmp$1.psl $CACHE/chrom.sizes $CACHE/rmsk.bed.gz $CACHE/$NET1.txt.gz $CACHE/$NET2.txt.gz $CACHE/simpleRepeat.bed.gz $CACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl $TMP/pseudoMrnaLink$1.txt /dev/null $CACHE/S1.lst $CACHE/mrna.2bit $CACHE/$GENE2.tab.gz $CACHE/$GENE3.tab.gz $CACHE/$GENE1.tab.gz $CACHE/$NET3.txt.gz $OUT/ortho$1.txt > $TMP/pseudo$1.log
if [ $? == 0 ]; then
ulimit -a >> $TMP/pseudo$1.log
# For each line of pseudoMrnaLink*.txt, if field 32 is 1 or -2, then print 
# line. field 32 is label: 1=pseudogene, -1 not pseudogene -2 expressed 
# retrogene therefore selecting pseudogenes and expressed retrogenes.  
 /bin/awk '$32 == 1 || $32==-2 {print $0}' $TMP/pseudoMrnaLink$1.txt > $TMP/pseudoGeneLink$1.bed
 # gzip log file and move from temp directory to log directory
 /bin/gzip $TMP/pseudo$1.log
 /bin/mv $TMP/pseudo$1.log.gz $LOG/pseudo$1.log.gz
# gzip the pseudoMrnaLink*txt file in the temp directory
 /bin/rm -f $TMP/pseudoMrnaLink$1.txt.gz
 /bin/gzip -f $TMP/pseudoMrnaLink$1.txt
# /bin/rm -f $TMP/pseudoMrnaLink$1.txt
# Move the file (pseudoGeneLink*.bed) of selected pseudogenes and expressed
# retrogenes from the temp directory to the result directory.
 /bin/mv $TMP/pseudoGeneLink$1.bed $RESULT
# /bin/rm -f $TMP/pseudo$1.axt.gz
# /bin/gzip $TMP/pseudo$1.axt
# Copy the gzipped pseudoMraLink*.txt file from the temp directory to the 
# result directory
 /bin/cp -fp $TMP/pseudoMrnaLink$1.txt.gz $RESULT
# gzip the pseudo*.axt file in the temp directory.
 /bin/gzip $TMP/pseudo$1.axt
# Move the pseudo*.axt file to the result/axt directory.
 /bin/mv $TMP/pseudo$1.axt.gz $RESULT/axt
# /bin/rm -rf $TMP


# /bin/awk '{OFS="	";print $14,$16,$17,$10,$1*3-$2,$9}' pseudoMrna.psl > pseudoMrna.txt
#    hgsql $DB -B < export.sql > export.txt
#    /bin/awk -f pseudoToHtml.awk < export.txt  > pseudoMrna.html
else
 /bin/gzip $TMP/pseudo$1.log
 /bin/mv $TMP/pseudo$1.log.gz $LOG/pseudo$1.log.gz
exit 3
fi
