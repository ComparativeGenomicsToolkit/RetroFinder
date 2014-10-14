#!/bin/bash 
#
# Summary: Create a 2bit file of mRNA sequences. For Blat alignments of mRNAs, 
# filter to remove mRNA with multiple hits that are ambigouous by picking the best
# hits. 
#
set -beEu -o pipefail
source $1
echo "starting filterMrna.sh $1 for $DB"
# Get the qName (field 10) from the PSL of Blat alignments of the native mRNAs,
# and sort and uniq the list  
tawk '{print $10}' $MRNABASE/gbMrnaOnly.psl |sort |uniq > $OUTDIR/all_mrna.list
# Sort the mRNA Blat alignments on the qName (accession)
cat $MRNABASE/gbMrnaOnly.psl |sort -k10,10> $OUTDIR/all_mrna.qName.psl
# Get PSL stats on the sorted PSL of mRNA Blat alignments 
pslStats -queryStats $OUTDIR/all_mrna.qName.psl $OUTDIR/all_mrna.stats 
# if number of fields is > 1 and field 3 is > 1 then print field 1 from 
# the PSL stats file, so if alnCnt is > 1 then print qName (field 1). Prints
# accessions of sequences that have more than one Blat alignment to the genome. 
awk 'NF > 1 && $3 > 1{print $1}' $OUTDIR/all_mrna.stats > $OUTDIR/all_mrna.doublehit.list
#awk 'NF > 1 && $3 == 1{print $1}' all_mrna.stats  > all_mrna.singlehit.list
# Select those mRNAs from the list from the PSL of Blat alignmente of mRNAs
# and remove haplotypes and random chroms and sort on qName
pslSelect -queries=$OUTDIR/all_mrna.doublehit.list $OUTDIR/all_mrna.qName.psl stdout | grep -v "_hap" |grep -v random |sort -k10,10> $OUTDIR/all_mrna_double.sort.psl
# Get field 10 (qName) of this PSL and sort and uniq and output to file
cut -f 10 $OUTDIR/all_mrna_double.sort.psl |sort |uniq > $OUTDIR/all_mrna.doublehit.list
# select non matching lines to field 1 (accessions/qNames) of 
# all_mrna.doublehit.list from $MRNABASE/gbMrnaOnly.psl using field 10 (qName) of 
# the PSL. The output is a PSL of selected Blat alignments of mRNAs with only a 
# single alignent with those to haplotyhpe and random chroms removed.  
$SCRIPT/selectById -not 1 $OUTDIR/all_mrna.doublehit.list 10 $MRNABASE/gbMrnaOnly.psl > $OUTDIR/all_mrna_kept.psl
# Make a 2bit file of the mRNA sequences in mrna.fa with versions removed from ids
faToTwoBit -stripVersion $MRNABASE/mrna.fa $OUTDIR/mrnaNoversion.2bit

# Input PSL is set of single hit mRNAs and the output of filtered alignments for
# best matches and cases with no filtering are in all_mrna_double.filter.psl
# -notAlignPenalty=N score non-aligning bases with 1/N (default 8)
# -minDiff=N minimum difference in score to filter out 2nd best hit (default 5)
# -bedOut=bed output file of mismatches.
pslCDnaGenomeMatch $OUTDIR/all_mrna_double.sort.psl $OUTDIR/S1.len $OUTDIR/mrnaNoversion.2bit $NIB $OUTDIR/all_mrna_double.filter.psl -score=$OUTDIR/all_mrna_double.score -bedOut=$OUTDIR/all_mrna_double.bed  -verbose=3 -minDiff=4 -notAlignPenalty=3 > $OUTDIR/all.log 
# Cat together the single hit Blat alignments and filtered PSL (for multiple hit
# mRNAs) for mRNA Blat alignments and gzip file
cat $OUTDIR/all_mrna_kept.psl $OUTDIR/all_mrna_double.filter.psl > $OUTDIR/all_mrnaFiltered.psl
gzip $OUTDIR/all_mrnaFiltered.psl

