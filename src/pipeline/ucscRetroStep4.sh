#!/bin/bash 
## Summary: Filters out low scoring hits based on the retrogene score and the 
# lastz alignment score (axtScore) and creates one pseudogeneMrnaLink file. 
# Cats pseudo[0-9]*.psl files and ortho files into single files. Splits the 
# pseudogeneMrnaLink file into smaller BED files, one per chromsome, and runs 
# bedOverlap on these files on the cluster to remove overlapping BEDs, the 
# highest scoring one being kept. 

# runs small cluster job to remove overlaps by picking highest scoring retro in each locus
#
set -bevEu -o pipefail
source $1
echo "---------------------------------------------------------------------------------------"
echo "Starting ucscRetroStep4.sh $1 on $HOST- catting output from retro pipeline cluster run"
echo "---------------------------------------------------------------------------------------"
# cd $OUTDIR
# Checks for S1.len file in $OUTDIR. This is the file of chromosome names 
# and sizes (no random chroms or chrM). If the file does not exist, prints
# an error message and aborts. 
if [[ -s $OUTDIR/S1.len ]] ; then
    echo "$OUTDIR/S1.len exists with `wc -l S1.len` rows"
else
    echo "Please create $OUTDIR/S1.len from chrom.sizes without random chroms or chrM."
    exit 3
fi
# Checks the job count from the cluster run from ucscRetroStep3.sh and the 
# number of output files. Code that checks the numbers is commented out.  
# NOTE: doesn't parasol check that all jobs are completed and output files
# exist? Also it is checking the numnber of pseudoMrnaLink*.bed and in the 
# (unlikely) case that a pseudoMrnaLink*txt output file has no pseudogenes
# or expressed retrogenes then the number of these files will not match the 
# original job list count.  
wc -l $OUTDIR/run.0/jobList | awk '{print $1}'> $OUTDIR/jobs.cnt
pushd $RESULT ; ls pseudoGeneLink[0-9]*.bed | wc -l | awk '{print $1}'> $OUTDIR/jobs.lst
popd
#echo Check Job Count from cluster run
#diff jobs.cnt jobs.lst 
#if [ $? != 0 ]; then
#  echo missing jobs aborting
#  exit 3
#fi

# Filter out low scoring hits as defined by retroscore and alignment score 
# to parent gene
rm -f $OUTDIR/pseudoGeneLinkSortFilter.bed.gz
echo catting Sorting and Filtering pseudoGeneLinkSortFilter.bed
# pushd $RESULT ; 
# Cat the pslPseudo output files and if field 5 is greater than 300 
# and field 14 is either > 10,000 or equal to -1 then this line is printed 
# and input to the removeTandemDups script. 
# field 5 is the retrogene score (isn't this already thresholded elsewhere?
# CHECK THIS) and field 14 is the axtScore (lastz score of the parent mRNA 
# aligned to the retrogene locus).
# removeTandemDups removes retros that have parents on the same chromosome are 
# in 200 kb of the retro to remove false positives due to tandemly duplicated 
# gene families.
# The script consists of one line: tawk '$1==$16{g=($17+$18)/2;d=$3-g;dd=sqrt(d*d)}$1!=$16 || ($1==$16 && dd > 200000){print $0}'
# The output is piped to sort and sorted on chrom, (start, end), retrogene id 
# and the output is in pseudoGeneLinkSortFilter.bed. 
cat $RESULT/pseudoGeneLink[0-9]*.bed | tawk '$5 > 300 && ($14 > 10000 || $14 == -1) {OFS="\t";print $0}'| $SCRIPT/removeTandemDups | sort -k1,1 -k2,3n -k4,4 -T /dev/shm/retroSort > $OUTDIR/pseudoGeneLinkSortFilter.bed; #/bin/rm $RESULT/pseudoGeneLink[0-9]*.bed
# popd
# Get count of lines in pseudoGeneLinkSortFilter.bed
wc -l $OUTDIR/pseudoGeneLinkSortFilter.bed
# pushd $OUT
# cat together pseudo[0-9]*.psl (WHERE DID THESE COME FROM?) into pseudo.psl
cat $OUT/pseudo[0-9]*.psl > $OUTDIR/pseudo.psl #;/bin/rm $OUT/pseudo[0-9]*.psl &
# Cat together ortho*.txt files and remove ".txt" from files, output 
# in ortho.txt
cat $OUT/ortho*.txt | sed -e 's/.txt//'> $OUTDIR/ortho.txt
# popd
echo Removing Overlaps
# gzip pseudoGeneLinkSortFilter.bed file
gzip $OUTDIR/pseudoGeneLinkSortFilter.bed
rm -rf $RESULTSPLIT
# Split the pseudoGeneLinkSortFilter.bed.gz BED file into one BED file per 
# chromosome in the $RESULTSPLIT directory (which is not defined). 
bedSplitOnChrom $OUTDIR/pseudoGeneLinkSortFilter.bed.gz $RESULTSPLIT
#speed up bedOverlap on chr19 (human) by splittting at the centromere
#tawk '$2 < 28500000{print $0}' chr19.bed > chr19p.bed &
#tawk '$2 >= 28500000{print $0}' chr19.bed > chr19q.bed 
#$SCRIPT/doSplit $OUTDIR
# Make $OVERLAPDIR: /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/retro/$DB/version$VERSION/run.o
mkdir -p $OVERLAPDIR

# Get list of bed files in the $RESULTSSPLIT directory and remove any for 
# random chroms and put list of resulting files in $OVERLAPDIR/results.lst
ls $RESULTSPLIT/*.bed | grep -v random > $OVERLAPDIR/results.lst
# Make a template to run bedOverlap as cluster jobs
echo "#LOOP" > $OVERLAPDIR/template
echo "${BINDIR}/bedOverlap -noBin \$(path1) {check out exists ../\$(root1)_NoOverlap.bed}" >> $OVERLAPDIR/template
echo "#ENDLOOP" >> $OVERLAPDIR/template
# Substitute the file list into the template to create a list of jobs for 
# the cluster run.
gensub2 $OVERLAPDIR/results.lst single $OVERLAPDIR/template $OVERLAPDIR/jobList

#echo clean up old files
#rm ../chr*_NoOverlap.bed
# Run jobs on the cluster using parasol.
echo "start cluster job"
echo "end of ucscRetroStep4.sh , check parasol status then run ucscRetroStep5.sh"
ssh -T $CLUSTER "cd $OVERLAPDIR ; /parasol/bin/para make jobList "
pwd
