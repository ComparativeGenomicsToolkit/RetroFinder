#!/bin/bash 
#
# runs small cluster job to remove overlaps by picking highest scoring retro in each locus
#
set -bevEu -o pipefail
source $1
echo "---------------------------------------------------------------------------------------"
echo "Starting ucscRetroStep4.sh $1 on $HOST- catting output from retro pipeline cluster run"
echo "---------------------------------------------------------------------------------------"
# cd $OUTDIR
if [[ -s $OUTDIR/S1.len ]] ; then
    echo "$OUTDIR/S1.len exists with `wc -l S1.len` rows"
else
    echo "Please create $OUTDIR/S1.len from chrom.sizes without random chroms or chrM."
    exit 3
fi
wc -l $OUTDIR/run.0/jobList | awk '{print $1}'> $OUTDIR/jobs.cnt
pushd $RESULT ; ls pseudoGeneLink[0-9]*.bed | wc -l | awk '{print $1}'> $OUTDIR/jobs.lst
popd
#echo Check Job Count from cluster run
#diff jobs.cnt jobs.lst 
#if [ $? != 0 ]; then
#  echo missing jobs aborting
#  exit 3
#fi

# filter out low scoring hits as defined by retroscore and alignment score to parent gene
rm -f $OUTDIR/pseudoGeneLinkSortFilter.bed.gz
echo catting Sorting and Filtering pseudoGeneLinkSortFilter.bed
# pushd $RESULT ; 

cat $RESULT/pseudoGeneLink[0-9]*.bed | tawk '$5 > 300 && ($14 > 10000 || $14 == -1) {OFS="\t";print $0}'| $SCRIPT/removeTandemDups | sort -k1,1 -k2,3n -k4,4 -T /scratch > $OUTDIR/pseudoGeneLinkSortFilter.bed; #/bin/rm $RESULT/pseudoGeneLink[0-9]*.bed
# popd
wc -l $OUTDIR/pseudoGeneLinkSortFilter.bed
# pushd $OUT
cat $OUT/pseudo[0-9]*.psl > $OUTDIR/pseudo.psl #;/bin/rm $OUT/pseudo[0-9]*.psl &
cat $OUT/ortho*.txt | sed -e 's/.txt//'> $OUTDIR/ortho.txt
# popd
echo Removing Overlaps
gzip $OUTDIR/pseudoGeneLinkSortFilter.bed
rm -rf $RESULTSPLIT
bedSplitOnChrom $OUTDIR/pseudoGeneLinkSortFilter.bed.gz $RESULTSPLIT
#speed up bedOverlap on chr19 (human) by splittting at the centromere
#tawk '$2 < 28500000{print $0}' chr19.bed > chr19p.bed &
#tawk '$2 >= 28500000{print $0}' chr19.bed > chr19q.bed 
#$SCRIPT/doSplit $OUTDIR
mkdir -p $OVERLAPDIR

ls $RESULTSPLIT/*.bed | grep -v random > $OVERLAPDIR/results.lst
echo "#LOOP" > $OVERLAPDIR/template
echo "${BINDIR}/bedOverlap -noBin \$(path1) {check out exists ../\$(root1)_NoOverlap.bed}" >> $OVERLAPDIR/template
echo "#ENDLOOP" >> $OVERLAPDIR/template
gensub2 $OVERLAPDIR/results.lst single $OVERLAPDIR/template $OVERLAPDIR/jobList

#echo clean up old files
#rm ../chr*_NoOverlap.bed
echo "start cluster job"
echo "end of ucscRetroStep4.sh , check parasol status then run ucscRetroStep5.sh"
ssh -T $CLUSTER "cd $OVERLAPDIR ; /parasol/bin/para make jobList "
pwd
