#!/bin/bash 
# post process retro scores
set -beEu -o pipefail
source $1
echo "---------------------------------------------------------------------------------------"
echo "Starting ucscRetroStep4.sh $1 on $HOST- catting output from retro pipeline cluster run"
echo "---------------------------------------------------------------------------------------"
wc -l run.0/jobList | awk '{print $1}'> jobs.cnt
pushd $RESULT ; ls pseudoGeneLink[0-9]*.bed | wc -l | awk '{print $1}'> $OUTDIR/jobs.lst
popd
echo Check Job Count from cluster run
diff jobs.cnt jobs.lst 
if [ $? != 0 ]; then
  echo missing jobs aborting
  exit 3
fi
cp $GENOME/$DB/chrom.sizes $OUTDIR/S1.len

rm -f pseudoGeneLinkSortFilter.bed.gz
echo catting Sorting and Filtering pseudoGeneLinkSortFilter.bed
echo "pushd $RESULT ; cat pseudoGeneLink[0-9]*.bed | tawk 'col5 > 10 && (col14 > 10000 || col14 == -1) && col35 > 620 {OFS="\t";print $0}'| sort -k1,1 -k2,3n -k4,4 to  $OUTDIR/pseudoGeneLinkSortFilter.bed; "
pushd $RESULT ; cat pseudoGeneLink[0-9]*.bed | tawk '$5 > 10 && ($14 > 10000 || $14 == -1) && $35 > 620 {OFS="\t";print $0}'| sort -k1,1 -k2,3n -k4,4 -T /scratch > $OUTDIR/pseudoGeneLinkSortFilter.bed; #/bin/rm $RESULT/pseudoGeneLink[0-9]*.bed
popd
wc -l pseudoGeneLinkSortFilter.bed
pushd $OUT
echo creating pseudo.psl
cat pseudo[0-9]*.psl > $OUTDIR/pseudo.psl & #;/bin/rm $OUT/pseudo[0-9]*.psl &
popd
echo Removing Overlaps
echo splitting
RESULTSPLIT=$OUTDIR/resultSplit
gzip pseudoGeneLinkSortFilter.bed
rm -rf $RESULTSPLIT
bedSplitOnChrom pseudoGeneLinkSortFilter.bed.gz $RESULTSPLIT
#speed up bedOverlap on chr19 (human) by splittting at the centromere
#tawk '$2 < 28500000{print $0}' chr19.bed > chr19p.bed &
#tawk '$2 >= 28500000{print $0}' chr19.bed > chr19q.bed 
#$SCRIPT/doSplit $OUTDIR
pushd $OUTDIR
mkdir -p $OVERLAPDIR
cd $OVERLAPDIR

ls $RESULTSPLIT/*.bed | grep -v random > results.lst
echo "#LOOP" > template
echo "bedOverlap -noBin \$(path1) {check out exists ../\$(root1)_NoOverlap.bed}" >> template
echo "#ENDLOOP" >> template
gensub2 results.lst single template jobList

echo clean up old files
rm ../chr*_NoOverlap.bed
echo "start cluster job"
echo "end of ucscRetroStep4.sh , check parasol status then run ucscRetroStep5.sh"
ssh -T $CLUSTER "cd $OVERLAPDIR ; para make jobList "
popd
pwd
$SCRIPT/ucscRetroStep5.sh $1
