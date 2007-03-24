# run from kkstore02
DB=hg18
BASE=/cluster/data/hg18/bed/pseudo
#OUTDIR=/cluster/bluearc/$DB/pseudo
#OUTDIR=/cluster/panasas/home/store/pseudo/$DB/
OUTDIR=/san/sanvol1/scratch/pseudo/$DB/    
RESULT=$OUTDIR/result
OUT=/cluster/bluearc/$DB/pseudo/out
echo catting output
wc -l run.1/spec | awk '{print $1}'> jobs.cnt
pushd $OUTDIR/result ; ls pseudoGeneLink[0-9]*.bed | wc -l | awk '{print $1}'> $BASE/jobs.lst
popd
echo Check Count 
diff jobs.cnt jobs.lst 
if [ $? != 0 ]; then
  echo missing jobs aborting
  exit 3
fi

echo catting Sorting and Filtering pseudoGeneLinkSortFilter.bed
#pushd $OUTDIR/result ; cat pseudoGeneLink[0-9]*.bed | tawk '$5 > 10 && $15 > 10000 && $35 > 650 {OFS="\t";print $0}'| sort -k1,1 -k2,3n -k4,4 > $BASE/pseudoGeneLinkSortFilter.bed ; #/bin/rm $OUTDIR/result/pseudoGeneLink[0-9]*.bed
echo "pushd $OUTDIR/result ; cat pseudoGeneLink[0-9]*.bed | tawk 'col5 > 10 && (col14 > 10000 || col14 == -1) && col35 > 620 {OFS="\t";print $0}'| sort -k1,1 -k2,3n -k4,4 to  $BASE/pseudoGeneLinkSortFilter.bed; "
pushd $OUTDIR/result ; cat pseudoGeneLink[0-9]*.bed | tawk '$5 > 10 && ($14 > 10000 || $14 == -1) && $35 > 620 {OFS="\t";print $0}'| sort -k1,1 -k2,3n -k4,4 -T /scratch > $BASE/pseudoGeneLinkSortFilter.bed; #/bin/rm $OUTDIR/result/pseudoGeneLink[0-9]*.bed
popd
wc -l pseudoGeneLinkSortFilter.bed
pushd $OUT
echo creating pseudo.psl
cat pseudo[0-9]*.psl > $BASE/pseudo.psl & #;/bin/rm $OUT/pseudo[0-9]*.psl &
popd
echo Removing Overlaps
echo splitting
doSplit $OUTDIR
rm -f pseudoGeneLinkSortFilter.bed.gz
gzip pseudoGeneLinkSortFilter.bed &
pushd $OUTDIR
cd run.o 
echo clean up old files
rm ../chr*pseudoNoOverlap.bed
echo "start cluster job"
ssh -T pk "cd $OUTDIR/run.o ; para make spec.check  "
#echo start spec.overlap
#spec.overlap
popd
pwd
echo buildSort.part2.sh $OUTDIR/result $DB $BASE
buildSort.part2.sh $OUTDIR/result $DB $BASE
#buildSort.part2.sh /san/sanvol1/scratch/pseudo/hg18/result/ hg18 /cluster/data/hg18/bed/pseudo/

