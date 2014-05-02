#!/bin/bash 
# this script runs the retroFinder program (pslPseudo) that looks at mRNA alignments and creates a score file.
# a cluster job is used to run the program after extracting gene predictions, repeats and other annotation from the database.
#
# est alignments are extracted at the end for use in expression analysis in further stages.
#
set -bevEu -o pipefail
source $1
echo "starting ucscRetroStep3.sh $1 for $DB"
# $OUTDIR was already created by ucscRetroStep2.sh
cp $1 $OUTDIR
# cd $OUTDIR
mkdir -p $OUTDIR/result
mkdir -p $OUTDIR/result/axt
mkdir -p $OUTDIR/log
mkdir -p $OUTDIR/out

find $NIB -name \*.nib > $OUTDIR/S1.lst
cp $GENOME/$DB/chrom.sizes $OUTDIR

#cat $RMSK/*.out |awk '{OFS="\t";print $5,$6,$7}' | grep -v position|grep -v sequence | tawk 'length($0)>2{print $0}' > rmsk.bed
if [[ -s $OUTDIR/rmsk.bed.gz ]] ; then
    echo "$OUTDIR/rmsk.bed.gz not refreshed"
else
    rm -f $OUTDIR/rmsk.bed
    if [ $RMSK == "rmsk" ]; then hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from rmsk" >> $OUTDIR/rmsk.bed ;
    else 
        for i in `cut -f 1 chrom.sizes` ;do hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from ${i}_rmsk" >> $OUTDIR/rmsk.bed ; done  
    fi
fi
if [ $NET1 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET1 where tName not like '%hap%'" > \
$OUTDIR/$NET1.txt;
fi
if [ $NET2 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET2 where tName not like '%hap%'" > \
$OUTDIR/$NET2.txt ;
fi
if [ $NET3 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET3 where tName not like '%hap%'" > \
$OUTDIR/$NET3.txt ;
fi
if [[ -s $OUTDIR/simpleRepeat.bed.gz ]] ; then
    echo "$OUTDIR/simpleRepeat.bed.gz not refreshed"
else
hgsql $DB -N -B -e "select chrom, chromStart, chromEnd from simpleRepeat" > $OUTDIR/simpleRepeat.bed ;
fi
if [ $GENE1 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 order by chrom, txStart, txEnd" | sort -k2,2 -k4,4n > $OUTDIR/$GENE1.tab;
fi
if [ $GENE2 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2" > $OUTDIR/$GENE2.tab ;
fi
if [ $GENE3 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3" > $OUTDIR/$GENE3.tab ;
fi
rm -f $OUTDIR/$GENE1.tab.gz $OUTDIR/$GENE2.tab.gz $OUTDIR/$GENE3.tab.gz $OUTDIR/*.txt.gz $OUTDIR/*.bed.gz
gzip $OUTDIR/*.tab &
gzip $OUTDIR/*.txt &
gzip $OUTDIR/*.bed &

rm -f $OUTDIR/mrna.2bit
ln $MRNABASE/mrna.2bit $OUTDIR -s
ln $MRNABASE/all_mrna.psl.gz $OUTDIR -s
ln $MRNABASE/cds.tab.gz $OUTDIR -s

####
# run retro pipeline on cluster
####
# cd $OUTDIR
mkdir -p $OUTDIR/run.0
# cd run.0
wc -l $OUTDIR/split/*.psl |grep -v total | grep -v acc.lst| sort -nr | awk '{print $2}' |sed -e 's,$OUTDIR/split/tmp,,'|sed -e 's/.psl//' > $OUTDIR/run.0/list
cp $OUTDIR/$1 $OUTDIR/run.0
echo "#LOOP" > $OUTDIR/run.0/gsub
echo "$SCRIPT/doBuildpk.sh \$(path1) $OUTDIR/$1 {check out exists $RESULT/pseudoGeneLink\$(path1).bed} " >> $OUTDIR/run.0/gsub
echo "#ENDLOOP" >> $OUTDIR/run.0/gsub

gensub2 $OUTDIR/run.0/list single $OUTDIR/run.0/gsub $OUTDIR/run.0/jobList
echo "Job Count"
wc -l $OUTDIR/run.0/jobList
ssh -T $CLUSTER "cd $OUTDIR/run.0 ; /parasol/bin/para make jobList ram=4g"
echo "check parasol status and then run ucscRetroStep4.sh DEF"
