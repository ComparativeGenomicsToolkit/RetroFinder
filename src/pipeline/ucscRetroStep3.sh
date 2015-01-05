#!/bin/bash 
# This script runs the retroFinder program (pslPseudo) that looks at mRNA 
# alignments and creates a score file.
# A cluster job is used to run the program after extracting gene predictions, 
# repeats and other annotations from the database.
#
## Summary: Make output directories for pslPseudo, get input annotation files 
# for pslPseudo (RepeatMasker, simple repeats, nets, gene annotations) then 
# create list of jobs to run doBuildpk.sh for split PSL files in $OUTDIR/split 
# directory. doBuild.pk runs pslPseudo which looks at the alignments and 
# annotations and scores the retrogene candidates; it compiles the information
# for the ucscRetroInfoXX table. Run the jobs on the cluster.  

set -bevEu -o pipefail
source $1
echo "starting ucscRetroStep3.sh $1 for $DB"
# $OUTDIR was already created by ucscRetroStep2.sh
cp $1 $OUTDIR
# cd $OUTDIR
# Create result, axt, log and out directories.
mkdir -p $OUTDIR/result
mkdir -p $OUTDIR/result/axt
mkdir -p $OUTDIR/log
mkdir -p $OUTDIR/out

# Now using genome 2bit file as input to pslPseudo so no longer need to 
# create a list of nibs as before. 

# Copy file of chromosome sizes here - isn't there already a file of chrom
# sizes created before the script is run? 
cp $GENOME/$DB/chrom.sizes $OUTDIR

#cat $RMSK/*.out |awk '{OFS="\t";print $5,$6,$7}' | grep -v position|grep -v sequence | tawk 'length($0)>2{print $0}' > rmsk.bed
# Get RepeatMasker repeat alignments from the database
if [[ -s $OUTDIR/rmsk.bed.gz ]] ; then
    echo "$OUTDIR/rmsk.bed.gz not refreshed"
else
    rm -f $OUTDIR/rmsk.bed
    if [ $RMSK == "rmsk" ]; then hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from rmsk" >> $OUTDIR/rmsk.bed ;
    else 
        for i in `cut -f 1 chrom.sizes` ;do hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from ${i}_rmsk" >> $OUTDIR/rmsk.bed ; done  
    fi
fi

# Get the net alignments for three genomes specified in the DEF file. 
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
# Get simple repeat annotations (from TRF run) from the database.
if [[ -s $OUTDIR/simpleRepeat.bed.gz ]] ; then
    echo "$OUTDIR/simpleRepeat.bed.gz not refreshed"
else
hgsql $DB -N -B -e "select chrom, chromStart, chromEnd from simpleRepeat" > $OUTDIR/simpleRepeat.bed ;
fi
# Get gene annotations from the database for the three genesets specified 
# in the DEF file. 
if [ $GENE1 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 order by chrom, txStart, txEnd" | sort -k2,2 -k4,4n > $OUTDIR/$GENE1.tab;
fi
if [ $GENE2 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2" > $OUTDIR/$GENE2.tab ;
fi
if [ $GENE3 != "/dev/null" ]; 
then hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3" > $OUTDIR/$GENE3.tab ;
fi
# Remove previous gzipped files and then gzip these annotation files. 
rm -f $OUTDIR/$GENE1.tab.gz $OUTDIR/$GENE2.tab.gz $OUTDIR/$GENE3.tab.gz $OUTDIR/*.txt.gz $OUTDIR/*.bed.gz
gzip $OUTDIR/*.tab &
gzip $OUTDIR/*.txt &
gzip $OUTDIR/*.bed &

# Remove mRNA 2bit file and symlink mrna.2bit file here. Symlink all_mrna PSL
# file (from ucscRetroStep1.sh) and cds.tab file of CDS regions of 
# Genbank sequences or other input sequences (from ucscRetroStep1.sh). 
rm -f $OUTDIR/mrna.2bit
ln $MRNABASE/mrna.2bit $OUTDIR -s
ln $MRNABASE/all_mrna.psl.gz $OUTDIR -s
ln $MRNABASE/cds.tab.gz $OUTDIR -s

####
# run retro pipeline on cluster
####
# cd $OUTDIR
# Make run.0 directory for cluster run of pslPseudo
mkdir -p $OUTDIR/run.0
# cd run.0
# Create a list of the numbers after tmp and before .psl for the file 
# names in the $OUTDIR/split directory.
wc -l $OUTDIR/split/*.psl |grep -v total | grep -v $OUTDIR/split/acc.lst| sort -nr | awk '{print $2}' | sed -e "s,$OUTDIR/split/tmp,," | sed -e 's/.psl//' > $OUTDIR/run.0/list
cp $OUTDIR/$1 $OUTDIR/run.0
echo "#LOOP" > $OUTDIR/run.0/gsub
echo "$SCRIPT/doBuildpk.sh \$(path1) $OUTDIR/$1 {check out exists $RESULT/pseudoGeneLink\$(path1).bed} " >> $OUTDIR/run.0/gsub
echo "#ENDLOOP" >> $OUTDIR/run.0/gsub

# Substitute in the numbers from list into gsub to make a list of jobs
gensub2 $OUTDIR/run.0/list single $OUTDIR/run.0/gsub $OUTDIR/run.0/jobList
# Get job count
echo "Job Count"
wc -l $OUTDIR/run.0/jobList
# Run jobs on the cluster using parasol and 4Gb of RAM
ssh -T $CLUSTER "cd $OUTDIR/run.0 ; /parasol/bin/para make jobList ram=4g"
echo "check parasol status and then run ucscRetroStep4.sh DEF"
