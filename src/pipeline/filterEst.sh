#!/bin/bash 
#
# Summary: Create a 2bit file of EST sequences. For Blat alignments of all ESTs, 
# filter to remove EST alignments with multiple hits that are ambigouous by 
# picking the best hits so the output is a filtered Blat aligned EST PSL file. 
# 
set -bevEu -o pipefail
source $1
echo "Starting estFilter.sh $1 for $DB"
#cd $OUTDIR
# Get a list of ids for problem sequences in Genbank
hgsql $DB -N -B -e "select acc from gbWarn" > $OUTDIR/gbWarn.id
# Get the PSL of Blat alignments for ESTs from the database and exclude those in 
# then list of problem ids. 
if [[ $SPLIT_SPLICED_EST == 1 ]] ; then
    for chr in `cut -f 1 $OUTDIR/S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$SPLICED_EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/splicedEst.psl ; done
else
    hgsql $DB -N -B -e "select * from $SPLICED_EST" |cut -f2-22 |grep -v -F -f gbWarn.id> $OUTDIR/splicedEst.psl 
fi
rm -f $OUTDIR/splicedEst.psl.gz
gzip $OUTDIR/splicedEst.psl
rm -f $OUTDIR/est.psl

if [[ $SPLIT_EST == 1 ]] ; then
    for chr in `cut -f 1 S1.len |grep -v random` ; do echo $chr ; hgsql $DB -N -B -e "select * from ${chr}_$EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/est.psl ; done
else
    hgsql $DB -N -B -e "select * from $EST" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id> $OUTDIR/est.psl 
fi
hgsql $DB -N -B -e "select * from all_est" |cut -f2-22 |grep -v -F -f $OUTDIR/gbWarn.id>> $OUTDIR/est.psl 
# Get the qNames (accessions) from the EST PSL with all ESTs and sort and uniq  
tawk '{print $10}' $OUTDIR/est.psl|sort |uniq > $OUTDIR/est.list &
# Sort the PSL of all ESTs on qName
cat $OUTDIR/est.psl |sort -k10,10> $OUTDIR/est.qName.psl
# Get PSL stats on this PSL file
pslStats -queryStats $OUTDIR/est.qName.psl $OUTDIR/est.stats
# For lines where the number of fields is > 1, select those ids with alnCount > 1 
awk 'NF > 1 && $3 > 1{print $1}' $OUTDIR/est.stats > $OUTDIR/est.doublehit.list
# For lines where the number of fields is > 1, select those ids with alnCount = 1 
awk 'NF > 1 && $3 == 1{print $1}' $OUTDIR/est.stats > $OUTDIR/est.singlehit.list
# Select the PSLs for those ids with single hits 
pslSelect -queries=$OUTDIR/est.singlehit.list $OUTDIR/est.qName.psl $OUTDIR/est.qName.single.psl &
# Get a FASTA file of EST sequences with multiple hits 
/cluster/data/genbank/bin/x86_64/gbGetSeqs -gbRoot=/cluster/data/genbank -db=$GBDB -native -inclVersion genbank est $OUTDIR/est.fa  -accFile=$OUTDIR/est.doublehit.list
# Convert the EST sequences to a 2bit file and remove ids from versions.
faToTwoBit -stripVersion $OUTDIR/est.fa $OUTDIR/est.2bit
rm -f $OUTDIR/est.fa
# Select qNames from the list of sequences with more than one hit from est.qName.psl
# and write the PSLs to est.qName.filter.psl
pslSelect -queries=$OUTDIR/est.doublehit.list $OUTDIR/est.qName.psl $OUTDIR/est.qName.filter.psl
rm -rf $OUTDIR/est
# Create est directory for output of filtering step
mkdir -p $OUTDIR/est
# Split PSL file for ids with >1 hit into files of 12000 unique qNames in the 
# est directory 
${BINDIR}/pslSplit nohead $OUTDIR/est $OUTDIR/est.qName.filter.psl -chunkSize=12000

# Make output and run directories
mkdir -p $OUTDIR/estOutput
mkdir -p $OUTDIR/estLog
mkdir -p $OUTDIR/run.est

# cd run.est
# Create a list of input PSLS for filtering
ls $OUTDIR/est/*psl > $OUTDIR/run.est/list
# Create a template for batch run of pslCDnaGenomeMatch 
echo "#LOOP" > $OUTDIR/run.est/template
echo "pslCDnaGenomeMatch \$(path1) S1.len $OUTDIR/est.2bit $TWOBIT $OUTDIR/estOutput/\$(file1).filter.psl -score=$OUTDIR/estLog/\$(file1).mrnaMatch.tab -bedOut=$OUTDIR/estLog/\$(file1).mrnaMis.bed -minDiff=4 -notAlignPenalty=3" >> $OUTDIR/run.est/template
echo "#ENDLOOP" >> $OUTDIR/run.est/template
# Substitute list of PSLs into template and then make the batch on the cluster 
# and run the jobs there. Output PSL files are in $OUTDIR/estOutput
gensub2 $OUTDIR/run.est/list single $OUTDIR/run.est/template $OUTDIR/run.est/spec
ssh $CLUSTER -T "cd $OUTDIR/run.est ; /parasol/bin/para make spec"

# Now run these commands in $OUTDIR
# Cat together the PSLs for single hits and the output of multiple hit id PSLs after
# the filtering step and cat together and sort on fields 14 and 16 (tName, tStart)
find $OUTDIR/est.qName.single.psl $OUTDIR/estOutput -name '*.psl' |xargs cat | sort -k14,14 -k16,16n > $OUTDIR/estFiltered.psl
gzip $OUTDIR/estFiltered.psl
#mkdir -p $OUTDIR/estSplit
#pslSplitOnTarget est.psl estSplit
