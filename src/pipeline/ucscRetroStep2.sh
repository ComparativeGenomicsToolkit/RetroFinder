#!/bin/bash 
## Summary: The script checks that all lastz jobs were completed. The script
# removes duplicate PSL rows and puts output in pslFilter directory. It creates
# chains from the axt files, filters the chains and converts the chained
# alignments to PSL and output is one per chromosome in the psl directory.
# polyA tails are reattached to these PSLs and these are concatenated and 
# sorted on qName, tName, tStart, qStart. pslFilter filters these alignments 
# and this file is split into smaller files in a split directory, and the 
# qNames are made unique by adding suffixes e.g. -1, -2 etc. The original 
# input FASTA files (mRNA and RefSeq or other sequences) are symlinked to a 
# /gbdb/$DB/blastzRetro$VERSION directory and the sequences are loaded into 
# the database in ucscRetroSeq$VERSION and ucscRetroExtFile$VERSION tables.  
set -bevEu -o pipefail
source $1
# cd $TMPMRNA
echo "working directory is $TMPMRNA."
# Check that the number of output files is the same as the number of jobs 
# to be run - doesn't parasol check that all jobs have run and completed? 
wc -l $TMPMRNA/run.0/jobList | awk '{print $1}'> $TMPMRNA/jobs.cnt
# jobs.lst will then contain full path
find $TMPMRNA/lastz -name \*.psl| wc -l | awk '{print $1}'> $TMPMRNA/jobs.lst
echo Check Job Count from cluster run
diff $TMPMRNA/jobs.cnt $TMPMRNA/jobs.lst 
if [ $? != 0 ]; then
  echo missing jobs aborting
  wc -l $TMPMRNA/jobs.cnt $TMPMRNA/jobs.lst
  exit 3
fi
#concatenate, sort (by name, score) and de-dup psls by chrom
# Make pslFilter directory.
mkdir -p $TMPMRNA/pslFilter
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do echo $i ; cat $TMPMRNA/lastz/$i/*.psl | awk '{print $0, $1*3-$2}' | \
 sort -k 10,10 -k 22nr -T /dev/shm/retroSort | awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | \
 ${BINDIR}/pslFilterDups stdin $TMPMRNA/pslFilter/$i.psl  ; done 

#chain blocks
# Make a psl directory 
mkdir -p $TMPMRNA/psl
# For each chromosome in the S1.len file, run doChain script (created in 
# ucscRetroStep1.sh) which was created in ucscRetroStep1.sh
# this runs axtChain, then filters the chains then creates a PSL from the 
# chained alignments. The resulting PSLs are put in the psl directory, one PSL
# file per chromosome.  
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do nohup $TMPMRNA/doChain $i ; done 

#reattach polyA tail and fix alignments
# Make pslLift directoryy. For each chromosome in the S1.len do a liftup of 
# each PSL file in the psl directory. Use the mrna.lft file as to lift up
# the coordinates. Lifted PSL files are put in the pslLift directory. This 
# essentially re-attaches the polyA tails by lifting coordinates to those 
# for the sequences plus polyA tails.  
mkdir -p $TMPMRNA/pslLift
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do liftUp $TMPMRNA/pslLift/$i.psl $TMPMRNA/mrna.lft warn $TMPMRNA/psl/$i.psl -pslQ -nohead ; done

# Concatenate the lifted PSL files in to one file, mrnaBlastz.psl
# cd $TMPMRNA/pslLift
echo "MRNABASE $MRNABASE"
pslCat -nohead $TMPMRNA/pslLift/*psl > $MRNABASE/mrnaBlastz.psl

# cd $MRNABASE
# Sort the PSL file by qName, tName, tStart, qStart
sort -k10,10 -k14,14 -k16,16n -k12,12n $MRNABASE/mrnaBlastz.psl > $MRNABASE/mrnaBlastz.sort.psl
#pslCDnaFilter -minCover=0.05 -minId=0.65 mrnaBlastz.sort.psl mrnaBlastz.65.psl
##                       seqs    aligns
#             total:     223377  11187645
#     drop minIdent:     103075  6082119
#     drop minCover:     51568   457653
#              kept:     222565  4647873
#pslCDnaFilter -minCover=0.05 -minId=0.62 mrnaBlastz.sort.psl mrnaBlastz.psl
#                        seqs    aligns
#             total:     223377  11187645
#     drop minIdent:     74887   3622214
#     drop minCover:     52526   567156
#              kept:     223112  6998275

# Then filter PSLs for minCover =0.05 and minimum identity = 0.58 and then 
# gzip the sorted PSL file (unfiltered).
pslCDnaFilter -minCover=0.05 -minId=0.58 $MRNABASE/mrnaBlastz.sort.psl $MRNABASE/mrnaBlastz.psl
rm -f $MRNABASE/mrnaBlastz.sort.psl.gz
gzip $MRNABASE/mrnaBlastz.sort.psl &

#split for pipeline cluster run 
# Create a split directory and split the PSLs with no header into the split 
# directory where there are 120 unique qNames per file.
mkdir -p $OUTDIR/split
${BINDIR}/pslSplit nohead $OUTDIR/split $MRNABASE/mrnaBlastz.psl -chunkSize=120

# cd $OUTDIR/split
# For each of the split PSL files, run pslQueryUniq and put output in temp.psl
# which is then renamed with original file name. pslQueryUniq makes the qNames
# unique by adding a dash and a number as a suffix e.g. -1, -2 etc. 
for i in `ls $OUTDIR/split/tmp*.psl` ; do $SCRIPT/pslQueryUniq $i > $OUTDIR/split/temp.psl ; mv $OUTDIR/split/temp.psl $i ; done
# For each tmp PSL file in the split directory, get the accession and file name
# and make an acc.lst file of the file name and accession (qName) for each file 
# and accession in those files. 
grep chr $OUTDIR/split/tmp* | awk '{print $1,$10}' | awk -F":" '{print $1,$2}'|awk '{print $1,$3}'| sort -u > $OUTDIR/split/acc.lst

#export TMPDIR=/scratch/tmp
# TMPDIR is an environment variable that points to /data/tmp
echo TMPDIR
#load mrna alignment track into browser, load into mrnaBlastz table with 
# the qName version numbers stripped off.
awk -f $SCRIPT/stripversion.awk $MRNABASE/mrnaBlastz.psl | hgLoadPsl $DB stdin -table=mrnaBlastz

echo "ucscRetroStep2.sh mrna alignments complete"
# cd $MRNABASE
cp $MRNABASE/S1.len $OUTDIR
# cd $OUTDIR
# pwd

#load mrna sequences into browser (with version numbers)
# Symlink the FASTA file of all mRNA input sequences to gbdb
mkdir -p /gbdb/$DB/blastzRetro${VERSION}
rm -f /gbdb/$DB/blastzRetro${VERSION}/mrna.fa 
ln $MRNABASE/mrna.fa /gbdb/$DB/blastzRetro${VERSION}/ -s
# Load sequences into the database and replace existing sequences with the 
# same id. Tables created are ucscRetroSeq$VERSION and ucscRetroExtFile$VERSION
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro${VERSION}/mrna.fa  -seqTbl=ucscRetroSeq${VERSION} -extFileTbl=ucscRetroExtFile${VERSION}
rm -f /gbdb/$DB/blastzRetro${VERSION}/refseq.fa
# Symlink FASTA file of RefSeq input sequences to gbdb 
ln $MRNABASE/refseq.fa /gbdb/$DB/blastzRetro${VERSION}/ -s
# Then add the RefSeq sequences to the database tables.
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro${VERSION}/refseq.fa  -seqTbl=ucscRetroSeq${VERSION} -extFileTbl=ucscRetroExtFile${VERSION}

echo "run ucscRetroStep3.sh DEF to run retro pipeline"
