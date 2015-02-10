# STEPS FOR RUNNING THE RETROFINDER PIPELINE V2 - 2015-02-10, Rachel Harte
# NOTE: do not use -skipBlatMerge which assumes this step was done manually
# and generally it isn't.
# Variables below are defined in the DEF file.
# Go to your working directory
cd /hive/groups/gencode/pseudogenes/retroFinder
# Make direcotry with database and date in format YYYYMMDD e.g. hg38.20150112
# DB and DATE are defined in the DEF file
mkdir $DB.$DATE
cd $DB.$DATE
# Copy the DEF file here, there are sample DEF files used for RetroFinder run
# with different assemblies in the RetroFinder SVN source tree. 
# Also create a README.txt file like this one to document the steps done and 
# problems encountered and how they were resolved.
chmod +x DEF

# Edit DEF file to change relevant variables:
# DATE=YYYYMMDD
# RUNDATE="YYYY-MM-DD"
# DB=<assembly name>
# Version is incremented each run for each assembly, 
# VERSION=9
# ROOTDIR="~/public_html/retro/hg38Jun14" # use database, month and year
# PDB=proteins140122 # Latest PDB in this proteins database.
# SPECIES="hg38 mm10"
# This flags that sequences other than GenBank or RefSeq are not being used
# USEALTSEQS=0 
# Should be set to the branch of the source tree
# BINDIR=/hive/users/hartera/GencodeWG/retroFinder/branches/version2/bin
# SCRIPT=/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline
# Make $RETRODIR:
mkdir -p /hive/data/genomes/$DB/bed/retro
cp DEF /hive/data/genomes/$DB/bed/retro/DEF 
# Once DEF file is ready to go, then prepare directories and data
# required for the RetroFinder run. 
# Create the $TMPMRNA directory:
mkdir -p /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/mrnaBlastz
cd /hive/groups/gencode/pseudogenes/retroFinder/$DB.$DATE/mrnaBlastz
cp ../DEF .
# Get the chrom.sizes without random chroms or chrM, there are many 
# alt loci also in the current human assembly hg (GRCh38) that were not in hg19 # so include those chroms (285 chroms total). 
# Get chromosome sizes and remove randoms, chrUn and chrM. Can also query the
# database to get these:
# hgsql -Ne 'select chrom, size from chromInfo where chrom not like "%random%"
# or chrom not like "chrUn%" or chrom != "chrM";' $DB
cat /hive/data/genomes/$DB/chrom.sizes | grep -v random \
   | grep -v chrUn | grep -v chrM > S1.len
# Also copy the chromosome sizes file to MRNABASE:
mkdir -p /hive/data/genomes/$DB/bed/mrnaBlastz.$VERSION
mkdir -p /hive/data/genomes/$DB/bed/mrnaBlastz.$VERSION
cp S1.len /hive/data/genomes/$DB/bed/mrnaBlastz.$VERSION

# A grad student once asked for the axt files to identify the 
# positions of differences between pseudogenes and reference. Therefore now 
# retain the axt files that used to be removed after running lastz. 
# This line was commented out in
# /trunk/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/lastz.sh and move them to another directory. 
# The axt files are put in $TMPMRNA/lastz/axt 
# # Before the template to run lastz.sh, add this line:
# Added --ambiguous=iupac as a flag to the lastz program in lastz.sh 
# as there are a number of sequences with non-ACTGN characters.
# Run step 1:
screen
cd /hive/groups/gencode/pseudogenes/retroFinder/$DB.DATE/mrnaBlastz
time \
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep1.sh DEF >& step1hg38.log &
# Check cluster jobs:
ssh ku
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz/hg38/run.0
/parasol/bin/para check
# Cluster jobs were crashing. The problem was there there were a Q in one of
# the RefSeq sequences and this is not an IUPAC code. There were also extra
# characters in the sequence that look like a binary file, I have seen this
# error before with gbGetSeqs. I reproduced this twice more in the same
# sequence when fetching the RefSeqs. Corrected this sequence in refseq.fa and
# continued running ucscRetroStep1.sh from the point where 
# # Concatenate the mRNA and RefSeq sequence files.
# cat $MRNABASE/mrna.fa $MRNABASE/refseq.fa > $MRNABASE/raw.fa
time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/step1.sh DEF >>& step1hg38.log &
105.623u 7.780s 6:15:49.04 0.5% 0+0k 0+936io 0pf+0w
Tue Jan 13 09:02:30 PST 2015
# Check cluster jobs:
ssh ku
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz/hg38/run.0
para check
# Jobs all ran ok. Run step2. 
ssh hgwdev
screen 
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep2.sh DEF \
    >& step2hg38.log &
# 122156.444u 49512.146s 46:45:52.26 101.9%       0+0k 3912+201408928io 4pf+0w
Thu Jan 15 09:39:36 PST 2015
# step2 has finished ok so start step3.
# in doBuildpk.sh, use $TMPDIR as root for temp directory, this is an 
# environment variable defined as /data/tmp on hgwdev but /scratch/tmp on the
# ku cluster so TMP=$TMPDIR/retro.$$.$USER
ssh hgwdev
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep3.sh DEF \
    >& step3hg38.log &
ssh ku
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38/run.0
/parasol/bin/para check
# Now run Step 4:
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
time
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep4.sh
 DEF >& step4hg38.log &

# Check cluster run finished ok. 
ssh ku
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38/run.o
/parasol/bin/para check
# All 243 jobs ran ok. 
# Run step 5
ssh hgwdev
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep5.sh DEF \
    >& step5hg38.log &
1147.857u 10.285s 19:08.38 100.8%       0+0k 6488+33128io 7pf+0w
# This program needs to be built from the kent source tree:
cd ~/kent/src/hg/pslCDnaGenomeMatch; make
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38
time
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/filterMrna.sh \
DEF >& filterMrna.log &
Fri Jan 16 09:34:42 PST 2015
# This ran fine.
# in filterEst.sh, also add this path for pslCDnaGenomeMatch
time
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/filterEst.sh \
DEF >& filterEst.log &
# 6265.105u 206.318s 2:30:27.65 71.6%     0+0k 0+18864192io 0pf+0w
# This ran ok. Checked final output and cluster run. 
ssh ku
cd
/hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38/run.est
para check
# jobs ran ok.
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38
# Then need to run analyseExpress.sh to update expression levels of retros
time 
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/analyseExpress.sh
\ DEF >& analyseExpress.log &
Tue Jan 20 15:51:20 PST 2015
# Seems that the analyseExpress.sh script has stalled so check the cluster
# jobs.
ssh ku
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38/estSplit
# These all ran ok. Just seems to have stalled so just run the rest of the 
# script and kill this instance of it running.
ssh hgwdev 
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/retro/version9/hg38
/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/analyseExp.sh DEF >>& analyseExpress.log &
# Crashed on all joins, says the files are not sorted even though they 
# have been. That happened earlier too. Fix by adding sort -k1,1. 
# Previously, in analyseExpress.sh, added option to pslSplitOnTarget to 
# increase the maximum allowed targets to 500. There are 455 for hg38.
# Add trackDb.retro entry to ~/kent/src/hg/makeDb/trackDb/human/hg38/trackDb.ra.  
