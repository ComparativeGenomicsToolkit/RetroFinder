# STEPS FOR RUNNING THE RETROFINDER PIPELINE V2
# NOTE: do not use -skipBlatMerge which assumes this step was done manually
# and generally it isn't.
# used -skipBlatMerge which assumes this was done manually and it wasn't
# Go to your working directory
cd /hive/groups/gencode/pseudogenes/retroFinder
# Make direcotry with database and date in format YYYYMMDD e.g. hg38.20150112
# DB and DATE are defined in the DEF file
mkdir DB.DATE
cd DB.DATE
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
# TWOBIT=$GENOME/$DB/$DB.2bit (Rather than using it from /scratch/data/$DB
# ROOTDIR="~/public_html/retro/hg38Jun14"
# PDB=proteins140122 # Latest PDB in this proteins database.
# SPECIES="hg38 mm10"
# Comment out all the ARRAY variables as there is no gnfAtlas data on hg38.
# Check that latest assembly nets are being used - there are ok, netMm10,
# netCanFam3, netRheMac3. Already changed altGraphX to different table:
# ALTSPLICE=sibTxGraph - comment this out as not in hg38. 
# In DEF file, add the name of the new cluster, ku
# CLUSTER=ku
# Add and commit this DEF.hg38 file to SVN. Copy the DEF file over
# Also requires the USEALTSEQS variable for step1 script so add this line:
# USEALTSEQS=0 
# Update SCRIPT AND BINDIR as now using a branch of the SVN source tree for
# RetroFinder:
# BINDIR=/hive/users/hartera/GencodeWG/retroFinder/branches/version2/bin
# SCRIPT=/hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline
#
mkdir -p /hive/data/genomes/hg38/bed/retro
cp DEF /hive/data/genomes/hg38/bed/retro/DEF 
# once DEF file is ready to go, then start RetroFinder run. Check first that
# ucscRetroStep1.sh has correct paths for programs defined. 
# Create directory for TMPMRNA
mkdir -p /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
cp ../DEF .
# get chrom.sizes without random chroms or chrM, there are many alt loci also
# in hg38 that were not in hg19 so 285 chroms total. 
cat /hive/data/genomes/hg38/chrom.sizes | grep -v random \
   | grep -v chrUn | grep -v chrM > S1.len
# Also copy to MRNABASE:
rm -r mkdir -p /hive/data/genomes/hg38/bed/mrnaBlastz.9
mkdir -p /hive/data/genomes/hg38/bed/mrnaBlastz.9
cp S1.len /hive/data/genomes/hg38/bed/mrnaBlastz.9

# Amie Radenbaugh (Haussler lab) would like the axt files to identify the 
# positions of differences between pseudogenes and reference. Need to retain
# the axt files that are removed after running lastz. 
# Instead, comment this out in
# /hive/users/hartera/GencodeWG/retroFinder/trunk/src/pipeline/lastz.sh and
# move them to another directory. 
# # Added extra input file: AXTOUT=$6
# # at end instead of rm -f $AXT, add mv $AXT $AXTOUT
# # then as input to lastz.sh in ucscRetroStep2.sh, add $TMPMRNA/lastz/axt as
# 6th input.
# # Before the template to run lastz.sh, add this line:
# mkdir -p $TMPMRNA/lastz/axt
# Then run the first step. Added --ambiguous=iupac as flag for lastz in
# lastz.sh as there are a number of sequences with non-ACTGN characters.
# Run step 1:
screen
cd /hive/groups/gencode/pseudogenes/retroFinder/hg38.20150112/mrnaBlastz
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
