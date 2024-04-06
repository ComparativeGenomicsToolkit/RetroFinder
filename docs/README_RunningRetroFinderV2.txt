# STEPS FOR RUNNING THE RETROFINDER PIPELINE V2 - 2015-02-10, Rachel Harte
# For more details on RetroFinder and what each step does, see 
# retroFinder/trunk/docs/README_AboutRetroFinderV2.txt
#
# NOTE: do not use -skipBlatMerge which assumes this step was done manually
# and generally it isn't.
# Variables below are defined in the DEF file.
# Go to your working directory

    cd /hive/groups/gencode/pseudogenes/retroFinder

# Make direcotry with database and date in format YYYYMMDD e.g. hg38.20150112
# DB and DATE are defined in the DEF file

    mkdir $DB/$DATE
    cd $DB/$DATE

# Copy the DEF file here, there are sample DEF files used for RetroFinder run
# with different assemblies in the RetroFinder SVN source tree. 
# Also create a README.txt file like this one to document the steps done and 
# problems encountered and how they were resolved.

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

   mkdir -p /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
   cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
   cp ../DEF .

# Get the chrom.sizes without random chroms or chrM, there are many 
# alt loci also in the current human assembly hg (GRCh38) that were not in hg19 # so include those chroms (285 chroms total). 
# Get chromosome sizes and remove randoms, chrUn and chrM. Can also query the
# database to get these:
# hgsql -Ne 'select chrom, size from chromInfo where chrom not like "%random%"
# or chrom not like "chrUn%" or chrom != "chrM";' $DB

    cat /hive/data/genomes/$DB/chrom.sizes | grep -v random | grep -v chrUn | grep -v chrM > S1.len

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
# Using --ambiguous=iupac as a flag to the lastz program in lastz.sh 
# as there are a number of sequences with non-ACTGN characters.
# Run step 1:

    screen

# cd to $RUNDIR/mrnaBlastz
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep1.sh DEF >& step1hg38.log &

# Check step1 cluster jobs:
    ssh ku
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz/hg38/run.0
    para check

# If jobs all ran ok, then run step2. 
    ssh hgwdev
    screen 

    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep2.sh DEF >& step2hg38.log &
# Start step3
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep3.sh DEF >& step3hg38.log &

# Check step3 cluster jobs
    ssh ku
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB/run.0
    para check

# If jobs ran ok, then run Step 4:
# cd to 
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep4.sh DEF >& step4hg38.log &

# Check step4 cluster jobs
    ssh ku

# cd to $OUTDIR/run.o
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB/run.o
    para check

# If cluster jobs ran ok, then run step 5
    ssh hgwdev

# This program needs to be built from the kent source tree:
    cd ~/kent/src/hg/pslCDnaGenomeMatch; make
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/mrnaBlastz
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/ucscRetroStep5.sh DEF >& step5hg38.log &

# Get expression information. cd to $OUTDIR
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/filterMrna.sh DEF >& filterMrna.log &
    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/filterEst.sh DEF >& filterEst.log &

# Check the filterEst.sh cluster run
    ssh ku
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB/run.est
    para check
    
# If the filterEst.sh cluster jobs ran ok, run analyseExpress.sh to update the
# expression levels of retros
    ssh hgwdev
# cd to $OUTDIR
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB

    time /hive/users/hartera/GencodeWG/retroFinder/branches/version2/src/pipeline/analyseExpress.sh DEF >& analyseExpress.log &

# Check the analseExpress.sh cluster jobs
    ssh ku

# cd to $OUTDIR/estSplit
    cd /hive/groups/gencode/pseudogenes/retroFinder/$DB/$DATE/retro/version$VERSION/$DB/estSplit
    para check

# If all cluster jobs ran ok, then can add the track to the UCSC Genome
# Browser
# Add $OUTDIR/trackDb.retro entry to 
# ~/kent/src/hg/makeDb/trackDb/<organism>/$DB/trackDb.ra.
#
# Tables built are:
# ucscRetroAli
# ucscRetroCds
# ucscRetroCount
# ucscRetroExpressed
# ucscRetroExtFile
# ucscRetroInfo
# ucscRetroOrtho
# ucscRetroSeq
# Version number is appended to these table names. 
# The sequence files are in:
# /gbdb/$DB/blastzRetro$VERSION
# symlinks point to $MRNABASE/mrna.fa and $MRNABASE/refseq.fa   
