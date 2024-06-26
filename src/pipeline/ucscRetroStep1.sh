#!/bin/bash 
### Summary:
# Gets input sequences, and PSLs, removes the polyA tails from the sequences, 
# aligns the trimmed sequences to the genome using lastz. 
set -bveEu -o pipefail
if [ $BASH_ARGC ]; then
       echo "Parameter file is $1"
else
       echo "usage $0 DEF"
       exit 3
fi
source $1
mkdir -p $OUTDIR
cp -f $1 $OUTDIR
echo running mrna alignment against $DB
#align mrnas to $DB using lastz
mkdir -p $TMPMRNA
mkdir -p $MRNABASE
cd $MRNABASE

if [[ -s $MRNABASE/S1.len ]] ; then
    echo "$MRNABASE/S1.len exists with `wc -l S1.len` rows"
else
    echo "Fatal error: Please create $MRNABASE/S1.len from chrom.sizes without random chroms or chrM, and rerun."
    exit 3
fi
#Pull the latest mRNA and RefSeq sequences from genbank or 
# use an alternative sequence source such as Ensembl.
if [[ -s $MRNABASE/mrna.fa ]] ; then
    echo "$MRNABASE/mrna.fa exists, extraction from genbank skipped"
elif [[ $USEALTSEQS == 1 ]] && [[ $ALTSEQSOURCE = "ensGene" ]]; then
    echo "Using Ensembl transcript sequences, extract from Ensembl database"
    echo "$SCRIPT/getEnsemblTxSeqsWithVersions.pl --all $GENOMENAME all $MRNABASE/ens.$DB.fa"
    $SCRIPT/getEnsemblTxSeqsWithVersions.pl --all $ORGNAME all $MRNABASE/ens.$DB.fa
    ln -s $MRNABASE/mrna.fa $MRNABASE/ens.$DB.fa
else
    echo "extracting mRNA and refseq from Genbank "
    echo "change non standard lower case characters to N's being careful not to change headers that are upper case"
# Change some non-standard characters to Ns e.g. P and I and convert sequence
# to upper case if not already. 
    /cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
       genbank mrna stdout | awk '/^>/{print $0}!/^>/{$1=toupper($0);s1=gsub("P","N");s2=gsub("I","N",$s1);print $s2}' > tmp.fa ; mv tmp.fa $MRNABASE/mrna.fa 
fi
if [[ -s $MRNABASE/refseq.fa ]] ; then
    echo "$MRNABASE/refseq.fa exists, extraction from genbank skipped"
else
    /cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
       refSeq mrna stdout | awk '/^>/{print $0}!/^>/{$1=toupper($0);s1=gsub("P","N");s2=gsub("I","N",$s1);print $s2}' > $MRNABASE/refseq.fa 
fi
# Get PSL files for mRNAs and RefSeqs or for alternative sequence source. 
if [[ -s $MRNABASE/all_mrna.psl.gz ]] ; then
    echo "$MRNABASE/all_mrna.psl.gz not refreshed, must be extracted at the same time as mrna and refseq sequences. "
else
    if [[ $USEALTSEQS == 1 ]]; then
        hgsql $DB -N -B -e "select * from $ALTSEQSOURCE;" | cut -f 2-> $MRNABASE/$ALTSEQSOURCE.$DB.gp
        genePredToPsl /hive/data/genomes/$DB/chrom.sizes $MRNABASE/$ALTSEQSOURCE.$DB.gp $MRNABASE/all_mrna.psl
    else
        hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from all_mrna" > $MRNABASE/all_mrna.psl 
    fi
    hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from refSeqAli" >> $MRNABASE/all_mrna.psl
    rm -f $MRNABASE/all_mrna.psl.gz
    gzip $MRNABASE/all_mrna.psl
fi
# Create a file of CDS regions for sequences. 
if [[ -s $MRNABASE/cds.tab.gz ]] ; then
    echo "$MRNABASE/cds.tab.gz not refreshed, must be extracted at the same time as all_mrna.psl.gz, mRNA and RefSeq sequences. "
else
    hgsql $DB -N -B -e "select acc, version, name, type from refSeqAli a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" > $MRNABASE/cds.tab
    if [[ $USEALTSEQS == 1 ]]; then
        hgsql $DB -N -B -e "select name, cdsStart, cdsEnd from $ALTSEQSOURCE;" \
           > $MRNABASE/$DB.$ALTSEQSOURCE.cds.txt
        tawk '{print $1, "", $2+1, $3}' $MRNABASE/$DB.$ALTSEQSOURCE.cds.txt \
           >> $MRNABASE/cds.tab
    else
    hgsql $DB -N -B -e "select acc, version, name, type from all_mrna a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" >> $MRNABASE/cds.tab
    fi
    gzip $MRNABASE/cds.tab
fi
# For the filterMrna.sh step, need a PSL file of only GenBank mRNAs
# GenBank tables are updated daily so all GenBank data should be collected
# at the same time otherwise files are out of sync.
if [[ -s $MRNABASE/gbMrnaOnly.psl ]] ; then
    echo "$MRNABASE/gbMrnaOnly.psl not refreshed, must be extracted at the same time as all_mrna.psl.gz, cds.tab.gz, mRNA and refSeq sequences."
else
    if [[ $USEALTSEQS == 1 ]]; then
       cp $MRNABASE/all_mrna.psl.gz $MRNABASE/gbMrnaOnly.psl.gz
    else 
       hgsql $DB -N -B -e "select * from all_mrna" | cut -f2-22 > $MRNABASE/gbMrnaOnly.psl
    fi
fi
# Concatenate the mRNA and RefSeq sequence files.
cat $MRNABASE/mrna.fa $MRNABASE/refseq.fa > $MRNABASE/raw.fa

# Remove polyA tail before aligning
faTrimPolyA $MRNABASE/raw.fa $MRNABASE/trim.fa 
mkdir -p $TMPMRNA 
cp -p $MRNABASE/trim.fa $TMPMRNA 

pwd

# Convert the mRNA and RefSeq sequences file (or other input sequence file) to
# 2bit format. 
faToTwoBit $MRNABASE/raw.fa $MRNABASE/mrna.2bit
# Then convert back to FASTA then back to twoBit format stripping version 
# number. 
twoBitToFa $MRNABASE/mrna.2bit stdout | faToTwoBit stdin -stripVersion $MRNABASE/mrnaNoversion.2bit

# Get sequence sizes for polyA trimmed sequences and sequences with polyA tails.
# Copy trim.len to S2.len.
# Then sort these files by sequence id and join them together - both.len.
faSize $MRNABASE/trim.fa -detailed > $MRNABASE/trim.len
cp $MRNABASE/trim.len $MRNABASE/S2.len
faSize $MRNABASE/raw.fa -detailed > $MRNABASE/raw.len
sort $MRNABASE/raw.len > $MRNABASE/x
mv $MRNABASE/x $MRNABASE/raw.len
sort $MRNABASE/trim.len > $MRNABASE/x
mv $MRNABASE/x $MRNABASE/trim.len
join $MRNABASE/raw.len $MRNABASE/trim.len  > $MRNABASE/both.len
# Make mRNA lift file from both.len
awk '{print 0,$1,$3,$1,$2}' $MRNABASE/both.len > $MRNABASE/mrna.lft
cp -p $MRNABASE/mrna.lft $TMPMRNA
cp -p $MRNABASE/S1.len $TMPMRNA
cp -p $MRNABASE/S2.len $TMPMRNA
cp -p $MRNABASE/raw.len $TMPMRNA

# Create temp work area for cluster run
mkdir -p $TMPMRNA/split
mkdir -p $TMPMRNA/lastz
mkdir -p $TMPMRNA/run.0
#rm -rf $TMPMRNA/run.0
# mkdir -p $TMPMRNA/run.0
# Split trim.fa into smaller FASTA files of about 1000000 bytes each by record.
faSplit about $MRNABASE/trim.fa 1000000 $TMPMRNA/split/mrna

#cd $TMPMRNA/split
# Creates a list of the mrn*.fa files after the split. 
ls $TMPMRNA/split/mrna*.fa > $TMPMRNA/S2.lst

# Create a doChain script that runs axtChain and pipes output to chainFilter
# then to chainToPsl.
# NOTE: this script is actually used in a later step so it can probably be 
# moved to where it is used.  
echo "#!/bin/bash" > $TMPMRNA/doChain
echo "BASE=$TMPMRNA" >> $TMPMRNA/doChain
echo "axtChain -linearGap=loose -verbose=0 -psl \$BASE/pslFilter/\$1.psl $TWOBIT -faQ \$BASE/trim.fa stdout | chainFilter -minScore=4000 stdin | chainToPsl stdin \$BASE/S1.len \$BASE/S2.len $TWOBIT \$BASE/trim.fa \$BASE/psl/\$1.psl" >> $TMPMRNA/doChain
chmod +x $TMPMRNA/doChain
# Make script to create chromosome directories in the lastz directory.
awk '{print "mkdir -p $TMPMRNA/lastz/"$1}' $TMPMRNA/S1.len > $TMPMRNA/create.dirs
source $TMPMRNA/create.dirs
# Create a list of chromosomes
awk '{print $1}' $TMPMRNA/S1.len > $TMPMRNA/S1.lst
# cd $TMPMRNA/run.0

# Make an axt directory. 
mkdir -p $TMPMRNA/lastz/axt
#Cluster job to run lastz to align mRNAs to genome
echo "#LOOP" > $TMPMRNA/run.0/template
echo "$SCRIPT/lastz.sh $TWOBIT/\$(path1) \$(path2) 10 62 {check out line $TMPMRNA/lastz/\$(root1)/\$(root2).psl} $TMPMRNA/lastz/axt/\$(root1) $TMPMRNA/S1.len $TMPMRNA/S2.len $LASTZPROG" >> $TMPMRNA/run.0/template
echo "#ENDLOOP" >> $TMPMRNA/run.0/template

# Copy over list of chromosomes and list of mRNA FASTA files
cp $TMPMRNA/S1.lst $TMPMRNA/run.0
cp $TMPMRNA/S2.lst $TMPMRNA/run.0
# Substitute contents of these files into template to create list of jobs
gensub2 $TMPMRNA/run.0/S1.lst $TMPMRNA/run.0/S2.lst $TMPMRNA/run.0/template $TMPMRNA/run.0/jobList

# Submit jobs to the cluster. 
ssh -T $CLUSTER "cd $TMPMRNA/run.0 ; /parasol/bin/para make jobList -ram=4g"
#    para create jobList
#    para try
#    para check
#    para push
#    para check
#Checking finished jobs

echo "run $SCRIPT/ucscRetroStep2.sh DEF after cluster job is finished"
#cd $OUTDIR
