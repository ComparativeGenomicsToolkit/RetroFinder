#align human mrnas to hg18 using lastz
DB=hg18
GBDB=hg

mkdir -p /hive/data/genomes/$DB/bed/mrnaBlastz/
cd /hive/data/genomes/$DB/bed/mrnaBlastz/

#pull latest mrna and refseq from genbank 
/cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
   genbank mrna stdout | tr acgt ACGT > mrna.fa 
/cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
   refSeq mrna stdout | tr acgt ACGT > refseq.fa 
cat mrna.fa refseq.fa > raw.fa

#remmove polyA tail before aligning
faTrimPolyA raw.fa trim.fa 
mkdir -p /san/sanvol1/scratch/mrnaBlastz/$DB
cp -p trim.fa /san/sanvol1/scratch/mrnaBlastz/$DB
faToTwoBit raw.fa mrna.2bit
twoBitToFa mrna.2bit stdout |faToTwoBit stdin -stripVersion mrnaNoversion.2bit
faSize trim.fa -detailed > trim.len
cp /hive/data/genomes/$DB/chrom.sizes S1.len
cp trim.len S2.len
faSize raw.fa -detailed > raw.len
sort raw.len > x
mv x raw.len
sort trim.len > x
mv x trim.len
join raw.len trim.len  > both.len
awk '{print 0,$1,$3,$1,$2}' both.len >  mrna.lft
cp -p mrna.lft /san/sanvol1/scratch/mrnaBlastz/$DB
cp -p S1.len /san/sanvol1/scratch/mrnaBlastz/$DB
cp -p S2.len /san/sanvol1/scratch/mrnaBlastz/$DB
cp -p raw.len /san/sanvol1/scratch/mrnaBlastz/$DB


# create temp work area for cluster run
mkdir -p /san/sanvol1/scratch/mrnaBlastz/$DB/split
mkdir -p /san/sanvol1/scratch/mrnaBlastz/$DB/lastz/
mkdir -p /san/sanvol1/scratch/mrnaBlastz/$DB/run.0
faSplit about trim.fa 1000000 /san/sanvol1/scratch/mrnaBlastz/$DB/split/mrna

cd /san/sanvol1/scratch/mrnaBlastz/$DB/split
ls mrna*.fa |awk '{print "/san/sanvol1/scratch/mrnaBlastz/'$DB'/split/"$1}' > /san/sanvol1/scratch/mrnaBlastz/$DB/S2.lst

cd ..
echo "#\!/bin/bash -fe" > doChain
echo "BASE=/san/sanvol1/scratch/mrnaBlastz/$DB" >> doChain
echo "axtChain -linearGap=loose -psl \$BASE/pslFilter/\$1.psl /scratch/data/$DB/$DB.2bit -faQ \$BASE/trim.fa stdout | chainFilter -minScore=4000 stdin | chainToPsl stdin S1.len S2.len nib.lst trim.fa psl/\$1.psl" >> doChain
chmod +x doChain
awk '{print "mkdir -p /san/sanvol1/scratch/mrnaBlastz/$DB/lastz/"$1}' S1.len |grep -v random > create.dirs
source create.dirs
grep -v "random" /scratch/data/$DB/chrom.sizes  |cut -f 1 > S1.lst
cd /san/sanvol1/scratch/mrnaBlastz/$DB/run.0

#cluster job to run lastz to align mRNAs to genome
echo "#LOOP" > template
echo "/cluster/home/baertsch/baertsch/scripts/lastz.sh /scratch/data/$DB/$DB.2bit/\$(path1) \$(path2) 10 62 {check out line+ /san/sanvol1/scratch/mrnaBlastz/$DB/lastz/\$(root1)/\$(root2).psl} /san/sanvol1/scratch/mrnaBlastz/$DB/S1.len /san/sanvol1/scratch/mrnaBlastz/$DB/S2.len" >> template
echo "#ENDLOOP" >> template

cp ../S1.lst .
cp ../S2.lst .
gensub2 S1.lst S2.lst template jobList
ssh -T pk "cd /san/sanvol1/scratch/mrnaBlastz/$DB/run.0 ; para make jobList"
#    para create jobList
#    para try
#    para check
#    para push
#    para check
#Checking finished jobs

#run ~/baertsch/scripts/ucscRetroStep2.sh next
