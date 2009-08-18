#!/bin/bash 
set -beEu -o pipefail
if [ $BASH_ARGC ]; then
       echo "Parameter file is $1"
else
       echo "usage $0 DEF"
       exit 3
fi
source $1
echo running mrna alignment against $DB
#align mrnas to $DB using lastz
mkdir -p $TMPMRNA
mkdir -p $MRNABASE
cd $MRNABASE

#pull latest mrna and refseq from genbank 
if [[ -s mrna.fa ]] ; then
    echo "mrna.fa exists, extraction from genbank skipped"
else
    echo "extracting mRNA and reqseq from Genbank "
    echo "/cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
       genbank mrna stdout | tr acgt ACGT "
    /cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
       genbank mrna stdout | tr acgt ACGT > tmp.fa ; mv tmp.fa mrna.fa 
fi
if [[ -s refseq.fa ]] ; then
    echo "refseq.fa exists, extraction from genbank skipped"
else
    /cluster/data/genbank/bin/x86_64/gbGetSeqs -db=$GBDB -inclVersion -native -gbRoot=/cluster/data/genbank \
       refSeq mrna stdout | tr acgt ACGT > refseq.fa 
fi
echo "cat mrna.fa refseq.fa output-to raw.fa"
cat mrna.fa refseq.fa > raw.fa

#remove polyA tail before aligning
echo "faTrimPolyA raw.fa trim.fa "
faTrimPolyA raw.fa trim.fa 
mkdir -p $TMPMRNA 
cp -p trim.fa $TMPMRNA 
faToTwoBit raw.fa mrna.2bit
twoBitToFa mrna.2bit stdout |faToTwoBit stdin -stripVersion mrnaNoversion.2bit
faSize trim.fa -detailed > trim.len
grep -v chrM $LOCAL/chrom.sizes > S1.len
cp trim.len S2.len
faSize raw.fa -detailed > raw.len
sort raw.len > x
mv x raw.len
sort trim.len > x
mv x trim.len
join raw.len trim.len  > both.len
awk '{print 0,$1,$3,$1,$2}' both.len >  mrna.lft
cp -p mrna.lft $TMPMRNA
cp -p S1.len $TMPMRNA
cp -p S2.len $TMPMRNA
cp -p raw.len $TMPMRNA


# create temp work area for cluster run
mkdir -p $TMPMRNA/split
mkdir -p $TMPMRNA/lastz
mkdir -p $TMPMRNA/run.0
faSplit about trim.fa 1000000 $TMPMRNA/split/mrna

cd $TMPMRNA/split
ls mrna*.fa |awk '{print "'$TMPMRNA'/split/"$1}' > $TMPMRNA/S2.lst

cd ..
echo "#!/bin/bash" > doChain
echo "BASE=$TMPMRNA" >> doChain
echo "axtChain -linearGap=loose -psl \$BASE/pslFilter/\$1.psl $LOCAL/$DB.2bit -faQ \$BASE/trim.fa stdout | chainFilter -minScore=4000 stdin | chainToPsl stdin S1.len S2.len nib.lst trim.fa psl/\$1.psl" >> doChain
chmod +x doChain
awk '{print "mkdir -p $TMPMRNA/lastz/"$1}' S1.len |grep -v random > create.dirs
source create.dirs
grep -v "random" $LOCAL/chrom.sizes |grep -v chrM |cut -f 1 > S1.lst
cd $TMPMRNA/run.0

#cluster job to run lastz to align mRNAs to genome
echo "#LOOP" > template
echo "$SCRIPT/lastz.sh $LOCAL/$DB.2bit/\$(path1) \$(path2) 10 62 {check out line+ $TMPMRNA/lastz/\$(root1)/\$(root2).psl} $TMPMRNA/S1.len $TMPMRNA/S2.len" >> template
echo "#ENDLOOP" >> template

cp ../S1.lst .
cp ../S2.lst .
gensub2 S1.lst S2.lst template jobList

ssh -T $CLUSTER "cd $TMPMRNA/run.0 ; /parasol/bin/para make jobList"
#    para create jobList
#    para try
#    para check
#    para push
#    para check
#Checking finished jobs

cd $TMPMRNA
echo run $SCRIPT/ucscRetroStep2.sh DEF after cluster job is finished
$SCRIPT/ucscRetroStep2.sh DEF 
