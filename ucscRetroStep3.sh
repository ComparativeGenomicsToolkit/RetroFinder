#!/bin/bash 
# this script runs the retroFinder program (pslPseudo) that looks at mRNA alignments and creates a score file.
# a cluster job is used to run the program after extracting gene predictions, repeats and other annotation from the database.
#
# est alignments are extracted at the end for use in expression analysis in further stages.
#
set -beEu -o pipefail
source $1
echo "starting ucscRetroStep3.sh $1 for $DB"
mkdir -p $OUTDIR
cd $OUTDIR
mkdir -p result
mkdir -p result/axt
mkdir -p log
mkdir -p out

ls $NIB/ > $OUTDIR/S1.lst
cp $GENOME/$DB/chrom.sizes .

hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from all_mrna" > all_mrna.psl 
hgsql $DB -N -B -e "select matches,misMatches,repMatches,nCount,qNumInsert,qBaseInsert,tNumInsert,tBaseInsert,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts from refSeqAli" >> all_mrna.psl
rm -f all_mrna.psl.gz
gzip all_mrna.psl
#cat $RMSK/*.out |awk '{OFS="\t";print $5,$6,$7}' | grep -v position|grep -v sequence | tawk 'length($0)>2{print $0}' > rmsk.bed
rm -f rmsk.bed
if [ $RMSK == "rmsk" ]; then hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from rmsk" >> rmsk.bed ;
else 
for i in `cut -f 1 chrom.sizes` ;do hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from ${i}_rmsk" >> rmsk.bed ; done  ; fi
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET1 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET1.txt
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET2 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET2.txt 
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN from $NET3 where type <> 'gap' or (type = 'gap' and qN*100/(qEnd-qStart) > 75)" > \
$NET3.txt 
hgsql $DB -N -B -e "select chrom, chromStart, chromEnd from simpleRepeat" > simpleRepeat.bed 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 order by chrom, txStart, txEnd" | sort -k2,2 -k4,4n > $GENE1.tab
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2" > $GENE2.tab 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3" > $GENE3.tab 
hgsql $DB -N -B -e "select acc, version, name, type from refSeqAli a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" > cds.tab
hgsql $DB -N -B -e "select acc, version, name, type from all_mrna a , gbCdnaInfo g , cds c where qName = acc and cds = c.id" >> cds.tab
rm -f *.tab.gz *.txt.gz *.bed.gz
gzip *.tab &
gzip *.txt &
gzip *.bed 

#filter mrna and est with pslCDnaGenomeMatch
mkdir -p $TMPEST/mrna

rm -f mrna.2bit
ln $MRNABASE/mrna.2bit . -s

####
# run retro pipeine on cluster
####
mkdir -p $RETRODIR
cd $OUTDIR
mkdir -p run.0
cd run.0
wc -l ../split/*.psl |grep -v total | grep -v acc.lst| sort -nr | awk '{print $2}' |sed -e 's,../split/tmp,,'|sed -e 's/.psl//' > list
cp $OUTDIR/$1 .
echo "#LOOP" > gsub
echo "$SCRIPT/doBuildpk.sh \$(path1) $OUTDIR/$1 {check out exists $RESULT/pseudoGeneLink\$(path1).bed} " >> gsub
echo "#ENDLOOP" >> gsub

gensub2 list single gsub jobList
echo "Job Count"
wc -l spec
ssh -T $CLUSTER "cd $OUTDIR/run.0 ; para make jobList"
echo "check parasol status and then run ucscRetroStep4.sh DEF"
cd ..

