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

ls $NIB/*.nib > $OUTDIR/S1.lst
cp $GENOME/$DB/chrom.sizes .

#cat $RMSK/*.out |awk '{OFS="\t";print $5,$6,$7}' | grep -v position|grep -v sequence | tawk 'length($0)>2{print $0}' > rmsk.bed
if [[ -s rmsk.bed.gz ]] ; then
    echo "rmsk.bed.gz not refreshed"
else
    rm -f rmsk.bed rmsk.bed.gz
    if [ $RMSK == "rmsk" ]; then hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from rmsk" >> rmsk.bed ;
    else 
        for i in `cut -f 1 chrom.sizes` ;do hgsql $DB -N -B -e "select genoName , genoStart , genoEnd from ${i}_rmsk" >> rmsk.bed ; done  ; fi
fi
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET1 where tName not like '%hap%'" > \
$NET1.txt
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET2 where tName not like '%hap%'" > \
$NET2.txt 
hgsql $DB -N -B -e "select tName, tStart, tEnd, level, qName, qStart, qEnd, type, qN, tR+tTrf from $NET3 where tName not like '%hap%'" > \
$NET3.txt 
hgsql $DB -N -B -e "select chrom, chromStart, chromEnd from simpleRepeat" > simpleRepeat.bed 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 order by chrom, txStart, txEnd" | sort -k2,2 -k4,4n > $GENE1.tab
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2" > $GENE2.tab 
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3" > $GENE3.tab 
rm -f $GENE1.tab.gz $GENE2.tab.gz $GENE3.tab.gz *.txt.gz *.bed.gz
gzip *.tab &
gzip *.txt &
gzip *.bed 

rm -f mrna.2bit
ln $MRNABASE/mrna.2bit . -s

####
# run retro pipeline on cluster
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
wc -l jobList
ssh -T $CLUSTER "cd $OUTDIR/run.0 ; /parasol/bin/para make jobList"
echo "check parasol status and then run ucscRetroStep4.sh DEF"
cd ..
