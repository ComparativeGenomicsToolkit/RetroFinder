#!/bin/bash
set -bevEu -o pipefail
source $1
echo '-------- script analyseExpress.sh -------------------'
####################
# analyze expression 
####################
#make exon list
# cd $OUTDIR
# Checks that output of filterMrna.sh and filterEst.sh has been produced 
# otherwise prints error. 
if [ ! -s $OUTDIR/all_mrnaFiltered.psl.gz ] 
then
   echo "$OUTDIR/all_mrnaFiltered.psl.gz not found or empty, check for errors when running filterMrna.sh"
   exit 3
fi
 
if [ ! -s $OUTDIR/estFiltered.psl.gz ] 
then
   echo "estFiltered.psl.gz not found or empty, check for errors when running script filterEst.sh"
   exit 3
fi
# Create exp directory
mkdir -p $OUTDIR/$EXPDIR
echo -n "Working directory: "
echo $OUTDIR/$EXPDIR
# Get genePred file for transcripts with more than one exon for GENE1 
# (typically knownGene)
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE1 where exonCount > 1" > $OUTDIR/$EXPDIR/$GENE1.multiExon.genePred
# Get genePred file for transcripts with more than one exon for GENE2 
# (typically refGene)
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE2 where exonCount > 1" > $OUTDIR/$EXPDIR/$GENE2.multiExon.genePred
# Get genePred file for transcripts with more than one exon for GENE2 
# (typically ensGene)
hgsql $DB -N -B -e "select name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from $GENE3 where exonCount > 1" > $OUTDIR/$EXPDIR/$GENE3.multiExon.genePred
echo "genePredFilter $OUTDIR/$EXPDIR/$GENE1.multiExon.genePred $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred -cdsExons=2"
# Filter genePred files and require at least 2 CDS exons
${BINDIR}/genePredFilter $OUTDIR/$EXPDIR/$GENE1.multiExon.genePred $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred -cdsExons=2
echo "genePredFilter $OUTDIR/$EXPDIR/$GENE2.multiExon.genePred $OUTDIR/$EXPDIR/$GENE2.multiCDSExon.genePred -cdsExons=2"
${BINDIR}/genePredFilter $OUTDIR/$EXPDIR/$GENE2.multiExon.genePred $OUTDIR/$EXPDIR/$GENE2.multiCDSExon.genePred -cdsExons=2
echo "genePredFilter $OUTDIR/$EXPDIR/$GENE3.multiExon.genePred $OUTDIR/$EXPDIR/$GENE3.multiCDSExon.genePred -cdsExons=2"
${BINDIR}/genePredFilter $OUTDIR/$EXPDIR/$GENE3.multiExon.genePred $OUTDIR/$EXPDIR/$GENE3.multiCDSExon.genePred -cdsExons=2
# Get counts of lines in these output files
wc -l $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE2.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred
echo "cat $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE2.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE3.multiCDSExon.genePred to $OUTDIR/$EXPDIR/all.multiCds.gp"
# Cat together the files of genePreds with >1 exon and >=2 CDS exons. 
cat $OUTDIR/$EXPDIR/$GENE1.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE2.multiCDSExon.genePred $OUTDIR/$EXPDIR/$GENE3.multiCDSExon.genePred > $OUTDIR/$EXPDIR/all.multiCds.gp

# Load these files of genePred with >1 exon and >=2 CDS exons for GENE1 and 
# GENE2 into database tables.
ldHgGene $DB rb${GENE1}Multi $OUTDIR/$EXPDIR/${GENE1}.multiCDSExon.genePred -predTab;
ldHgGene $DB rb${GENE2}Multi $OUTDIR/$EXPDIR/${GENE2}.multiCDSExon.genePred -predTab;

# Extract CDS regions from genePred tables for GENE1 and GENE2
featureBits $DB ${GENE1}:cds -bed=$OUTDIR/$EXPDIR/$GENE1.cds.bed;
featureBits $DB ${GENE2}:cds -bed=$OUTDIR/$EXPDIR/$GENE2.cds.bed;
# Intersect the CDS regions with the multiCds tables using featureBits
featureBits $DB rb${GENE1}Multi $OUTDIR/$EXPDIR/$GENE1.cds.bed -bed=$OUTDIR/$EXPDIR/$GENE1.multiCds.bed ;
featureBits $DB rb${GENE2}Multi $OUTDIR/$EXPDIR/$GENE2.cds.bed -bed=$OUTDIR/$EXPDIR/$GENE2.multiCds.bed ;

echo "overlap with ests - use cluster jobs to speed this step"
# Make an estSplit directory
mkdir -p $OUTDIR/estSplit
# Default max target count is 300 as open a file handle for each target.
# The human GRCh38 (hg38) assembly has a number of alt alleles and there 
# are 455 targets so use option to increase the number allowed. 
pslSplitOnTarget -maxTargetCount=500 $OUTDIR/estFiltered.psl.gz estSplit
# Remove previous jobList
rm -f $OUTDIR/estSplit/jobList
# Create template for estStat.sh jobs
echo "#LOOP"> $OUTDIR/estSplit/template
echo "$SCRIPT/estStat.sh \$(root1) {check out line pseudoEst.\$(root1).bed}">>$OUTDIR/estSplit/template
echo "#ENDLOOP">> $OUTDIR/estSplit/template
#awk '{print $1}' $OUTDIR/S1.len |grep -v chrM | gensub2 stdin single $OUTDIR/estSplit/template $OUTDIR/estSplit/jobList
# Get a list of the chrom names from the PSLs in the estSplit directory, 
# remove chrM and substitute the names in the template to create a jobList
ls $OUTDIR/estSplit/*psl | sed -e "s,$OUTDIR/estSplit/,," | sed -e 's/.psl//' | grep -v chrM | gensub2 stdin single $OUTDIR/estSplit/template $OUTDIR/estSplit/jobList
# Run the jobs in jobList on the cluster
ssh -T $CLUSTER "cd $OUTDIR/estSplit ; /parasol/bin/para make jobList"

# Cat together the output files from estStat.sh into a single file
# and count the lines in this file 
cat $OUTDIR/estSplit/pseudoEst*.bed > $OUTDIR/$EXPDIR/pseudoEstAll.bed
wc -l $OUTDIR/$EXPDIR/pseudoEstAll.bed

# mrnaToGene converts PSL alignments with CDS annotation from genbank to 
# gene annotations in genePred format
# all_mrnaFiltered.psl is converted to genePred using CDS from the $DB
# and merging gaps in CDS no larger than 10 nts and keep sequences with 
# invalid CDS. 
echo "mrnaToGene $OUTDIR/all_mrnaFiltered.psl.gz -cdsMergeSize=10 all_mrna.gp -cdsDb=$DB -keepInvalid "
mrnaToGene $OUTDIR/all_mrnaFiltered.psl.gz -cdsMergeSize=10 $OUTDIR/$EXPDIR/all_mrna.gp -cdsDb=$DB -keepInvalid > $OUTDIR/$EXPDIR/all_mrna.log 2> $OUTDIR/$EXPDIR/all_mrna.err
# then do the same but merge gaps in the UTR no larger than 10 nts
mrnaToGene $OUTDIR/all_mrnaFiltered.psl.gz -cdsMergeSize=10 -utrMergeSize=10 $OUTDIR/$EXPDIR/all_mrna_utr.gp -cdsDb=$DB -keepInvalid > $OUTDIR/$EXPDIR/all_mrna_utr.log 2> $OUTDIR/$EXPDIR/all_mrna_utr.err

# From these genePreds, select those that are multi exonic so exons >1 
awk '$8>1{print }' $OUTDIR/$EXPDIR/all_mrna.gp > $OUTDIR/$EXPDIR/all_mrna_multiExon.gp
awk '$8>1{print }' $OUTDIR/$EXPDIR/all_mrna_utr.gp > $OUTDIR/$EXPDIR/all_mrna_multiExonUTR.gp

#expression according to knownGene
# Get the entries from hgncXref which is gene symbol, refSeq id, 
# UniProt id, HGNC id, Entrez and description.
hgsql $PDB -N -B -e "select * from hgncXref " > $OUTDIR/$EXPDIR/hgncXref.tab
# if $GENE1 is knownGene:
if [[ $GENE1 == "knownGene" ]] 
then
# join knownGene and kgTxInfo tables and get output where category is coding
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'coding'" > $OUTDIR/$EXPDIR/kgCoding.gp
# join knownGene and kgTxInfo tables and get output where category is noncoding
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'noncoding'" > $OUTDIR/$EXPDIR/kgNoncoding.gp
# join knownGene and kgTxInfo tables and get output where category is nearcoding
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'nearcoding'" > $OUTDIR/$EXPDIR/kgNearcoding.gp
# join knownGene and kgTxInfo tables and get output where category is antisense
hgsql $DB -N -B -e "select kg.* from knownGene kg, kgTxInfo i where kg.name = i.name and category = 'antisense'" > $OUTDIR/$EXPDIR/kgAntisense.gp
# Get overlap where each subset of known genes is the select file and 
# ucscRetroInfo.tab is the inFile
for i in Coding Noncoding Nearcoding Antisense ; do overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/kg${i}.gp $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/retroKg${i}.bed; done
# unzip knownGene file and sort
zcat $OUTDIR/knownGene.tab.gz |sort > $OUTDIR/knownGene.sort.gp
# select knownGene ids that represent genes in UniProt with evidence type <=1
hgsql $DB -N -B -e "select distinct kgID from kgSpAlias, uniProt.proteinEvidence where spId = acc and proteinEvidenceType <=1" |sort> $OUTDIR/$EXPDIR/kgProtEvidence.id
# select knownGene ids that represent genes in UniProt with evidence type <=2
hgsql $DB -N -B -e "select distinct kgID from kgSpAlias, uniProt.proteinEvidence where spId = acc and proteinEvidenceType <=2" |sort> $OUTDIR/$EXPDIR/kgTransEvidence.id
# select rows from sorted knownGene file by matching ids from column 1 with 
# the ids in column 1 of the protein evidence files. 
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/kgProtEvidence.id 1 $OUTDIR/knownGene.sort.gp > $OUTDIR/$EXPDIR/knownGeneProt.gp
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/kgTransEvidence.id 1 $OUTDIR/knownGene.sort.gp > $OUTDIR/$EXPDIR/knownGeneTrans.gp

# use the set of known genes with protein evidence to intersect with 
# GENE1 (knownGene) CDS file
overlapSelect $OUTDIR/$EXPDIR/knownGeneProt.gp $OUTDIR/$EXPDIR/${GENE1}.cds.bed $OUTDIR/$EXPDIR/knownGeneProtCds.bed -selectFmt=genePred
# use the set of known genes with protein and transcript evidence to 
# intersect with GENE1 (knownGene) CDS file
overlapSelect $OUTDIR/$EXPDIR/knownGeneTrans.gp $OUTDIR/$EXPDIR/${GENE1}.cds.bed $OUTDIR/$EXPDIR/knownGeneTransCds.bed -selectFmt=genePred
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/knownGeneProtCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/ucscRetroProtein.bed
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/knownGeneTransCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/ucscRetroTranscript.bed
fi

echo "Exon Shuffling"
# exon shuffling
##cat $OUTDIR/$TABLE.bed $OUTDIR/retroMrnaInfoZnf.bed > ucscRetroInfoZnf.bed
# select records from ucscRetroInfo.bed based on overlap with those in the file
# GENE1 annotations thah have >1 exon and >=2 CDS exons
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE1}.multiCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE1}.multiCds.bed $OUTDIR/$TABLE.bed -statsOutput $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.out
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE1}.multiCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/pseudo${GENE1}Cds50.bed -overlapThreshold=0.50
# select records from ucscRetroInfo.bed based on overlap with those in the file
# GENE2 annotations thah have >1 exon and >=2 CDS exons
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE2}.multiCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE2}.multiCds.bed $OUTDIR/$TABLE.bed -statsOutput $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.out
# select records from ucscRetroInfo.bed based on overlap with those in the file
# GENE2 annotations thah have >1 exon and >=2 CDS exons and 
# overlapThreshold = 0.5 so each infile record must be overlapped at least
# 50% by the select record. 
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/${GENE2}.multiCds.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/pseudo${GENE2}Cds50.bed -overlapThreshold=0.50

# Select those retro records that overlap GENE2 annotations with >1 exon and
# >=2 CDS exons and overlap $GENE2 annotations
overlapSelect -inCoordCols=0,1,2,5,3 -selectFmt=genePred $OUTDIR/${GENE2}.tab.gz $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed $OUTDIR/$EXPDIR/shuffleEns.bed
# select the retro records that overlap GENE2 annotations with >1 exon and
# >=2 CDS exons and overlap the GENE2 genePreds with >1 exon and >=2 CDS exons
overlapSelect -inCoordCols=0,1,2,5,3 -selectFmt=genePred $OUTDIR/$EXPDIR/${GENE2}.multiCDSExon.genePred $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed $OUTDIR/$EXPDIR/shuffleEnsMulti.bed

# select retros overlapping ESTs that do not overlap GENE1 annotations with
# >1 exon and >=2 exon
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed $OUTDIR/$EXPDIR/pseudoEstAll.bed $OUTDIR/$EXPDIR/pseudoEstAllNotShuffle.bed -nonOverlapping

if [[ $GENE1 == "knownGene" ]] 
then
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/knownGeneProt.gp $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed $OUTDIR/$EXPDIR/pseudo${GENE1}ProtCds.bed
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/knownGeneProt.gp $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed $OUTDIR/$EXPDIR/pseudo${GENE2}ProtCds.bed
fi

ls -l $OUTDIR/$EXPDIR/pseudo${GENE1}*.bed
wc -l $OUTDIR/$EXPDIR/pseudo${GENE1}*.bed
cut -f 4 $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed > $OUTDIR/$EXPDIR/shuffle.list
awk '{print "update '$TABLE' set type = \"expressed shuffle\" where name = \""$1"\" ;"}' $OUTDIR/$EXPDIR/shuffle.list> $OUTDIR/$EXPDIR/updateShuffle.sql
echo "Updating mysql database with updateShuffle.sql `wc -l $OUTDIR/$EXPDIR/updateShuffle.sql` records."
hgsql $DB < $OUTDIR/$EXPDIR/updateShuffle.sql

# expression analysis based on mrna and est overlap
# grab results of cluster job
cat $OUTDIR/estSplit/est.*.id > $OUTDIR/$EXPDIR/est.id
cat $OUTDIR/estSplit/stat.*.out > $OUTDIR/$EXPDIR/stat.out
cat $OUTDIR/estSplit/statagg.*.out |sort |grep -v "#"> $OUTDIR/$EXPDIR/statagg.out
awk '{print $1}' $OUTDIR/$EXPDIR/stat.out |uniq -c |awk '{print $2,$1}' |sort> $OUTDIR/$EXPDIR/estCount.out
join $OUTDIR/$EXPDIR/estCount.out $OUTDIR/$EXPDIR/statagg.out > $OUTDIR/$EXPDIR/estCoverage.out
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/all_mrnaFiltered.psl.gz $OUTDIR/$TABLE.bed -statsOutput stdout |sort > $OUTDIR/$EXPDIR/mrna.out
zcat $OUTDIR/splicedEst.psl.gz |cut -f 10 |sort |uniq> $OUTDIR/$EXPDIR/splicedEst.id
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/splicedEst.psl.gz $OUTDIR/$TABLE.bed -statsOutput stdout | sort > $OUTDIR/$EXPDIR/splicedstat.out
awk '{print $1}' $OUTDIR/$EXPDIR/splicedstat.out |uniq -c |awk '{print $2,$1}' |sort> $OUTDIR/$EXPDIR/splicedEstCount.out
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/splicedEst.psl.gz $OUTDIR/$TABLE.bed -statsOutput -aggregate stdout | sort > $OUTDIR/$EXPDIR/splicedagg.out
join $OUTDIR/$EXPDIR/splicedEstCount.out $OUTDIR/$EXPDIR/splicedagg.out > $OUTDIR/$EXPDIR/splicedEstCoverage.out

awk '$3>0.50{print $1}' $OUTDIR/$EXPDIR/mrna.out |sort > $OUTDIR/$EXPDIR/mrna.id
awk '$2>=10 && $3>0.50{print $0}' $OUTDIR/$EXPDIR/estCoverage.out > $OUTDIR/$EXPDIR/est10.out
awk '$2>=5 && $3>0.50{print $0}' $OUTDIR/$EXPDIR/estCoverage.out > $OUTDIR/$EXPDIR/est5.out
awk '$2>=1 && $3>0.50{print $0}' $OUTDIR/$EXPDIR/estCoverage.out > $OUTDIR/$EXPDIR/est.out
awk '$2>=5 && $3>0.50{print $0}' $OUTDIR/$EXPDIR/splicedEstCoverage.out > $OUTDIR/$EXPDIR/splicedEst5.out
awk '$2>=10 && $3>0.50{print $0}' $OUTDIR/$EXPDIR/splicedEstCoverage.out > $OUTDIR/$EXPDIR/splicedEst10.out
awk '{print $1}' $OUTDIR/$EXPDIR/est.out |sort > $OUTDIR/$EXPDIR/est1.id
awk '{print $1}' $OUTDIR/$EXPDIR/est5.out |sort > $OUTDIR/$EXPDIR/est5.id
awk '{print $1}' $OUTDIR/$EXPDIR/est10.out |sort > $OUTDIR/$EXPDIR/est10.id
awk '{print $1}' $OUTDIR/$EXPDIR/splicedEst5.out |sort > $OUTDIR/$EXPDIR/splicedEst5.id
awk '{print $1}' $OUTDIR/$EXPDIR/splicedEst10.out |sort > $OUTDIR/$EXPDIR/splicedEst10.id

join $OUTDIR/$EXPDIR/est10.id $OUTDIR/$EXPDIR/mrna.id > $OUTDIR/$EXPDIR/est10Mrna.id
join $OUTDIR/$EXPDIR/est5.id $OUTDIR/$EXPDIR/mrna.id > $OUTDIR/$EXPDIR/est5Mrna.id
grep -F -f $OUTDIR/$EXPDIR/est5.id $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst5.bed
cat $OUTDIR/$EXPDIR/mrna.id $OUTDIR/$EXPDIR/est10.id |sort |uniq > $OUTDIR/$EXPDIR/mrnaEst10.id
cat $OUTDIR/$EXPDIR/mrna.id $OUTDIR/$EXPDIR/est5.id |sort |uniq > $OUTDIR/$EXPDIR/mrnaEst5.id
cat $OUTDIR/$EXPDIR/mrna.id $OUTDIR/$EXPDIR/est1.id |sort |uniq > $OUTDIR/$EXPDIR/mrnaEst1.id
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/est5Mrna.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/est10Mrna.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst10AndMrna.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/mrnaEst1.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEstMrna.filter.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/mrnaEst5.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst5Mrna.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/mrnaEst10.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst10Mrna.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/est5.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst5.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/est10.id 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst10.bed
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/pseudoEst10Mrna.bed $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed $OUTDIR/$EXPDIR/est10MrnaRefSeq.bed

overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/retroKgCoding.bed $OUTDIR/$EXPDIR/pseudoEst5Mrna.bed $OUTDIR/$EXPDIR/pseudoEst5MrnaNotKg.bed -nonOverlapping
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/retroKgCoding.bed $OUTDIR/$EXPDIR/pseudoEst10Mrna.bed $OUTDIR/$EXPDIR/pseudoEst10MrnaNotKg.bed -nonOverlapping
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/retroKgCoding.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrnaNotKg.bed -nonOverlapping

#sort -k2,2nr estCoverage.out |head -50|awk '{print $1}' > top50
#echo $SCRIPT/selectById -tsv 1 top50 4 $OUTDIR/$TABLE.bed   top50.bed
#$SCRIPT/selectById -tsv 1 top50 4 $OUTDIR/$TABLE.bed  > top50.bed

echo "expressed retro stats"
wc -l $OUTDIR/$EXPDIR/*.id |sort -n
overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/${GENE2}.tab.gz $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/pseudo${GENE2}.bed -selectFmt=genePred

#regenerate gene predictions from expressed retros

#set with 1 or more EST (tight) not shuffle events
#look for orf in window 200 bp around retro
tawk '$2>200{$2=$2-200; $3=$3+200;print $0}$2<=200{$3=$3+200;print $0}' $OUTDIR/$EXPDIR/pseudoEstAllNotShuffle.bed > $OUTDIR/$EXPDIR/pseudoEstAll.bed
echo "orfBatch $DB pseudoEstAll.bed pseudoEstAll.out pseudoEstAll.gp to borfEst.out"
${BINDIR}/orfBatch $DB $OUTDIR/$EXPDIR/pseudoEstAll.bed $OUTDIR/$EXPDIR/pseudoEstAll.out $OUTDIR/$EXPDIR/pseudoEstAll.gp > $OUTDIR/$EXPDIR/borfEst.out
echo gene-check -genome-seqs $NIB $OUTDIR/$EXPDIR/pseudoEstAll.gp $OUTDIR/$EXPDIR/checkEst.rdb
gene-check  -genome-seqs $NIB $OUTDIR/$EXPDIR/pseudoEstAll.gp $OUTDIR/$EXPDIR/checkEst.rdb
awk '$7=="ok"&& $25=="noStart" || $25==""{print $1}' $OUTDIR/$EXPDIR/checkEst.rdb > $OUTDIR/$EXPDIR/goodOrf.list

$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/goodOrf.list 4 $OUTDIR/$EXPDIR/pseudoEstAll.out > $OUTDIR/$EXPDIR/pseudoEstAll.good.out

echo "histogram of best.orf scores"
awk '{print $52}' $OUTDIR/$EXPDIR/pseudoEstAll.good.out |textHistogram stdin -maxBinCount=40 -binSize=10
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed weak\", posConf = \""$52"\" where name = \""$4"\" ;"}' $OUTDIR/$EXPDIR/pseudoEstAll.good.out > $OUTDIR/$EXPDIR/updateExp.sql
echo "hgsql $DB $OUTDIR/$EXPDIR/updateExp.sql"
echo "Updating mysql database with $OUTDIR/$EXPDIR/updateExp.sql `wc -l $OUTDIR/$EXPDIR/updateExp.sql` records."
hgsql $DB < $OUTDIR/$EXPDIR/updateExp.sql

#set with 5 Est and 1 spliced mRNA
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrnaNotShuffle.bed -nonOverlapping
tawk '$2>200{$2=$2-200; $3=$3+200;print $0}$2<=200{$3=$3+200;print $0}' $OUTDIR/$EXPDIR/pseudoEst5AndMrnaNotShuffle.bed > $OUTDIR/$EXPDIR/pseudoEst5AndMrna.window.bed
${BINDIR}/orfBatch $DB $OUTDIR/$EXPDIR/pseudoEst5AndMrna.window.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrna.out $OUTDIR/$EXPDIR/pseudoEst5AndMrna.gp > $OUTDIR/$EXPDIR/borfEst5AndMrna.out
echo gene-check -genome-seqs $NIB $OUTDIR/$EXPDIR/pseudoEst5AndMrna.gp $OUTDIR/$EXPDIR/checkEst5AndMrna.rdb
gene-check -genome-seqs $NIB $OUTDIR/$EXPDIR/pseudoEst5AndMrna.gp $OUTDIR/$EXPDIR/checkEst5AndMrna.rdb
tawk '$7=="ok" && $25=="noStart" || $25==""{print $1}' $OUTDIR/$EXPDIR/checkEst5AndMrna.rdb > $OUTDIR/$EXPDIR/goodOrf5AndMrna.list
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/goodOrf5AndMrna.list 4 $OUTDIR/$EXPDIR/pseudoEst5AndMrna.out > $OUTDIR/$EXPDIR/pseudoEst5AndMrna.good.out
awk '{print "update '$TABLE' set thickStart = "$7", thickEnd = "$8" , type = \"expressed strong\", posConf = \""$52"\" where name = \""$4"\" ;"}' $OUTDIR/$EXPDIR/pseudoEst5AndMrna.good.out > $OUTDIR/$EXPDIR/updateExp5.sql
echo "hgsql $DB $OUTDIR/$EXPDIR/updateExp5.sql"
echo "Updating mysql database with $OUTDIR/$EXPDIR/updateExp5.sql `wc -l $OUTDIR/$EXPDIR/updateExp5.sql` records."
hgsql $DB < $OUTDIR/$EXPDIR/updateExp5.sql

#bad orfs
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' $OUTDIR/$EXPDIR/checkEst.rdb > $OUTDIR/$EXPDIR/badOrf.list
awk '{print "update '$TABLE' set type = \"expressed weak noOrf\" where name = \""$1"\" ;"}' $OUTDIR/$EXPDIR/badOrf.list> $OUTDIR/$EXPDIR/updateNoOrf.sql
echo "Updating mysql database with $OUTDIR/$EXPDIR/updateNoOrf.sql `wc -l $OUTDIR/$EXPDIR/updateNoOrf.sql` records."
hgsql $DB < $OUTDIR/$EXPDIR/updateNoOrf.sql
awk '$6!="ok"|| $7!="ok" || ($7=="ok" && $25!="noStart"&& $25!=""){print $1}' $OUTDIR/$EXPDIR/checkEst5AndMrna.rdb > $OUTDIR/$EXPDIR/badOrf5.list
awk '{print "update '$TABLE' set type = \"expressed strong noOrf\" where name = \""$1"\" ;"}' $OUTDIR/$EXPDIR/badOrf5.list> $OUTDIR/$EXPDIR/updateStrongNoOrf.sql
echo "Updating mysql database with $OUTDIR/$EXPDIR/updateStrongNoOrf.sql `wc -l $OUTDIR/$EXPDIR/updateStrongNoOrf.sql` records."
hgsql $DB < $OUTDIR/$EXPDIR/updateStrongNoOrf.sql
wc -l $OUTDIR/$EXPDIR/update*.sql

$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/goodOrf.list 1 $OUTDIR/$EXPDIR/pseudoEstAll.gp > $OUTDIR/$EXPDIR/pseudoExpressed.gp
echo "ldHgGene $DB pseudoExpressed $OUTDIR/$EXPDIR/pseudoExpressed.gp -genePredExt -predTab"
ldHgGene $DB ucscRetroExpressed$VERSION $OUTDIR/$EXPDIR/pseudoExpressed.gp -genePredExt -predTab
genePredToBed $OUTDIR/$EXPDIR/pseudoExpressed.gp $OUTDIR/$EXPDIR/pseudoTmp.bed
$SCRIPT/selectById -tsv 1 $OUTDIR/$EXPDIR/goodOrf.list 4 $OUTDIR/$EXPDIR/pseudoEstAll.bed > $OUTDIR/$EXPDIR/pseudoExpressed.bed
echo "length histogram of coding region"
awk '{print $7-$6}' $OUTDIR/$EXPDIR/pseudoExpressed.gp|textHistogram stdin -maxBinCount=100 -binSize=50

tawk '($8-$7)>300{print }' $OUTDIR/$EXPDIR/pseudoTmp.bed> $OUTDIR/$EXPDIR/pseudoTmp2.bed
$SCRIPT/selectById -tsv 4 $OUTDIR/$EXPDIR/pseudoTmp2.bed 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudoEst100AA.bed
echo tawk '$8-$7>300{print }' $OUTDIR/$EXPDIR/pseudoExpressed.bed redirect $OUTDIR/$EXPDIR/pseudoEst100AA.bed
$SCRIPT/selectById -tsv 4 $OUTDIR/$EXPDIR/pseudoEst5AndMrna.good.out 4 $OUTDIR/$TABLE.bed > $OUTDIR/$EXPDIR/pseudo5Est100AA.bed
wc -l $OUTDIR/$EXPDIR/pseudoExpressed.bed $OUTDIR/$EXPDIR/pseudoEst100AA.bed $OUTDIR/$EXPDIR/goodOrf.list $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrna.good.out
overlapSelect -selectCoordCols=0,1,2,5,3 -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed $OUTDIR/$EXPDIR/pseudo5Est100AA.bed $OUTDIR/$EXPDIR/pseudoEst5AndMrna100AA.bed

#retros involved in alt-splicing

if [[ -n $ALTSPLICE ]] 
then
echo "hgsql $DB -N -B -e select * from $ALTSPLICE " 
hgsql $DB -N -B -e "select * from $ALTSPLICE " |cut -f2-19 | agxToBed stdin $OUTDIR/$EXPDIR/altSplice.bed
fi
if [ -f $OUTDIR/$EXPDIR/altSplice.bed ] 
then
 echo "overlapSelect $OUTDIR/$EXPDIR/altSplice.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/retroSplice.bed"
 overlapSelect -inCoordCols=0,1,2,5,3 $OUTDIR/$EXPDIR/altSplice.bed $OUTDIR/$TABLE.bed $OUTDIR/$EXPDIR/retroSplice.bed
fi

echo '-------- END script analyseExpress.sh -------------------'
echo '------ run script analyzeAge.sh DEF to compute age or retro -------------------'
echo '------ run script ucscRetroStep6.sh DEF to make Html pages -------------------'
