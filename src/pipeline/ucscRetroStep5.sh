#!/bin/bash 
# post process cluster jobs that runs retroFinder
# and then filters out zinc fingers and immunoglobin retros if PFAM is defined in the DEF file.
# the retro results, alignments and cds coding of the parent mapped to the retro is loaded into the database.
#
# script analyzeExpress.sh is called at the end to check for expression level of retros.
#
DEF=$1
source $DEF
set -beEu -o pipefail
CWD=`pwd`
## called from ~baertsch/baertsch/scripts/ucscRetroStep4.sh
echo "-----------------------------------"
echo "Starting ucscRetroStep5.sh $DEF on $HOST"
echo "-----------------------------------"
date
# cd $OUTDIR
pwd
# Cat the BED files with overlaps removed to pseudoGeneLinkNoOverlap.bed
echo cat $OUTDIR/chr*_NoOverlap.bed to $OUTDIR/pseudoGeneLinkNoOverlap.bed
cat $OUTDIR/chr*_NoOverlap.bed | $SCRIPT/removeTandemDups > $OUTDIR/pseudoGeneLinkNoOverlap.bed
# Get number of output files from bedOverlap
wc -l $OUTDIR/chr*_NoOverlap.bed
# Create four files based on filtering retrogene score. Filter to keep only
# retrogenes scoring 350 or higher - pseudoGeneLinkNoOverlapFilter.bed. Filter
# to keep only retrogenes scoring 425 or higher - pseudogGeneLink425.bed.
# Filter to keep retrogenes scoring >= 510 - retroMrnaInfo.raw.bed. Filter to
# keep retrogenes scoring >= 650 - retroMrnaInfo650.bed. 
tawk '$5>=350{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/pseudoGeneLinkNoOverlapFilter.bed                       
tawk '$5>=425{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/pseudoGeneLink425.bed
tawk '$5>=510{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/retroMrnaInfo.raw.bed
tawk '$5>=650{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/retroMrnaInfo650.bed
#cut -f 1-12 retroMrnaInfo650.bed > retroMrnaInfo.12.bed
# Count the number of lines in each of the above files - original 
# (pseudoGeneLinkNoOverlap.bed) and filtered.
wc -l $OUTDIR/pseudoGeneLinkNoOverlap.bed $OUTDIR/pseudoGeneLinkNoOverlapFilter.bed $OUTDIR/pseudoGeneLink425.bed $OUTDIR/retroMrnaInfo.raw.bed $OUTDIR/retroMrnaInfo650.bed 

# If a PFAM table is specified in the DEF file, e.g. knownToPfam, then use
# this to remove retrogenes predicted for very large gene families such 
# as immunoglobulins and zinc finger genes and NBPF (DNA-binding TF). 

# Get the known genes (UCSC Genes) that have gene names that start with 
# "IGH", "IGK" or "IGL". Get the GenBank sequences whose gene name is 
# like "IGH", "IKL" or "IGL". Also select name from knownGene where proteinID
# is like "IGH" or "IKL" or "IGL".  
if [[ -n $PFAM ]] 
then
   #remove immunoglobulins
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where (a1.alias like 'IGH%' or a1.alias like 'IGK%' or a1.alias like 'IGL%') and a1.kgID = a2.kgID"  > $OUTDIR/kgImmuno.lst
    hgsql $DB -N -B -e "     select gbCdnaInfo.acc
          from gbCdnaInfo, geneName, description 
             where gbCdnaInfo.geneName=geneName.id and gbCdnaInfo.description = description.id and (geneName.name like 'IGH%' or geneName.name like 'IKL%' or geneName.name like 'IGL%')" >> $OUTDIR/kgImmuno.lst
    hgsql $DB -N -B -e "select name from knownGene where proteinID like 'IGH%' or proteinID like 'IKL%' or proteinID like 'IGL%'" >> $OUTDIR/kgImmuno.lst
#remove znf and NBPF
# Get other aliases from kgAlias where alias is like "ZNF". 
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'ZNF%' and a1.kgID = a2.kgID" |sort -u > $OUTDIR/kgZnf.lst
# Get other aliases from kgAlias where alias is like "NBPF". 
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'NBPF%' and a1.kgID = a2.kgID" |sort -u >> $OUTDIR/kgZnf.lst

#    cat kgZnf.lst refZnf.lst kgImmuno.lst > bothZnf.lst
# Concatenate the zinc finger gene list and the immunoglobulin gene list 
# into a file, bothZnf.lst. 
cat $OUTDIR/kgZnf.lst $OUTDIR/kgImmuno.lst > $OUTDIR/bothZnf.lst

# grep genes with pfam domains (zinc finger, immunoglobin, NBPF, and olfactory receptor
# Get UCSC Genes (known gene) that have PFAM domains for zinc finger, 
# immunoglobulin, NBPF and olfactory receptor genes. 
    hgsql $DB -N -B -e "select  k.name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, \
                exonStarts, exonEnds from $GENEPFAM k, $PFAM p \
                where k.name = p.$PFAMIDFIELD and p.$PFAMDOMAIN in (\
            'PF00096', 'PF01352', 'PF06758', 'PF00047', 'PF07654', 'PF00642'  );" > $OUTDIR/zincKg.gp
    echo zcat $OUTDIR/$GENEPFAM.tab.gz to $OUTDIR/$GENEPFAM.tab
# unzip the knownGene.tab.gz file to knownGene.tab
    zcat $OUTDIR/$GENEPFAM.tab.gz > $OUTDIR/$GENEPFAM.tab
# Use overlapSelect with the zincKg.gp file as the select file and 
# knownGene.tab as the infile and resulting records are in zincKg2.gp
# So this is extracting all known genes that overlap those genes with
# immunoglobulin, zinc finger, NBPF and olfactory receptor domains.
# NOTE: Since these are being selected from the knownGene table, wouldn't you 
# already get all those transcripts encoding proteins that have these domains. 
# Also, no strand is specified with overlapSelect so any genes overlapping on 
# the opposite strand are also being picked up.  
    overlapSelect $OUTDIR/zincKg.gp $OUTDIR/$GENEPFAM.tab $OUTDIR/zincKg2.gp -inFmt=genePred
# Select the name field from the zincKg2.gp file and output to zincKg.lst
    cut -f 1 $OUTDIR/zincKg2.gp | sort | uniq > $OUTDIR/zincKg.lst
    echo "Zinc fingers excluded"
# Get number of genes excluded. 
    wc -l $OUTDIR/zincKg2.gp $OUTDIR/zincKg.lst

# Use the selectById script (column numbers are 1-based), select id in 
# zincKg.lst from column 1 and use to it to select rows by the id in column 1 
# of knownGene.tab and output the resulting genePred records to kgZnf.gp.
# NOTE: Shouldn't the kgZnf.gp file be the same as zincKg.gp or zincKg2.gp?
    echo "$SCRIPT/selectById 1 $OUTDIR/zincKg.lst 1 $OUTDIR/$GENEPFAM.tab to $OUTDIR/kgZnf.gp"
    $SCRIPT/selectById 1 $OUTDIR/zincKg.lst 1 $OUTDIR/$GENEPFAM.tab > $OUTDIR/kgZnf.gp
else
echo "xxx" > $OUTDIR/zincKg.lst
echo "Skipping zinc finger and immunoglobin filtering"
fi

# Use this gene list, zincKg.lst, to extract those genes from the 
# retroMrnaInfo.raw.bed file. Use column 1 of zincKg.lst as id to select from 
# column 44 of retroMrnaInfo.raw.bed and print those rows to 
# retroMrnaInfoZnf.bed.
$SCRIPT/selectById 1 $OUTDIR/zincKg.lst 44 $OUTDIR/retroMrnaInfo.raw.bed > $OUTDIR/retroMrnaInfoZnf.bed
# Using ids from column 1 of zincKg.lst, select those rows that do not have
# this id in column 44 of retroMrnaInfo.raw.bed and output to 
# retroMrnaInfoLessZnF.bed. 
$SCRIPT/selectById -not 1 $OUTDIR/zincKg.lst 44 $OUTDIR/retroMrnaInfo.raw.bed > $OUTDIR/retroMrnaInfoLessZnf.bed
echo "Before and after zinc finger filtering"
rm -f $OUTDIR/$TABLE.bed
# Copy file retroMrnaInfoLessZnf.bed to ucscRetroInfo$VERSION.bed
cp -pf $OUTDIR/retroMrnaInfoLessZnf.bed $OUTDIR/$TABLE.bed
# Sort ucscRetroInfo$VERSION.bed file on the name field and output to 
# ucscRetroInfo$VERSION.sort.bed. 
sort -k4,4 $OUTDIR/$TABLE.bed > $OUTDIR/$TABLE.sort.bed
# Sort on the 4th field of ortho.txt (there is no 4th field) and pipe into 
# join on the 4th field of the sorted file, ucscRetroInfo$VERSION.sort.bed 
# with the 1st field of the sorted ortho.txt file (i.e. join on ids) with 
# output format -o 2.1, 2.2, 2.3 which is FILENUMBER.FIELD i.e the first 3 
# fields of file number 2 which is the sorted ortho.txt and pipe the output 
# of the join into awk and output a tab-separated file inot ortho.filter.txt.  
sort -k4,4 $OUTDIR/ortho.txt | join -1 4 -2 1 -o 2.1 2.2 2.3 $OUTDIR/$TABLE.sort.bed -|awk '{$1=$1;OFS="\t";print $0}' > $OUTDIR/ortho.filter.txt

wc -l $OUTDIR/retroMrnaInfo.raw.bed $OUTDIR/retroMrnaInfoLessZnf.bed $OUTDIR/retroMrnaInfoZnf.bed $OUTDIR/$TABLE.bed
cut -f 1-12 $OUTDIR/$TABLE.bed > $OUTDIR/retroMrnaInfo.12.bed
textHistogram -col=5 $OUTDIR/$TABLE.bed -binSize=50 -maxBinCount=50
echo Creating $OUTDIR/$ALIGN.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' $OUTDIR/$TABLE.bed > $OUTDIR/pseudoGeneLinkSelect.tab
pslSelect -qtStart=$OUTDIR/pseudoGeneLinkSelect.tab $OUTDIR/pseudo.psl $OUTDIR/$ALIGN.psl
wc -l $OUTDIR/$ALIGN.psl $OUTDIR/pseudoGeneLinkSelect.tab
hgLoadBed $DB -verbose=9 -renameSqlTable -allowNegativeScores -noBin ucscRetroInfoXX -sqlTable=$KENTDIR/src/hg/lib/ucscRetroInfo.sql $OUTDIR/$TABLE.bed
mkdir -p $RETRODIR
rm -f $RETRODIR/$TABLE.bed
cp -p $OUTDIR/$TABLE.bed $RETRODIR
hgsql $DB -e "drop table if exists $TABLE;"
hgsql $DB -e "alter table ucscRetroInfoXX rename $TABLE;"
hgLoadPsl $DB $OUTDIR/$ALIGN.psl
rm -f $RETRODIR/$ALIGN.psl
cp -p $OUTDIR/$ALIGN.psl $RETRODIR
hgLoadSqlTab $DB ${ORTHOTABLE} $KENTDIR/src/hg/lib/ucscRetroOrtho.sql $OUTDIR/ortho.filter.txt

zcat $OUTDIR/cds.tab.gz |tawk '{print $1"."$2,$3}' | sort | uniq > $OUTDIR/ucscRetroCds${VERSION}.tab 
hgLoadSqlTab $DB ucscRetroCds${VERSION} $KENTDIR/src/hg/lib/ucscRetroCds.sql $OUTDIR/ucscRetroCds${VERSION}.tab
rm -f $RETRODIR/ucscRetroCds${VERSION}.tab
cp -p $OUTDIR/ucscRetroCds${VERSION}.tab $RETRODIR
cp -p $OUTDIR/DEF $RETRODIR

# Get data for count table and create a new table, ucscRetroCountXX.
hgsql $DB -Ne "create table ucscRetroCount${VERSION} select geneSymbol, count(*) as retroCount, avg(score) as averageScore from ucscRetroInfo${VERSION}, kgXref where kgName = kgID and kgName <> 'noKg' and score > 650  group by geneSymbol having count(*) > 1 order by 2 desc"

#writing TrackDb.ra entry to temp file
$SCRIPT/makeTrackDb.sh $OUTDIR/$DEF > $OUTDIR/trackDb.retro
echo "Writing template trackDb.ra entry to $OUTDIR/trackDb.retro"

echo "Database loaded, update trackDb.ra entry"
echo "run $SCRIPT/analyseExpress.sh $DEF to update expression level of each retro"
