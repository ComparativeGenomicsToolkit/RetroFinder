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
echo cat $OUTDIR/chr*_NoOverlap.bed to $OUTDIR/pseudoGeneLinkNoOverlap.bed
cat $OUTDIR/chr*_NoOverlap.bed | $SCRIPT/removeTandemDups > $OUTDIR/pseudoGeneLinkNoOverlap.bed
wc -l $OUTDIR/chr*_NoOverlap.bed
tawk '$5>=350{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/pseudoGeneLinkNoOverlapFilter.bed                       
tawk '$5>=425{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/pseudoGeneLink425.bed
tawk '$5>=510{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/retroMrnaInfo.raw.bed
tawk '$5>=650{print $0}' $OUTDIR/pseudoGeneLinkNoOverlap.bed > $OUTDIR/retroMrnaInfo650.bed
#cut -f 1-12 retroMrnaInfo650.bed > retroMrnaInfo.12.bed
wc -l $OUTDIR/pseudoGeneLinkNoOverlap.bed $OUTDIR/pseudoGeneLinkNoOverlapFilter.bed $OUTDIR/pseudoGeneLink425.bed $OUTDIR/retroMrnaInfo.raw.bed $OUTDIR/retroMrnaInfo650.bed 

if [[ -n $PFAM ]] 
then
   #remove immunoglobulins
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where (a1.alias like 'IGH%' or a1.alias like 'IGK%' or a1.alias like 'IGL%') and a1.`kgID = a2.kgID"  > $OUTDIR/kgImmuno.lst
    hgsql $DB -N -B -e "     select gbCdnaInfo.acc
          from gbCdnaInfo, geneName, description 
             where gbCdnaInfo.geneName=geneName.id and gbCdnaInfo.description = description.id and (geneName.name like 'IGH%' or geneName.name like 'IKL%' or geneName.name like 'IGL%')" >> $OUTDIR/kgImmuno.lst
    hgsql $DB -N -B -e "select name from knownGene where proteinID like 'IGH%' or proteinID like 'IKL%' or proteinID like 'IGL%'" >> $OUTDIR/kgImmuno.lst
#remove znf and NBPF
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'ZNF%' and a1.kgID = a2.kgID" |sort -u > $OUTDIR/kgZnf.lst
    hgsql $DB -N -B -e "select a2.alias from kgAlias a1, kgAlias a2 where a1.alias like 'NBPF%' and a1.kgID = a2.kgID" |sort -u >> $OUTDIR/kgZnf.lst

#    cat kgZnf.lst refZnf.lst kgImmuno.lst > bothZnf.lst
cat $OUTDIR/kgZnf.lst $OUTDIR/kgImmuno.lst > $OUTDIR/bothZnf.lst

# grap genes with pfam domains (zinc finger, immunoglobin, NBPF, and olfactory receptor
    hgsql $DB -N -B -e "select  k.name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, \
                exonStarts, exonEnds from $GENEPFAM k, $PFAM p \
                where k.name = p.$PFAMIDFIELD and p.$PFAMDOMAIN in (\
            'PF00096', 'PF01352', 'PF06758', 'PF00047', 'PF07654', 'PF00642'  );" > $OUTDIR/zincKg.gp
    echo zcat $OUTDIR/$GENEPFAM.tab.gz to $OUTDIR/$GENEPFAM.tab
    zcat $OUTDIR/$GENEPFAM.tab.gz > $OUTDIR/$GENEPFAM.tab
    overlapSelect $OUTDIR/zincKg.gp $OUTDIR/$GENEPFAM.tab $OUTDIR/zincKg2.gp -inFmt=genePred
    cut -f 1 $OUTDIR/zincKg2.gp | sort | uniq > $OUTDIR/zincKg.lst
    echo "Zinc fingers excluded"
    wc -l $OUTDIR/zincKg2.gp $OUTDIR/zincKg.lst

    echo "$SCRIPT/selectById 1 $OUTDIR/zincKg.lst 1 $OUTDIR/$GENEPFAM.tab to $OUTDIR/kgZnf.gp"
    $SCRIPT/selectById 1 $OUTDIR/zincKg.lst 1 $OUTDIR/$GENEPFAM.tab > $OUTDIR/kgZnf.gp
else
echo "xxx" > $OUTDIR/zincKg.lst
echo "Skipping zinc finger and immunoglobin filtering"
fi


$SCRIPT/selectById 1 $OUTDIR/zincKg.lst 44 $OUTDIR/retroMrnaInfo.raw.bed > $OUTDIR/retroMrnaInfoZnf.bed
$SCRIPT/selectById -not 1 $OUTDIR/zincKg.lst 44 $OUTDIR/retroMrnaInfo.raw.bed > $OUTDIR/retroMrnaInfoLessZnf.bed
echo "Before and after zinc finger filtering"
rm -f $OUTDIR/$TABLE.bed
cp -pf $OUTDIR/retroMrnaInfoLessZnf.bed $OUTDIR/$TABLE.bed
sort -k4,4 $OUTDIR/$TABLE.bed > $OUTDIR/$TABLE.sort.bed
sort -k4,4 $OUTDIR/ortho.txt | join -1 4 -2 1 -o 2.1 2.2 2.3 $OUTDIR/$TABLE.sort.bed -|awk '{$1=$1;OFS="\t";print $0}' > $OUTDIR/ortho.filter.txt
wc -l $OUTDIR/retroMrnaInfo.raw.bed $OUTDIR/retroMrnaInfoLessZnf.bed $OUTDIR/retroMrnaInfoZnf.bed $OUTDIR/$TABLE.bed
cut -f 1-12 $OUTDIR/$TABLE.bed > $OUTDIR/retroMrnaInfo.12.bed
textHistogram -col=5 $OUTDIR/$TABLE.bed -binSize=50 -maxBinCount=50
echo Creating $OUTDIR/$ALIGN.psl
awk '{printf("%s\t%s\t%s\n", $4,$1,$2)}' $OUTDIR/$TABLE.bed > $OUTDIR/pseudoGeneLinkSelect.tab
pslSelect -qtStart=pseudoGeneLinkSelect.tab $OUTDIR/pseudo.psl $OUTDIR/$ALIGN.psl
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

#writing TrackDb.ra entry to temp file
$SCRIPT/makeTrackDb.sh $OUTDIR/$DEF > $OUTDIR/trackDb.retro
echo "Writing template trackDb.ra entry to $OUTDIR/trackDb.retro"

echo "Database loaded, update trackDb.ra entry"
echo "run $SCRIPT/analyseExpress.sh $DEF to update expression level of each retro"
