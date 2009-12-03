#!/bin/bash
set -beEu -o pipefail
source $1
DEF=$1
echo "-------- script ucscRetroStep6.sh make retro webpages for $SPECIES -------------------"
#extract extra columns for html pages
$SCRIPT/makeRetroExtraAttr.sh $DEF $OUTDIR/$EXPDIR
for i in `echo $SPECIES` ; do echo "make html $i"; done

cp $DEF $OUTDIR/$EXPDIR
cd $OUTDIR/$EXPDIR
#web pages for shuffling events
mkdir -p $ROOTDIR/retro/shuffle2009
cp $SCRIPT/header.html $ROOTDIR/retro/shuffle2009/index.html
echo "<TR><TH>data set</TH>" >> $ROOTDIR/retro/shuffle2009/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db Ancient</TH><TH>$db Recent</TH>" >> $ROOTDIR/retro/shuffle2009/index.html ; done
echo "</TR>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</THEAD><TBODY>" >> $ROOTDIR/retro/shuffle2009/index.html
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo${GENE1}Cds.bed shuffle2009/kgshuffle ${GENE1}Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo${GENE1}Cds50.bed shuffle2009/kgshuffle50 ${GENE1}50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo${GENE2}Cds.bed shuffle2009/refseqshuffle ${GENE2}Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo${GENE2}Cds50.bed shuffle2009/refseqshuffle50 ${GENE2}50%Shuffle;
echo "</tbody>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</table>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</body>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</html>" >> $ROOTDIR/retro/shuffle2009/index.html

#web pages for duplication (type1) events
mkdir -p $ROOTDIR/retro/type1_2009
cp $SCRIPT/header.html $ROOTDIR/retro/type1_2009/index.html
echo "<TR><TH>data set</TH>" >> $ROOTDIR/retro/type1_2009/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db Ancient</TH><TH>$db Recent</TH>" >> $ROOTDIR/retro/type1_2009/index.html ; done
echo "</TR>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</THEAD><TBODY>" >> $ROOTDIR/retro/type1_2009/index.html
for i in Coding Noncoding Nearcoding Antisense ; do $SCRIPT/makeHtmlRightDb.sh $DEF retroKg${i}.bed type1_2009/known${i} known${i}; done
$SCRIPT/makeHtmlRightDb.sh $DEF ../retroMrnaInfoLessZnf.bed type1_2009/all510 expressedAndNonExpressedScore510;
$SCRIPT/makeHtmlRightDb.sh $DEF ../retroMrnaInfo650.bed type1_2009/all650 expressedAndNonExpressedScore650;
#$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressedAll.bed type1_2009/expressedAll510 expressedAll510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstMrna.filter.bed type1_2009/estMrna atLeastOneEstOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstAll.bed type1_2009/est atLeastOneEst;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressed.bed type1_2009/estOrf atLeastOneEstGoodOrf
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5Mrna.bed type1_2009/est5Mrna 5EstOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5.bed type1_2009/est5 5Est;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrna.bed type1_2009/5estAndMrna 5EstAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst100AA.bed type1_2009/estOrf100 1EstGoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo5Est100AA.bed type1_2009/5estOrf100 5estGoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrna100AA.bed type1_2009/5estAndMrnaOrf100 5estAndMrnaGoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF retroArray.bed type1_2009/gnfAtlas2 gnfAtlas2;
$SCRIPT/makeHtmlRightDb.sh $DEF retroSplice.bed type1_2009/retroSplice Splice;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGene.bed type1_2009/refSeqAll refSeqAll;
echo "</tbody>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</table>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</body>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</html>" >> $ROOTDIR/retro/type1_2009/index.html

echo " add following entry to trackDb.ra"
echo ""
echo "track ucscRetroAli$VERSION"
echo "shortLabel Retroposed Genes $VERSION.0"
echo "longLabel Retroposed Genes, Including Pseudogenes August 2009 [ucscRetroAli$VERSION]"
echo "group genes"
echo "type psl"
echo "priority 37.14"
echo "color 20,0,250"
echo "visibility pack"
echo "nextItemButton on"
echo "retroMrnaInfo ucscRetroInfo$VERSION"
echo "baseColorDefault diffCodons"
echo "baseColorUseCds table ucscRetroCds"
echo "baseColorUseSequence extFile ucscRetroSeq ucscRetroExtFile"
echo "indelDoubleInsert on"
echo "indelQueryInsert on"
echo "showDiffBasesAllScales ."
echo "showDiffBasesMaxZoom 10000.0"
echo "showCdsAllScales ."
echo "showCdsMaxZoom 10000.0"


echo '-------- END script ucscRetroStep6.sh -------------------'
