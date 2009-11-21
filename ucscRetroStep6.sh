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
echo "$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds.bed shuffle2009/shuffle refSeqShuffle;"
$SCRIPT/makeHtmlRightDb.sh $DEF retroShuffleGood.bed shuffle2009/shuffleGood shuffleGood;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds.bed shuffle2009/shuffle refSeqShuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds50.bed shuffle2009/shuffle50 refSeq50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleEns.bed shuffle2009/shuffleEnd shuffleEns;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleEnsMulti.bed shuffle2009/shuffleEnsMulti shuffleEnsMulti;
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
$SCRIPT/makeHtmlRightDb.sh $DEF ../retroMrnaInfo650.bed type1_2009/all all;
#$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressedAll.bed type1_2009/expressedAll510 expressedAll510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstMrna.filter.bed type1_2009/estMrna expressedMrnaOrEst510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstAll.bed type1_2009/est expressedEst510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressed.bed type1_2009/estOrf estGoodOrf600;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5Mrna.bed type1_2009/est5Mrna expressedEst5orMoreOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5.bed type1_2009/est5 expressedEst5orMore;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrna.bed type1_2009/est5AndMrna expressedEst5orMoreAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst100AA.bed type1_2009/estOrf100 estGoodOrfBigger100aa.score600;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudo5Est100AA.bed type1_2009/5estOrf100 5estGoodOrfBigger100aa.score600;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrna100AA.bed type1_2009/5estAndMrnaOrf100 5estAndMrnaGoodOrfBigger100aa.score600;
$SCRIPT/makeHtmlRightDb.sh $DEF retroArray.bed type1_2009/gnfAtlas2 gnfAtlas2;
$SCRIPT/makeHtmlRightDb.sh $DEF retroSplice.bed type1_2009/retroSplice Splice;
echo "</tbody>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</table>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</body>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</html>" >> $ROOTDIR/retro/type1_2009/index.html

echo '-------- END script ucscRetroStep6.sh -------------------'
