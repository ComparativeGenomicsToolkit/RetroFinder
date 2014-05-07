#!/bin/bash
set -bevEu -x -o pipefail
source $1
DEF=$1
echo "-------- script ucscRetroStep6.sh make retro webpages for $SPECIES -------------------"
#extract extra columns for html pages
$SCRIPT/makeRetroExtraAttr.sh $DEF $OUTDIR/$EXPDIR
for i in `echo $SPECIES` ; do echo "make html $i"; done

cp $DEF $OUTDIR/$EXPDIR
#cd $OUTDIR/$EXPDIR
#web pages for shuffling events
#SHUFFLEDIR=shuffle
# SHUFFLEROOT=$WEBROOT/$SHUFFLEDIR
echo "mkdir -p $SHUFFLEROOT"
mkdir -p $SHUFFLEROOT
cp $SCRIPT/header.html $SHUFFLEROOT/index.html
echo "<TR><TH>data set</TH>" >> $SHUFFLEROOT/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db Recent</TH><TH>$db Ancient</TH>" >> $SHUFFLEROOT/index.html ; done
echo "</TR>" >> $SHUFFLEROOT/index.html
echo "</THEAD><TBODY>" >> $SHUFFLEROOT/index.html
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE1}Cds.bed $SHUFFLEDIR/kgshuffle ${GENE1}Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE1}Cds50.bed $SHUFFLEDIR/kgshuffle50 ${GENE1}50Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE2}Cds.bed $SHUFFLEDIR/refseqshuffle ${GENE2}Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE2}Cds50.bed $SHUFFLEDIR/refseqshuffle50 ${GENE2}50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE1}ProtCds.bed $SHUFFLEDIR/kgshuffleProt ${GENE1}ShuffleProtEvidence;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE2}ProtCds.bed $SHUFFLEDIR/refseqshuffleProt ${GENE2}ShuffleProtEvidence;
echo "</tbody>" >> $SHUFFLEROOT/index.html
echo "</table>" >> $SHUFFLEROOT/index.html
echo "</body>" >> $SHUFFLEROOT/index.html
echo "</html>" >> $SHUFFLEROOT/index.html

#web pages for duplication (type2) events
#DUPDIR=dups
#DUPROOT=$WEBROOT/$DUPDIR
mkdir -p $DUPROOT
cp $SCRIPT/header.html $DUPROOT/index.html
echo "<TR><TH>data set</TH>" >> $DUPROOT/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db Recent</TH><TH>$db Ancient</TH>" >> $DUPROOT/index.html ; done
echo "</TR>" >> $DUPROOT/index.html
echo "</THEAD><TBODY>" >> $DUPROOT/index.html
for i in Coding Noncoding Nearcoding Antisense ; do $SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/retroKg${i}.bed $DUPDIR/known${i} known${i}; done
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/ucscRetroProtein.bed $DUPDIR/SpProtEvidence SpProtEvidence;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/ucscRetroTranscript.bed $DUPDIR/SpTransEvidence SpTransEvidence;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/retroMrnaInfoLessZnf.bed $DUPDIR/all510 expressedAndNonExpressedScore510;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/retroMrnaInfo650.bed $DUPDIR/all650 expressedAndNonExpressedScore650;
#$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressedAll.bed $DUPDIR/expressedAll510 expressedAll510;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEstMrna.filter.bed $DUPDIR/estMrna atLeastOneEstOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEstAll.bed $DUPDIR/est atLeastOneEst;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoExpressed.bed $DUPDIR/estOrf atLeastOneEstGoodOrf
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.bed $DUPDIR/est5Mrna 5EstOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5MrnaNotKg.bed $DUPDIR/est5MrnaNotKg 5EstOrMrnaNotKg;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst10Mrna.bed $DUPDIR/est10Mrna 10EstOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst10MrnaNotKg.bed $DUPDIR/est10MrnaNotKg 10EstOrMrnaNotKg;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5.bed $DUPDIR/est5 5Est;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5AndMrna.bed $DUPDIR/5estAndMrna 5EstAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5AndMrnaNotKg.bed $DUPDIR/5estAndMrnaNotKg 5EstAndMrnaNotKg;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst100AA.bed $DUPDIR/estOrf100 1EstGoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo5Est100AA.bed $DUPDIR/5estOrf100 5estGoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5AndMrna100AA.bed $DUPDIR/5estAndMrnaOrf100 5estAndMrnaGoodOrfBigger100aa;
# retroArray.bed is made by the analyseArray.sh script which I don't run
#$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/retroArray.bed $DUPDIR/gnfAtlas2 gnfAtlas2;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/retroSplice.bed $DUPDIR/retroSplice Splice;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudo${GENE2}.bed $DUPDIR/refSeqAll refSeqAll;

# top50.bed used to be made by analyseExpress.sh but code to create it is now
# commented out. 
# $SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/top50.bed $DUPDIR/top50 top50;
echo "</tbody>" >> $DUPROOT/index.html
echo "</table>" >> $DUPROOT/index.html

echo "</body>" >> $DUPROOT/index.html
echo "</html>" >> $DUPROOT/index.html

#web pages for age breakdown 
# AGEDIR=age
# AGEROOT=$WEBROOT/$AGEDIR
mkdir -p $AGEROOT
cp $SCRIPT/header.html $AGEROOT/index.html
echo "<TR><TH>data set</TH>" >> $AGEROOT/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db Recent</TH><TH>$db Ancient</TH>" >> $AGEROOT/index.html ; done
echo "</TR>" >> $AGEROOT/index.html
echo "</THEAD><TBODY>" >> $AGEROOT/index.html
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/ucscRetroProtein.bed $AGEDIR/SpProtEvidence SpProtEvidence;
# These files are made by the analyseAge.sh script or doAge script called 
# by this script. 
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.0ancient.bed $AGEDIR/ancient ancient;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.1mouse.bed $AGEDIR/mouse mouse;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.2rhesus.bed $AGEDIR/rhesus rhesus;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.3ape.bed $AGEDIR/ape ape;
$SCRIPT/makeHtmlRightDb.sh $DEF $OUTDIR/$EXPDIR/pseudoEst5Mrna.4human.bed $AGEDIR/human humanSpecific;
echo "</tbody>" >> $AGEROOT/index.html
echo "</table>" >> $AGEROOT/index.html

echo "</body>" >> $AGEROOT/index.html
echo "</html>" >> $AGEROOT/index.html


echo '-------- END script ucscRetroStep6.sh -------------------'
