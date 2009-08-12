#!/bin/bash
set -beEu -o pipefail
source $1
DEF=$1
echo "-------- script ucscRetroStep6.sh make retro webpages for $SPECIES -------------------"
#make extract extra columns for html pages
$SCRIPT/makeRetroExtraAttr.sh $DEF $EXPDIR
for i in `echo $SPECIES` ; do echo "make html $i"; done

cp $DEF $EXPDIR
cd $EXPDIR
#web pages for shuffling events
cp $SCRIPT/header.html $ROOTDIR/retro/shuffle2009/index.html
echo "<TR><TH>data set</TH>" >> $ROOTDIR/retro/shuffle2009/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db</TH>" >> $ROOTDIR/retro/shuffle2009/index.html ; done
echo "</TR>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</THEAD><TBODY>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds.bed shuffle2009/shuffle refSeqShuffle;"
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds.bed shuffle2009/shuffle refSeqShuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds50.bed shuffle2009/shuffle50 refSeq50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleEns.bed shuffle2009/shuffleEnd shuffleEns;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleEnsMulti.bed shuffle2009/shuffleEnsMulti shuffleEnsMulti;
echo "</tbody>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</table>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</body>" >> $ROOTDIR/retro/shuffle2009/index.html
echo "</html>" >> $ROOTDIR/retro/shuffle2009/index.html

#web pages for type1 events
cp $SCRIPT/header.html $ROOTDIR/retro/type1_2009/index.html
echo "<TR><TH>data set</TH>" >> $ROOTDIR/retro/type1_2009/index.html
for db in `echo $SPECIES` ; do echo "<TH>$db</TH>" >> $ROOTDIR/retro/type1_2009/index.html ; done
echo "</TR>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</THEAD><TBODY>" >> $ROOTDIR/retro/type1_2009/index.html
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressedAll.bed type1_2009/expressedAll510 expressedAll510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstMrna.filter.bed type1_2009/estMrna expressedMrnaOrEst510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstAll.bed type1_2009/est expressedEst510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoExpressed.bed type1_2009/estOrf estGoodOrf600;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5Mrna.bed type1_2009/est5Mrna expressedEst5orMoreOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5.bed type1_2009/est5 expressedEst5orMore;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrna.bed type1_2009/est5AndMrna expressedEst5orMoreAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst100AA.bed type1_2009/estOrf100 estGoodOrfBigger100aa.score600;
$SCRIPT/makeHtmlRightDb.sh $DEF retroArray.bed type1_2009/gnfAtlas2 gnfAtlas2;
$SCRIPT/makeHtmlRightDb.sh $DEF retroSplice.bed type1_2009/retroSplice Splice;
echo "</tbody>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</table>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</body>" >> $ROOTDIR/retro/type1_2009/index.html
echo "</html>" >> $ROOTDIR/retro/type1_2009/index.html
exit
#$SCRIPT/makeHtmlRightDb.sh $DEF shufflingAnalysis134.bed shuffle/orig134 orig134;
#$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCdsAll.bed shuffle/all refSeqExpressed;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCds50.bed shuffle/shuffle refSeq50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleNew.bed shuffle/new refSeqShuffle-orig134;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleNew50.bed shuffle/new50 refSeqShuffle50-orig134;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleMissed.bed shuffle/miss retroLostFromOrig134;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleBig.final.bed shuffle/mrna50 shuffleKgMgcRefSeq50.UTR;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleBigKgMgcRefGene.bed shuffle/kgMgcRefGene shuffleKgMgcRefGene;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleKgMgcRefNew.bed shuffle/kgMgcRefGene50 shuffleKgMgcRefGene50;
$SCRIPT/makeHtmlRightDb.sh $DEF shuffleKgMgcRefNew.bed shuffle/kgMgcRefGeneNew shuffleKgMgc50-RefSeq;
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefNewNew.bed shuffle/kgMgcRefGeneNewNew shuffleBigKgMgc-RefSeq-orig134;
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefLate.bed shuffle/kgMgcRefGeneLate shuffleKgMgcRefSeq-211
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefLate2.bed shuffle/kgMgcRefGeneLate2 shuffleKgMgcRefSeq-211-58
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefLate3.bed shuffle/kgMgcRefGeneLate3 shuffleKgMgcRefSeq-211-58-100
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefRecent.bed shuffle/kgMgcRefGeneNotAncient shuffleKgMgcRefSeqNotAncient
$SCRIPT/makeHtmlRightSortParentDb.sh $DEF shuffleKgMgcRefRecentNew.bed shuffle/kgMgcRefGeneNotAncientNew shuffleKgMgcRefSeqNotAncientNew
#$SCRIPT/makeHtmlRightDb.sh $DEF type1MrnaNotShuffle.bed type1/notShuffle type1MrnaNotShuffle;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstNotShuffle.bed type1/est expressedEstNotShuff;
#$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstAll.filter.bed type1/estDedup expressedEstRemoveOverlap;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoMrnaNotShuff.bed type1/mrna expressedMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF estMrnaRefSeq.bed type1/estMrnaRefSeq expressedMrna+EstRefSeq.score510;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstOrf100NotRefSeq.bed type1/estOrf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEstOrf100RefSeq.bed type1/estOrf100RefSeq estGoodOrfBigger100aaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5AndMrnaMulti.bed type1/est5AndMrnaMulti expressedEst5orMoreAndMrnaMulti;
$SCRIPT/makeHtmlRightDb.sh $DEF est5MrnaRefSeq.bed type1/est5AndMrnaRefSeq expressedEst5orMoreAndMrnaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF est5MrnaNotRefSeq.bed type1/est5MrnaNotRefSeq expressedEst5orMoreAndMrnaNOTRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCdsAll.bed type1/refSeq type1RefSeqCoding;
$SCRIPT/makeHtmlRightDb.sh $DEF est10MrnaRefSeq.bed type1/est10AndMrnaRefSeq expressedEst10orMoreAndMrnaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5MrnaOrf100.bed type1/estOrf100 est5GoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5MrnaOrf100NotRefSeq.bed type1/est5Orf100NotRefSeq est5GoodOrfBigger100aaNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5MrnaOrf100RefSeq.bed type1/est5Orf100RefSeq est5GoodOrfBigger100aaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF est5Ens.bed type1/est5Ens est5MrnaEnsemblNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF est5NotEns.bed type1/est5NotEns est5MrnaNotEnsemblNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst5MrnaOrf100Spliced.bed type1/5estmrna100Spliced est5mrna100Splice;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst10Mrna.bed type1/est10Mrna expressedEst10orMoreOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst10.bed type1/est10 expressedEst10orMore;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoEst10AndMrna.bed type1/est10AndMrna expressedEst10orMoreAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF matchKass.bed kass/matchMrna MatchingMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF matchEst.bed kass/matchEst MatchingEst;
$SCRIPT/makeHtmlRightDb.sh $DEF stillExpMiss.bed kass/notK MissingFromMyExpSet;
$SCRIPT/makeHtmlRightDb.sh $DEF stillMiss.bed kass/notK MissingFromMyExpSet;
$SCRIPT/makeHtmlRightDb.sh $DEF kassMissed.bed kass/kMissed HitsMissedByKaessmannMrna;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoRefGeneCdsAll.bed kass/kMissedRefSeq HitsMissedByKaessmannRefSeq;
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DEF notMatchKass.bed /usr/local/apache/htdocs/retro/shuffle/notK
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DEF pseudoRefGeneCds.bed /usr/local/apache/htdocs/retro/shuffle/all 
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoHuman.bed kass/human humanSpecific;
$SCRIPT/makeHtmlRightDb.sh $DEF pseudoApe.bed kass/ape apeSpecific;
$SCRIPT/makeHtmlRightDb.sh $DEF apeNotKass.bed kass/apeNotK apeSpecificNotKaess;

#cp  finalType1Header.html /usr/local/apache/htdocs/retro/final/index.html
$SCRIPT/makeHtmlRightDb.sh $DEF table2.long.bed final/type1 type1shuffleTable2
echo '-------- END script analyseExpress.sh -------------------'
