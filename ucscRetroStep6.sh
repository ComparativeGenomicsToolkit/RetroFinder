#!/bin/bash
source $1
echo '-------- script ucscRetroStep6.sh make retro webpages for $SPECIES -------------------'
for i in `echo $SPECIES` ; do echo "make html $i"; done


echo "$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCds.bed shuffle/shuffle refSeqShuffle;"
#$SCRIPT/makeHtmlRightDb.sh $DB shufflingAnalysis134.bed shuffle/orig134 orig134;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCds.bed shuffle/shuffle refSeqShuffle;
#$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCdsAll.bed shuffle/all refSeqExpressed;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCds50.bed shuffle/shuffle refSeq50%Shuffle;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleNew.bed shuffle/new refSeqShuffle-orig134;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleNew50.bed shuffle/new50 refSeqShuffle50-orig134;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleMissed.bed shuffle/miss retroLostFromOrig134;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleBig.final.bed shuffle/mrna50 shuffleKgMgcRefSeq50.UTR;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleBigKgMgcRefGene.bed shuffle/kgMgcRefGene shuffleKgMgcRefGene;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleKgMgcRefNew.bed shuffle/kgMgcRefGene50 shuffleKgMgcRefGene50;
$SCRIPT/makeHtmlRightDb.sh $DB shuffleKgMgcRefNew.bed shuffle/kgMgcRefGeneNew shuffleKgMgc50-RefSeq;
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefNewNew.bed shuffle/kgMgcRefGeneNewNew shuffleBigKgMgc-RefSeq-orig134;
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefLate.bed shuffle/kgMgcRefGeneLate shuffleKgMgcRefSeq-211
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefLate2.bed shuffle/kgMgcRefGeneLate2 shuffleKgMgcRefSeq-211-58
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefLate3.bed shuffle/kgMgcRefGeneLate3 shuffleKgMgcRefSeq-211-58-100
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefRecent.bed shuffle/kgMgcRefGeneNotAncient shuffleKgMgcRefSeqNotAncient
$SCRIPT/makeHtmlRightSortParentDb.sh $DB shuffleKgMgcRefRecentNew.bed shuffle/kgMgcRefGeneNotAncientNew shuffleKgMgcRefSeqNotAncientNew
#$SCRIPT/makeHtmlRightDb.sh $DB type1MrnaNotShuffle.bed type1/notShuffle type1MrnaNotShuffle;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoExpressedAll.bed type1/expressedAll510 expressedAll510;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstMrna.filter.bed type1/estMrna expressedMrnaOrEst510;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstAll.bed type1/est expressedEst510;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstNotShuffle.bed type1/est expressedEstNotShuff;
#$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstAll.filter.bed type1/estDedup expressedEstRemoveOverlap;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoExpressed.bed type1/estOrf estGoodOrf600;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoMrnaNotShuff.bed type1/mrna expressedMrna;
$SCRIPT/makeHtmlRightDb.sh $DB estMrnaRefSeq.bed type1/estMrnaRefSeq expressedMrna+EstRefSeq.score510;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst100AA.bed type1/estOrf100 estGoodOrfBigger100aa.score600;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstOrf100NotRefSeq.bed type1/estOrf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEstOrf100RefSeq.bed type1/estOrf100RefSeq estGoodOrfBigger100aaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5Mrna.bed type1/est5Mrna expressedEst5orMoreOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5.bed type1/est5 expressedEst5orMore;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5AndMrna.bed type1/est5AndMrna expressedEst5orMoreAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5AndMrnaMulti.bed type1/est5AndMrnaMulti expressedEst5orMoreAndMrnaMulti;
$SCRIPT/makeHtmlRightDb.sh $DB est5MrnaRefSeq.bed type1/est5AndMrnaRefSeq expressedEst5orMoreAndMrnaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB est5MrnaNotRefSeq.bed type1/est5MrnaNotRefSeq expressedEst5orMoreAndMrnaNOTRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCdsAll.bed type1/refSeq type1RefSeqCoding;
$SCRIPT/makeHtmlRightDb.sh $DB est10MrnaRefSeq.bed type1/est10AndMrnaRefSeq expressedEst10orMoreAndMrnaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5MrnaOrf100.bed type1/estOrf100 est5GoodOrfBigger100aa;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5MrnaOrf100NotRefSeq.bed type1/est5Orf100NotRefSeq est5GoodOrfBigger100aaNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5MrnaOrf100RefSeq.bed type1/est5Orf100RefSeq est5GoodOrfBigger100aaRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB est5Ens.bed type1/est5Ens est5MrnaEnsemblNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB est5NotEns.bed type1/est5NotEns est5MrnaNotEnsemblNotRefSeq;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst5MrnaOrf100Spliced.bed type1/5estmrna100Spliced est5mrna100Splice;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst10Mrna.bed type1/est10Mrna expressedEst10orMoreOrMrna;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst10.bed type1/est10 expressedEst10orMore;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoEst10AndMrna.bed type1/est10AndMrna expressedEst10orMoreAndMrna;
$SCRIPT/makeHtmlRightDb.sh $DB matchKass.bed kass/matchMrna MatchingMrna;
$SCRIPT/makeHtmlRightDb.sh $DB matchEst.bed kass/matchEst MatchingEst;
$SCRIPT/makeHtmlRightDb.sh $DB stillExpMiss.bed kass/notK MissingFromMyExpSet;
$SCRIPT/makeHtmlRightDb.sh $DB stillMiss.bed kass/notK MissingFromMyExpSet;
$SCRIPT/makeHtmlRightDb.sh $DB kassMissed.bed kass/kMissed HitsMissedByKaessmannMrna;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoRefGeneCdsAll.bed kass/kMissedRefSeq HitsMissedByKaessmannRefSeq;
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB notMatchKass.bed /usr/local/apache/htdocs/retro/shuffle/notK
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB pseudoRefGeneCds.bed /usr/local/apache/htdocs/retro/shuffle/all 
$SCRIPT/makeHtmlRightDb.sh $DB pseudoHuman.bed kass/human humanSpecific;
$SCRIPT/makeHtmlRightDb.sh $DB pseudoApe.bed kass/ape apeSpecific;
$SCRIPT/makeHtmlRightDb.sh $DB apeNotKass.bed kass/apeNotK apeSpecificNotKaess;

#cp  finalType1Header.html /usr/local/apache/htdocs/retro/final/index.html
$SCRIPT/makeHtmlRightDb.sh $DB table2.long.bed final/type1 type1shuffleTable2
echo '-------- END script analyseExpress.sh -------------------'
