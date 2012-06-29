#
DEF=$1
source $DEF
cp -f shuffleHeader.html $WEBROOT/shuffle/index.html;
#cp -f shuffleHeader.html $WEBROOT/rtshuffle/index.html;
/bin/cp -f type1Header.html $WEBROOT/type1/index.html;
cp -f kaessHeader.html $WEBROOT/kass/index.html;
cp -f ageHeader.html $WEBROOT/age/index.html;

makeHtmlRight.sh $DEF shufflingAnalysis134.bed shuffle/orig134 orig134;
makeHtmlRight.sh $DEF pseudoRefGeneCds.bed shuffle/shuffle refSeqShuffle;
#makeHtmlRight.sh $DEF pseudoRefGeneCdsAll.bed shuffle/all refSeqExpressed;
makeHtmlRight.sh $DEF pseudoRefGeneCds50.bed shuffle/shuffle refSeq50%Shuffle;
makeHtmlRight.sh $DEF shuffleNew.bed shuffle/new refSeqShuffle-orig134;
makeHtmlRight.sh $DEF shuffleNew50.bed shuffle/new50 refSeqShuffle50-orig134;
makeHtmlRight.sh $DEF shuffleMissed.bed shuffle/miss retroLostFromOrig134;
makeHtmlRight.sh $DEF shuffleBig.final.bed shuffle/mrna50 shuffleKgMgcRefSeq50.UTR;
makeHtmlRight.sh $DEF shuffleBigKgMgcRefGene.bed shuffle/kgMgcRefGene shuffleKgMgcRefGene;
makeHtmlRight.sh $DEF shuffleKgMgcRefNew.bed shuffle/kgMgcRefGene50 shuffleKgMgcRefGene50;
makeHtmlRight.sh $DEF shuffleKgMgcRefNew.bed shuffle/kgMgcRefGeneNew shuffleKgMgc50-RefSeq;
makeHtmlRight.sh $DEF shuffleKgMgcRefNewNew.bed shuffle/kgMgcRefGeneNewNew shuffleBigKgMgc-RefSeq-orig134;
makeHtmlRight.sh $DEF shuffleKgMgcRefLate.bed shuffle/kgMgcRefGeneLate shuffleKgMgcRefSeq-211
makeHtmlRight.sh $DEF shuffleKgMgcRefRecent.bed shuffle/kgMgcRefGeneRecent shuffleKgMgcRefSeqRecent
#makeHtmlRight.sh $DEF type1MrnaNotShuffle.bed type1/notShuffle type1MrnaNotShuffle;
makeHtmlRight.sh $DEF pseudoExpressedAll.bed type1/expressedAll expressedAll;
makeHtmlRight.sh $DEF pseudoEstMrna.filter.bed type1/estMrna expressedMrnaOrEst;
makeHtmlRight.sh $DEF pseudoEstAll.bed type1/est expressedEst;
#makeHtmlRight.sh $DEF pseudoEstAll.filter.bed type1/estDedup expressedEstRemoveOverlap;
makeHtmlRight.sh $DEF pseudoEst5Mrna.bed type1/est5Mrna expressedEst5orMoreOrMrna;
makeHtmlRight.sh $DEF pseudoEst5.bed type1/est5 expressedEst5orMore;
makeHtmlRight.sh $DEF pseudoEst5AndMrna.bed type1/est5AndMrna expressedEst5orMoreAndMrna;
makeHtmlRight.sh $DEF est5MrnaRefSeq.bed type1/est5AndMrnaRefSeq expressedEst5orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DEF pseudoEst10Mrna.bed type1/est10Mrna expressedEst10orMoreOrMrna;
makeHtmlRight.sh $DEF pseudoEst10.bed type1/est10 expressedEst10orMore;
makeHtmlRight.sh $DEF pseudoEst10AndMrna.bed type1/est10AndMrna expressedEst10orMoreAndMrna;
makeHtmlRight.sh $DEF pseudoMrnaNotShuff.bed type1/mrna expressedMrna;
makeHtmlRight.sh $DEF est5MrnaNotRefSeq.bed type1/est5MrnaNotRefSeq expressedEst5orMoreAndMrnaNOTRefSeq;
makeHtmlRight.sh $DEF estMrnaRefSeq.bed type1/estMrnaRefSeq expressedMrna+EstRefSeq;
makeHtmlRight.sh $DEF pseudoRefGeneCdsAll.bed type1/refSeq type1RefSeq;
makeHtmlRight.sh $DEF est10MrnaRefSeq.bed type1/est10AndMrnaRefSeq expressedEst10orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DEF pseudoExpressed.bed type1/estOrf estGoodOrf;
makeHtmlRight.sh $DEF pseudoEst100AA.bed type1/estOrf100 estGoodOrfBigger100aa;
makeHtmlRight.sh $DEF pseudoEstOrf100NotRefSeq.bed type1/estOrf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DEF pseudoEstOrf100RefSeq.bed type1/estOrf100RefSeq estGoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DEF pseudoEst5MrnaOrf100.bed type1/estOrf100 est5GoodOrfBigger100aa;
makeHtmlRight.sh $DEF pseudoEst5MrnaOrf100NotRefSeq.bed type1/est5Orf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DEF pseudoEst5MrnaOrf100RefSeq.bed type1/est5Orf100RefSeq estGoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DEF est5Ens.bed type1/est5Ens est5MrnaEnsemblNotRefSeq;
makeHtmlRight.sh $DEF est5NotEns.bed type1/est5NotEns est5MrnaNotEnsemblNotRefSeq;
makeHtmlRight.sh $DEF pseudoEst5MrnaOrf100Spliced.bed type1/5estmrna100Spliced est5mrna100Splice;
makeHtmlRight.sh $DEF matchKass.bed kass/matchMrna MatchingMrna;
makeHtmlRight.sh $DEF matchEst.bed kass/matchEst MatchingEst;
makeHtmlRight.sh $DEF stillMiss.bed kass/notK MissingFromMySet;
makeHtmlRight.sh $DEF kassMissed.bed kass/kMissed HitsMissedByKaessmannMrna;
makeHtmlRight.sh $DEF kassMissedRefSeq.bed kass/kMissedRefSeq HitsMissedByKaessmannRefSeq;
makeHtmlRight.sh $DEF pseudoHuman.bed kass/human humanSpecific;
makeHtmlRight.sh $DEF pseudoApe.bed kass/ape apeSpecific;
makeHtmlRight.sh $DEF apeNotKass.bed kass/apeNotK apeSpecificNotKaess;
 makeHtmlRight.sh $DEF hostEst10Mrna.bed type1/hostMulti type1est10MrnaMultiExon
