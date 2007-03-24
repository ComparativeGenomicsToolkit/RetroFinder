#
DB=$1
cp -f shuffleHeader.html /usr/local/apache/htdocs/retro/shuffle/index.html;
#cp -f shuffleHeader.html /usr/local/apache/htdocs/retro/rtshuffle/index.html;
/bin/cp -f type1Header.html /usr/local/apache/htdocs/retro/type1/index.html;
cp -f kaessHeader.html /usr/local/apache/htdocs/retro/kass/index.html;
cp -f ageHeader.html /usr/local/apache/htdocs/retro/age/index.html;

makeHtmlRight.sh $DB shufflingAnalysis134.bed shuffle/orig134 orig134;
makeHtmlRight.sh $DB pseudoRefGeneCds.bed shuffle/shuffle refSeqShuffle;
#makeHtmlRight.sh $DB pseudoRefGeneCdsAll.bed shuffle/all refSeqExpressed;
makeHtmlRight.sh $DB pseudoRefGeneCds50.bed shuffle/shuffle refSeq50%Shuffle;
makeHtmlRight.sh $DB shuffleNew.bed shuffle/new refSeqShuffle-orig134;
makeHtmlRight.sh $DB shuffleNew50.bed shuffle/new50 refSeqShuffle50-orig134;
makeHtmlRight.sh $DB shuffleMissed.bed shuffle/miss retroLostFromOrig134;
makeHtmlRight.sh $DB shuffleBig.final.bed shuffle/mrna50 shuffleKgMgcRefSeq50.UTR;
makeHtmlRight.sh $DB shuffleBigKgMgcRefGene.bed shuffle/kgMgcRefGene shuffleKgMgcRefGene;
makeHtmlRight.sh $DB shuffleKgMgcRefNew.bed shuffle/kgMgcRefGene50 shuffleKgMgcRefGene50;
makeHtmlRight.sh $DB shuffleKgMgcRefNew.bed shuffle/kgMgcRefGeneNew shuffleKgMgc50-RefSeq;
makeHtmlRight.sh $DB shuffleKgMgcRefNewNew.bed shuffle/kgMgcRefGeneNewNew shuffleBigKgMgc-RefSeq-orig134;
makeHtmlRight.sh $DB shuffleKgMgcRefLate.bed shuffle/kgMgcRefGeneLate shuffleKgMgcRefSeq-211
makeHtmlRight.sh $DB shuffleKgMgcRefRecent.bed shuffle/kgMgcRefGeneRecent shuffleKgMgcRefSeqRecent
#makeHtmlRight.sh $DB type1MrnaNotShuffle.bed type1/notShuffle type1MrnaNotShuffle;
makeHtmlRight.sh $DB pseudoExpressedAll.bed type1/expressedAll expressedAll;
makeHtmlRight.sh $DB pseudoEstMrna.filter.bed type1/estMrna expressedMrnaOrEst;
makeHtmlRight.sh $DB pseudoEstAll.bed type1/est expressedEst;
#makeHtmlRight.sh $DB pseudoEstAll.filter.bed type1/estDedup expressedEstRemoveOverlap;
makeHtmlRight.sh $DB pseudoEst5Mrna.bed type1/est5Mrna expressedEst5orMoreOrMrna;
makeHtmlRight.sh $DB pseudoEst5.bed type1/est5 expressedEst5orMore;
makeHtmlRight.sh $DB pseudoEst5AndMrna.bed type1/est5AndMrna expressedEst5orMoreAndMrna;
makeHtmlRight.sh $DB est5MrnaRefSeq.bed type1/est5AndMrnaRefSeq expressedEst5orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DB pseudoEst10Mrna.bed type1/est10Mrna expressedEst10orMoreOrMrna;
makeHtmlRight.sh $DB pseudoEst10.bed type1/est10 expressedEst10orMore;
makeHtmlRight.sh $DB pseudoEst10AndMrna.bed type1/est10AndMrna expressedEst10orMoreAndMrna;
makeHtmlRight.sh $DB pseudoMrnaNotShuff.bed type1/mrna expressedMrna;
makeHtmlRight.sh $DB est5MrnaNotRefSeq.bed type1/est5MrnaNotRefSeq expressedEst5orMoreAndMrnaNOTRefSeq;
makeHtmlRight.sh $DB estMrnaRefSeq.bed type1/estMrnaRefSeq expressedMrna+EstRefSeq;
makeHtmlRight.sh $DB pseudoRefGeneCdsAll.bed type1/refSeq type1RefSeq;
makeHtmlRight.sh $DB est10MrnaRefSeq.bed type1/est10AndMrnaRefSeq expressedEst10orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DB pseudoExpressed.bed type1/estOrf estGoodOrf;
makeHtmlRight.sh $DB pseudoEst100AA.bed type1/estOrf100 estGoodOrfBigger100aa;
makeHtmlRight.sh $DB pseudoEstOrf100NotRefSeq.bed type1/estOrf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DB pseudoEstOrf100RefSeq.bed type1/estOrf100RefSeq estGoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100.bed type1/estOrf100 est5GoodOrfBigger100aa;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100NotRefSeq.bed type1/est5Orf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100RefSeq.bed type1/est5Orf100RefSeq estGoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DB est5Ens.bed type1/est5Ens est5MrnaEnsemblNotRefSeq;
makeHtmlRight.sh $DB est5NotEns.bed type1/est5NotEns est5MrnaNotEnsemblNotRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100Spliced.bed type1/5estmrna100Spliced est5mrna100Splice;
makeHtmlRight.sh $DB matchKass.bed kass/matchMrna MatchingMrna;
makeHtmlRight.sh $DB matchEst.bed kass/matchEst MatchingEst;
makeHtmlRight.sh $DB stillMiss.bed kass/notK MissingFromMySet;
makeHtmlRight.sh $DB kassMissed.bed kass/kMissed HitsMissedByKaessmannMrna;
makeHtmlRight.sh $DB kassMissedRefSeq.bed kass/kMissedRefSeq HitsMissedByKaessmannRefSeq;
makeHtmlRight.sh $DB pseudoHuman.bed kass/human humanSpecific;
makeHtmlRight.sh $DB pseudoApe.bed kass/ape apeSpecific;
makeHtmlRight.sh $DB apeNotKass.bed kass/apeNotK apeSpecificNotKaess;
 makeHtmlRight.sh $DB hostEst10Mrna.bed type1/hostMulti type1est10MrnaMultiExon
