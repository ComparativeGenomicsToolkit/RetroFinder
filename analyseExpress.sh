#!/bin/bash -e
DB=$1
TABLE=$2
VERSION=$3
echo '-------- script analyseExpress.sh -------------------'
overlapSelect refGeneMultiCds.bed retroMrnaInfo650.bed pseudoRefGeneCds.bed
overlapSelect refGeneMultiCds.bed retroMrnaInfo.12.bed -statsOutput pseudoRefGeneCds.out
overlapSelect refGeneMultiCds.bed retroMrnaInfo.12.bed pseudoRefGeneCds50.bed -overlapThreshold=0.50
#regenerate gene predictions from expressed retros
tawk '$5 > 600 {$2=$2-25; $3=$3+25;print $0}' ../retroMrnaInfo$VERSION.bed > pseudoExpressed.bed
orfBatch $DB pseudoExpressed.bed pseudoExpressed.out pseudoExpressed.gp >borf.out
genePredSingleCover pseudoExpressed.gp pseudoExpressed.single.gp
awk '{print "mrna."$1}' pseudoExpressed.single.gp |sort > pseudoExpressed.ids

overlapSelect estFiltered.psl.gz $TABLE.bed -idOutput stdout | sort > est.id
overlapSelect estFiltered.psl.gz $TABLE.bed -statsOutput stdout | sort > stat.out
overlapSelect all_mrna.psl.gz ucscRetroInfo4.bed -statsOutput stdout |sort > mrna.out
awk '{print $1}' est10Mrna.out |sort |uniq> est10Mrna.id


#makeHtmlRight.sh $DB shufflingAnalysis134.bed shuffle/orig134 orig134;
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
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefNewNew.bed shuffle/kgMgcRefGeneNewNew shuffleBigKgMgc-RefSeq-orig134;
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefLate.bed shuffle/kgMgcRefGeneLate shuffleKgMgcRefSeq-211
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefLate2.bed shuffle/kgMgcRefGeneLate2 shuffleKgMgcRefSeq-211-58
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefLate3.bed shuffle/kgMgcRefGeneLate3 shuffleKgMgcRefSeq-211-58-100
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefRecent.bed shuffle/kgMgcRefGeneNotAncient shuffleKgMgcRefSeqNotAncient
makeHtmlRightSortParent.sh $DB shuffleKgMgcRefRecentNew.bed shuffle/kgMgcRefGeneNotAncientNew shuffleKgMgcRefSeqNotAncientNew
#makeHtmlRight.sh $DB type1MrnaNotShuffle.bed type1/notShuffle type1MrnaNotShuffle;
makeHtmlRight.sh $DB pseudoExpressedAll.bed type1/expressedAll510 expressedAll510;
makeHtmlRight.sh $DB pseudoEstMrna.filter.bed type1/estMrna expressedMrnaOrEst510;
makeHtmlRight.sh $DB pseudoEstAll.bed type1/est expressedEst510;
makeHtmlRight.sh $DB pseudoEstNotShuffle.bed type1/est expressedEstNotShuff;
#makeHtmlRight.sh $DB pseudoEstAll.filter.bed type1/estDedup expressedEstRemoveOverlap;
makeHtmlRight.sh $DB pseudoExpressed.bed type1/estOrf estGoodOrf600;
makeHtmlRight.sh $DB pseudoMrnaNotShuff.bed type1/mrna expressedMrna;
makeHtmlRight.sh $DB estMrnaRefSeq.bed type1/estMrnaRefSeq expressedMrna+EstRefSeq.score510;
makeHtmlRight.sh $DB pseudoEst100AA.bed type1/estOrf100 estGoodOrfBigger100aa.score600;
makeHtmlRight.sh $DB pseudoEstOrf100NotRefSeq.bed type1/estOrf100NotRefSeq estGoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DB pseudoEstOrf100RefSeq.bed type1/estOrf100RefSeq estGoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DB pseudoEst5Mrna.bed type1/est5Mrna expressedEst5orMoreOrMrna;
makeHtmlRight.sh $DB pseudoEst5.bed type1/est5 expressedEst5orMore;
makeHtmlRight.sh $DB pseudoEst5AndMrna.bed type1/est5AndMrna expressedEst5orMoreAndMrna;
makeHtmlRight.sh $DB pseudoEst5AndMrnaMulti.bed type1/est5AndMrnaMulti expressedEst5orMoreAndMrnaMulti;
makeHtmlRight.sh $DB est5MrnaRefSeq.bed type1/est5AndMrnaRefSeq expressedEst5orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DB est5MrnaNotRefSeq.bed type1/est5MrnaNotRefSeq expressedEst5orMoreAndMrnaNOTRefSeq;
makeHtmlRight.sh $DB pseudoRefGeneCdsAll.bed type1/refSeq type1RefSeqCoding;
makeHtmlRight.sh $DB est10MrnaRefSeq.bed type1/est10AndMrnaRefSeq expressedEst10orMoreAndMrnaRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100.bed type1/estOrf100 est5GoodOrfBigger100aa;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100NotRefSeq.bed type1/est5Orf100NotRefSeq est5GoodOrfBigger100aaNotRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100RefSeq.bed type1/est5Orf100RefSeq est5GoodOrfBigger100aaRefSeq;
makeHtmlRight.sh $DB est5Ens.bed type1/est5Ens est5MrnaEnsemblNotRefSeq;
makeHtmlRight.sh $DB est5NotEns.bed type1/est5NotEns est5MrnaNotEnsemblNotRefSeq;
makeHtmlRight.sh $DB pseudoEst5MrnaOrf100Spliced.bed type1/5estmrna100Spliced est5mrna100Splice;
makeHtmlRight.sh $DB pseudoEst10Mrna.bed type1/est10Mrna expressedEst10orMoreOrMrna;
makeHtmlRight.sh $DB pseudoEst10.bed type1/est10 expressedEst10orMore;
makeHtmlRight.sh $DB pseudoEst10AndMrna.bed type1/est10AndMrna expressedEst10orMoreAndMrna;
makeHtmlRight.sh $DB matchKass.bed kass/matchMrna MatchingMrna;
makeHtmlRight.sh $DB matchEst.bed kass/matchEst MatchingEst;
makeHtmlRight.sh $DB stillExpMiss.bed kass/notK MissingFromMyExpSet;
makeHtmlRight.sh $DB stillMiss.bed kass/notK MissingFromMyExpSet;
makeHtmlRight.sh $DB kassMissed.bed kass/kMissed HitsMissedByKaessmannMrna;
makeHtmlRight.sh $DB pseudoRefGeneCdsAll.bed kass/kMissedRefSeq HitsMissedByKaessmannRefSeq;
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB notMatchKass.bed /usr/local/apache/htdocs/retro/shuffle/notK
#~markd/bin/bedToHtmlDir -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu $DB pseudoRefGeneCds.bed /usr/local/apache/htdocs/retro/shuffle/all 
makeHtmlRight.sh $DB pseudoHuman.bed kass/human humanSpecific;
makeHtmlRight.sh $DB pseudoApe.bed kass/ape apeSpecific;
makeHtmlRight.sh $DB apeNotKass.bed kass/apeNotK apeSpecificNotKaess;

cp  finalType1Header.html /usr/local/apache/htdocs/retro/final/index.html
makeHtmlRight.sh $DB table2.long.bed final/type1 type1shuffleTable2
echo '-------- END script analyseExpress.sh -------------------'
