#
DEF=$1
INPUT=$2
DIR=$3
NAME=$4
overlapSelect -inCoordCols=0,1,2,5,3 pnasretroExp.from.hg17.bed $INPUT kassMatch.$INPUT 
overlapSelect -inCoordCols=0,1,2,5,3 pnasretroExp.from.hg17.bed $INPUT kassMissed.$INPUT -nonOverlapping
overlapSelect -inCoordCols=0,1,2,5,3 $INPUT pnasretroExp.from.hg17.bed temp.bed -nonOverlapping
overlapSelect -inCoordCols=0,1,2,5,3 shuffleGood.bed temp.bed Missed.$INPUT -nonOverlapping
makeHtmlRight.sh $DEF kassMatch.$INPUT $DIR $NAME.Kass
makeHtmlRight.sh $DEF Missed.$INPUT $DIR.miss $NAME.Kass.miss
wc -l $INPUT kassMatch.$INPUT kassMissed.$INPUT
