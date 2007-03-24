#
DB=$1
INPUT=$2
DIR=$3
NAME=$4
overlapSelect pnasretroExp.from.hg17.bed $INPUT kassMatch.$INPUT 
overlapSelect pnasretroExp.from.hg17.bed $INPUT kassMissed.$INPUT -nonOverlapping
overlapSelect $INPUT pnasretroExp.from.hg17.bed temp.bed -nonOverlapping
overlapSelect shuffleGood.bed temp.bed Missed.$INPUT -nonOverlapping
makeHtmlRight.sh $DB kassMatch.$INPUT $DIR $NAME.Kass
makeHtmlRight.sh $DB Missed.$INPUT $DIR.miss $NAME.Kass.miss
wc -l $INPUT kassMatch.$INPUT kassMissed.$INPUT
