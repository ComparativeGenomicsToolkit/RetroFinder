#
DEF=$1
INPUT=$2
NAME=$3
source $DEF
tawk -f $SCRIPTS/ape.awk $INPUT |sort -u> ape.$INPUT
tawk -f $SCRIPTS/primate.awk $INPUT |sort -u> primate.$INPUT
tawk -f $SCRIPTS/euArch.awk $INPUT |sort -u> euArch.$INPUT
overlapSelect euArch.$INPUT $INPUT tmp1.bed  -nonOverlapping
overlapSelect primate.$INPUT tmp1.bed tmp.bed  -nonOverlapping
overlapSelect ape.$INPUT tmp.bed old.$INPUT -nonOverlapping
#tawk '{print $4,$35,$5,$32,$54,$31,$1":"$2"-"$3}' s100.old.bed|sort -k4,4nr -k6,6nr -k5,5nr
wc -l $INPUT ape.$INPUT primate.$INPUT old.$INPUT
makeHtmlRight.sh $DEF ape.$INPUT age/ape.${INPUT%%.bed} ape.$NAME;
makeHtmlRight.sh $DEF primate.$INPUT age/primate.${INPUT%%.bed} primate.$NAME;
makeHtmlRight.sh $DEF euArch.$INPUT age/euArch.${INPUT%%.bed} euArch.$NAME;
makeHtmlRight.sh $DEF old.$INPUT age/old.${INPUT%%.bed} old.$NAME;

