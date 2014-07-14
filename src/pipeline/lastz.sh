#!/bin/bash
set -beEu -o pipefail
TARGET=$1
QUERY=$2
COV=$3
ID=$4
OUT=$5
TMPOUT=$5.$(hostname)..$$tmp
AXTOUT=$6
S1=$7
S2=$8
# Location of lastz program
LASTZ=$9
FILE=${OUT%.*}
CHROM=${OUT%/*}
LAST=${OUT##*/}
MID=${CHROM##*/}
#echo MID $MID LAST $LAST $MID/$LAST
# $TMPDIR is an environment variable that points at /data/tmp
AXT=$TMPDIR/retro/$MID/$LAST.axt
mkdir -p $TMPDIR/retro/$MID
#echo mkdir -p $TMPDIR/retro/$MID >> lastz.log
#pwd >> lastz.log
#uname -n >> lastz.log
echo "AXT: $AXT AXTOUT: $AXTOUT"
$LASTZ $TARGET $QUERY --ambiguous=iupac --hspthresh=2000 --format=axt > $AXT
axtToPsl $AXT $S1 $S2 $TMPOUT
mv $TMPOUT $OUT
echo "Making directory $AXTOUT"
mkdir -p $AXTOUT
#rm -f $AXT
mv $AXT $AXTOUT
