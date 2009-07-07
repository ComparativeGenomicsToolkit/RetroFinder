#!/bin/bash
set -beEu -o pipefail
TARGET=$1
QUERY=$2
COV=$3
ID=$4
TMPOUT=$5.$(hostname)..$$tmp
OUT=$5
S1=$6
S2=$7
FILE=${OUT%.*}
CHROM=${OUT%/*}
LAST=${OUT##*/}
MID=${CHROM##*/}
#echo MID $MID LAST $LAST $MID/$LAST
AXT=/scratch/tmp/retro/$MID/$LAST.axt
mkdir -p /scratch/tmp/retro/$MID
#echo mkdir -p /scratch/tmp/retro/$MID >> lastz.log
#pwd >> lastz.log
#uname -n >> lastz.log
/cluster/bin/penn/x86_64/lastz $TARGET $QUERY --hspthresh=2000 --format=axt > $AXT
axtToPsl $AXT $S1 $S2 $TMPOUT
mv $TMPOUT $OUT
rm -f $AXT
