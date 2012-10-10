#!/bin/bash 
set -bveEu -o pipefail
CHROM=$1
overlapSelect -inCoordCols=0,1,2,5,3 $CHROM.psl ../retroMrnaInfoLessZnf.bed pseudoEst.$CHROM.bed
overlapSelect -inCoordCols=0,1,2,5,3 $CHROM.psl ../retroMrnaInfoLessZnf.bed -idOutput stdout |awk '{print $1}' | sort |uniq> est.$CHROM.id
overlapSelect -inCoordCols=0,1,2,5,3 $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput stat.$CHROM.out
overlapSelect -inCoordCols=0,1,2,5,3 $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput -aggregate statagg.$CHROM.out
