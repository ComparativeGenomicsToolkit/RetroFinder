#
CHROM=$1
overlapSelect -inCoorCols=chromCol,startCol,endCol,strandCol,name $CHROM.psl ../retroMrnaInfoLessZnf.bed pseudoEst.$CHROM.bed
overlapSelect -inCoorCols=chromCol,startCol,endCol,strandCol,name $CHROM.psl ../retroMrnaInfoLessZnf.bed -idOutput stdout |awk '{print $1}' | sort |uniq> est.$CHROM.id
overlapSelect -inCoorCols=chromCol,startCol,endCol,strandCol,name $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput stat.$CHROM.out
overlapSelect -inCoorCols=chromCol,startCol,endCol,strandCol,name $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput -aggregate statagg.$CHROM.out
