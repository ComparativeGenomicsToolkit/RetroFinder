#
CHROM=$1
overlapSelect $CHROM.psl ../retroMrnaInfoLessZnf.bed pseudoEst.$CHROM.bed
overlapSelect $CHROM.psl ../retroMrnaInfoLessZnf.bed -idOutput stdout |awk '{print $1}' | sort |uniq> est.$CHROM.id
overlapSelect $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput stdout | sort > stat.$CHROM.out
overlapSelect $CHROM.psl ../retroMrnaInfoLessZnf.bed -statsOutput -aggregate stdout | sort > statagg.$CHROM.out
