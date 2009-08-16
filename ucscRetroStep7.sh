#calc expression 
source $1
#hgsql $DB -N -B -e "select * from gnfAtlas2Distance" > gnfAtlas2Distance.tab
overlapSelect knownGene.tab.gz retroMrnaInfo650.bed stdout -selectFmt=genePred -idOutput |sort >ucscRetroKnownGene.id 
sort -k4,4 retroMrnaInfo650.bed > retroMrnaInfo650.sort.bed
join -1 1 -2 4 -o 1.1 1.2 2.47 ucscRetroKnownGene.id retroMrnaInfo650.sort.bed |awk '$2!=$3{OFS="\t";$1=$1;print}'> ucscRetroExpression.tab
join -1 1 -2 4 -o 1.1 1.2 2.47 ucscRetroKnownGene.id retroMrnaInfo650.sort.bed |awk '$2!=$3{print "select * from gnfAtlas2Distance where query = \"" $2 "\" and target = \"" $3"\";"}' > ucscRetroKnownGene.sql

