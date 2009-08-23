#i!/bin/bash
set -beEu -o pipefail
set -o verbose
source $1
hgsql $DB -N -B -e "select chrom , chromStart , chromEnd   , name       , score      , strand     , thickStart , thickEnd   , reserved   , blockCount , blockSizes , chromStarts from $ARRAY" > array.bed
echo "SCORE Threshold is $SCORETHRESH"
tawk '$5>'$SCORETHRESH' {print $0}' ../$TABLE.bed | overlapSelect ../pseudoEstAll.bed stdin retroExpress.bed -inFmt=bed
#retroArray has retros with expression data
overlapSelect array.bed retroExpress.bed retroArray.bed 
#retroArrayKnown has retros overlapping knownGenes with expression data
overlapSelect ../knownGene.tab.gz retroArray.bed retroArrayKnown.bed -selectFmt=genePred
wc -l retroArray.bed retroArrayKnown.bed
#
# Tau stat - tissue specificity per probe using absolute expression level, output to gnfAtlas2_abs_median.tab 
#
hgTissueSpecific $DB hgFixed.$ARRAYMEDIAN hgFixed.$ARRAYEXP gnfAtlas2_abs_median -lookup=$ARRAYLOOKUP 
tawk '{print $47,$4}' retroArrayKnown.bed | sort > retroArrayKnown.sort.bed
# retroArrayKnown.sort.bed contains parent kgid and pro
sort gnfAtlas2_abs_median.tab > gnfAtlas2_abs_median.sort.tab
#retroParentTissue adds tau stat to parent of expressed retros
join retroArrayKnown.sort.bed gnfAtlas2_abs_median.sort.tab |awk '{OFS="\t";$1=$1;print $2,$1,$3,$4}' | sort > retroParentTissue.tab
#
# extract tissue specific expression (log ratio)
#
# output is bed with location of expressed retros with gnf expression
rm -f retroTissue.bed
for i in `cat retroArrayKnown.bed|tawk '{print $1":"$2"-"$3}'` ; do hgGetAnn $DB $ARRAY $i stdout ; done |grep -v x_at |grep -v s_at >> retroTissue.bed
#  filters cases with a uniquely mapped probe
tawk 'NF==15{print $0}' retroTissue.bed |sort -u > retroTissue.ok.bed
# map to full retro bed
overlapSelect retroTissue.ok.bed retroArrayKnown.bed stdout -idOutput |sort > retroTissue.tab 
# join gnf expression with tau stat
join -1 1 -1 1 retroParentTissue.tab retroTissue.tab|uniq > retroTissueBoth.tab
#add tau stat
tawk '{print $3,$2}' gnfAtlas2_abs_median.tab |sort > gnfProb.tab
join retroTissueBoth.tab gnfProb.tab
sort -k5 retroTissueBoth.tab  > x
join -1 5 -2 1 x gnfProb.tab|uniq |awk '{OFS="\t";$1=$1;print}'> retroTissueFinal.tab
#
# tissue specific analysis
#
cut -f 1-6,15 retroTissue.ok.bed |tawk '{print $1,$2,$3,$4,$5,$6,","$7}'> retroTissueExtract.bed
hgsql hgFixed -N -B -e "select name from hgFixed.gnfHumanAtlas2MedianExps " |sed -e "s/ /_/g"|sed -e "s/+//g" |sed -e "s/(/_/g"|sed -e "s/)//g" > expName.tab
hgsql hgFixed -N -B -e "select id+2, name from hgFixed.gnfHumanAtlas2MedianExps " |sed -e "s/ /_/g"|sed -e "s/+//g" |sed -e "s/(/_/g"|sed -e "s/)//g" > expName.tab
# id is col offset into tissue table
mkdir -p tissue
tawk '{print "awk -F, xx{print $1, tt\$"$1"}xx retroTissueExtract.bed |tawk xx{\$8\=\(\$8\+5)\*10;print $1,$2,$3,$4,$5,$6,$8}xx > tissue/retroTissue."$2".bed"} ' expName.tab |sed -e "s/xx/'/g" |sed -e 's/tt/"\t"/' > buildTissues.sh
chmod +x buildTissues.sh
./buildTissues.sh
mkdir -p retro.tissue
for i in `ls tissue/*.bed`; do cut -f 1-3 $i | overlapSelect -selectFmt=bed stdin ../$TABLE.bed retro.$i ; done

for t in `cut -f 2 expName.tab` ; do $SCRIPT/makeHtmlRightArray.sh DEF retro.tissue/retroTissue.$t.bed array/tissue/$DB/$t all_$t;done

#for i in `ls retroTissue.*.bed` ; do textHistogram $i -col=7 -autoScale=30 ; done > tissue.hist
wc -l retroTissueExtract.bed 
echo "run arrayCompare.sh for R analysis after all species are complete"
