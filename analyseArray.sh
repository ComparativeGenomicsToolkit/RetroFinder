#!/bin/bash
set -beEu -o pipefail
source $1
hgsql $DB -N -B -e "select chrom , chromStart , chromEnd   , name       , score      , strand     , thickStart , thickEnd   , reserved   , blockCount , blockSizes , chromStarts from $ARRAY" > array.bed
tawk '$5>650{print $0}' ../$TABLE.bed | overlapSelect pseudoExpressed.bed stdin retroExpress.bed -inFmt=bed
overlapSelect array.bed retroExpress.bed retroArray.bed 
overlapSelect ../knownGene.tab.gz retroArray.bed retroArrayKnown.bed -selectFmt=genePred
wc -l retroArray.bed retroArrayKnown.bed
hgTissueSpecific hg18 hgFixed.gnfHumanAtlas2Median hgFixed.gnfHumanAtlas2MedianExps gnfAtlas2_abs_median -lookup=knownToGnfAtlas2 
join retroArrayKnown.sort.bed gnfAtlas2_abs_median.sort.tab |awk '{OFS="\t";$1=$1;print $2,$1,$3,$4}' | sort > retroParentTissue.tab
rm -f retroTissue.bed
for i in `cat retroArrayKnown.bed|tawk '{print $1":"$2"-"$3}'` ; do echo $i ; hgGetAnn hg18 gnfAtlas2 $i stdout ; done |grep -v x_at retro.tab|grep -v s_at >> retroTissue.bed
tawk 'NF==15{print $0}' retroTissue.bed > retroTissue.ok.bed
tawk 'NF!=15{print $0}' retroTissue.bed > retroTissue.noexp.bed
overlapSelect retroTissue.ok.bed retroArrayKnown.bed stdout -idOutput |sort > retroTissue.tab 
join -1 1 -1 1 retroParentTissue.tab retroTissue.tab|uniq > retroTissueBoth.tab
tawk '{print $3,$2}' gnfAtlas2_abs_median.tab |sort > gnfProb.tab
join retroTissueBoth.tab gnfProb.tab
sort -k5 retroTissueBoth.tab  > x
join -1 5 -2 1 x gnfProb.tab|uniq |awk '{OFS="\t";$1=$1;print}'> retroTissueFinal.tab





