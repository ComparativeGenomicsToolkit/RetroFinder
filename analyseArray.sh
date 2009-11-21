#i!/bin/bash
set -beEu -o pipefail
set -o verbose
source $1
# get alignments for array data and find list of probes
pslCat $ARRAYPSLS -nohead > affy.psl 
cut -f 10 affy.psl |sort -u > array.id
#
#get list of experiments and header
#
hgsql $DB -N -B -e "select id, replace(name,' ','_') from $ARRAYEXPALL" > expNameAll.tab
echo -n "Gene Name" >  $DB.AbsAll.csv ; for i in `cut -f 2 expNameAll.tab ` ; do echo -n ","$i ; done >>  $DB.AbsAll.csv ; echo "">> $DB.AbsAll.csv
rm -f  $DB.LogRatioAll.csv
cp  $DB.AbsAll.csv  $DB.LogRatioAll.csv
#
#get expression data and filter controls and genes not present in the genome
#
hgsql $DB -N -B -e "select concat(name,',', expScores) from $ARRAYABS" | sed -e 's/,$//'>>  $DB.AbsAll.csv
awk -F, '{OFS="\t";$1=$1;print }'  $DB.AbsAll.csv |grep -F -f array.id >  $DB.AbsAll.tab
hgsql $DB -N -B -e "select concat(name,',', expScores) from $ARRAYRATIO" |grep -F -f array.id | sed -e 's/,$//'>>  $DB.LogRatioAll.csv
awk -F, '{OFS="\t";$1=$1;print }'  $DB.LogRatioAll.csv >  $DB.LogRatioAll.tab
exit
#
#find retros and parent genes with array expression data
#
hgMapMicroarray array.bed $ARRAYRATIO $ARRAYPSLS
grep -v _hap array.bed |grep -v "_random"> array.nohap.bed
tawk '{print $4}' array.nohap.bed|uniq -c|awk '$1==1{print $2}'> array.uniq
$SCRIPT/selectById 1 array.uniq 4 array.nohap.bed |sort -k4 > array.uniq.bed
wc -l array.nohap.bed array.uniq.bed
overlapSelect affy.psl ../$TABLE.bed retroAffy.bed
overlapSelect affy.psl ../$TABLE.bed retroAffy.id -idOutput
$SCRIPT/selectById 1 array.uniq 2 retroAffy.id > retroAffy.uniq.id
$SCRIPT/selectById 2 retroAffy.uniq.id 1  $DB.LogRatioAll.tab > retroAffyLogRatioAll.tab
hgsql $DB -N -B -e "select * from refGene" |cut -f 2-16 > refGene.genePred
overlapSelect retroAffy.bed refGene.genePred retroAffyGenePred.genePred
overlapSelect retroAffy.bed refGene.genePred retroAffyGenePred.id


#covariate file for edge
rm -f cov.$DB.csv
head -1  $DB.AbsAll.csv > cov.$DB.csv
head -1  $DB.AbsAll.csv|sed -e 's/_2//g'|sed -e 's/^Gene Name/Tissue/'|awk -F, '{OFS="\t";print}' >> cov.$DB.csv
awk 'BEGIN{printf("Individual")}{printf(",%d",1-NR%2+1)}END{print ""}' expNameAll.tab >> cov.$DB.csv
awk -F, '{OFS="\t";$1=$1;print }' cov.$DB.csv > cov.$DB.tab



hgsql $DB -N -B -e "select chrom, chromStart, chromEnd, value, substring_index(concat(replace(geneSymbol,' ','.'),'_',spDisplayID),'_',2) as gene, substring(description, 1,30 ) as description from knownToGnfAtlas2 a, knownCanonical k, kgXref x  where a.name = transcript and chrom not like '%_hap%' and chrom not like '%random' and a.name = kgID order by value " > gnfName.bed
cat gnfName.bed|cut -f 4,5 >gnfName.id 
join -1 4 -2 1  array.uniq.bed gnfName.id > array.$DB.txt


hgsql $DB -N -B -e " select distinct value, geneSymbol from kgXref x, knownToGnfAtlas2 k where x.kgID = k.name order by 1" > gnfGeneName.tab
tawk 'BEGIN{p=$1}p==$1{g=g","$2}p!=$1{print p,g;p=$1;g=$2}' gnfGeneName.tab > gnfGeneNameUniq.tab
hgsql $DB -N -B -e " select distinct value , spDisplayID from kgXref x, knownToGnfAtlas2 k where x.kgID = k.name and spDisplayID is not null order by 1 " |tawk '$2!=""{print}' > gnfGeneName2.tab

echo "SCORE Threshold is $SCORETHRESH"
tawk '$5>'$SCORETHRESH' {print $0}' ../$TABLE.bed | overlapSelect ../pseudoEstAll.bed stdin retroExpress.bed -inFmt=bed
#retroArray has retros with expression data
overlapSelect array.uniq.bed retroExpress.bed retroArray.bed 
overlapSelect array.uniq.bed ../$TABLE.bed retroAllArray.bed 
#retroArrayKnown has retros overlapping knownGenes with expression data
overlapSelect ../knownGene.tab.gz retroArray.bed retroArrayKnown.bed -selectFmt=genePred
overlapSelect ../knownGene.tab.gz array.uniq.bed geneTissue.bed -selectFmt=genePred
overlapSelect ../knownGene.tab.gz retroArrayKnown.bed  -selectFmt=genePred -idOutput retroKg.id
overlapSelect retroArrayKnown.bed array.uniq.bed retroTissue.bed

wc -l retroArray.bed retroArrayKnown.bed geneTissue.bed

#get experiment id and name
hgsql hgFixed -N -B -e "select name from $ARRAYEXP" |sed -e "s/ /_/g"|sed -e "s/+//g" |sed -e "s/(/_/g"|sed -e "s/)//g" > expName.tab
hgsql hgFixed -N -B -e "select id+2, name from $ARRAYEXP " |sed -e "s/ /_/g"|sed -e "s/+//g" |sed -e "s/(/_/g"|sed -e "s/)//g" > expName.tab
exit
#run edge 
awk -F\' '{OFS="\t";$1=$1;print $2,$3,$4,$1}' fullsetRetro.txt |tawk '$3<0.01{print $1}' > 1percentRetroFull.id
awk -F\' '{OFS="\t";$1=$1;print $2,$3,$4,$1}' fullsetRetro.txt |tawk '$3<0.05{print $1}' > 5percentRetroFull.id

$SCRIPT/selectById 1 5percentRetroFull.id 10 affy.psl > 5percentRetroFull.psl
$SCRIPT/selectById 1 1percentRetroFull.id 10 affy.psl > 1percentRetroFull.psl

overlapSelect 5percentRetroFull.psl ../$TABLE.bed 5percentRetroFull.bed
overlapSelect 1percentRetroFull.psl ../$TABLE.bed 1percentRetroFull.bed

#prepare set for clustering in genesis
$SCRIPT/selectById 1 1percentRetroFull.id 1 $DB.LogRatioAll.tab > $DB.LogRatioSig.txt



#
# Tau stat - tissue specificity per probe using absolute expression level, output to gnfAtlas2_abs_median.tab 
#
hgTissueSpecific $DB $ARRAYMEDIAN $ARRAYEXP gnfAtlas2_abs_median -lookup=$ARRAYLOOKUP 
#get parent kg id and sort by it
tawk '{print $47,$4}' retroArrayKnown.bed | sort > retroArrayParentKnown.sort.bed
# retroArrayParentKnown.sort.bed contains parent kgid 
sort gnfAtlas2_abs_median.tab > gnfAtlas2_abs_median.sort.tab
#custom column track for gene sorter
cut -f 2  gnfAtlas2_abs_median.sort.tab > tauAbsolute.tab
cp tauAbsolute $ROOTDIR/retro/array

#retroParentTissue adds tau stat to parent of expressed retros
join retroArrayParentKnown.sort.bed gnfAtlas2_abs_median.sort.tab |awk '{OFS="\t";$1=$1;print $2,$1,$3,$4}' | sort > retroParentTissue.tab
tawk '{print $2,$1}' retroKg.id |sort> retroKg.sort.id
join retroKg.sort.id gnfAtlas2_abs_median.sort.tab |awk '{OFS="\t";$1=$1;print $2,$1,$3,$4}' | sort > retroChildTissue.tab
#
# extract tissue specific expression (log ratio)
#
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
cut -f 1-6,15 geneTissue.bed |tawk '{print $1,$2,$3,$4,$5,$6,","$7}'> geneTissueExtract.bed

# id is col offset into tissue table
mkdir -p tissue
tawk '{print "awk -F, xx{print $1, tt\$"$1"}xx retroTissueExtract.bed |tawk xx{\$8\=\(\$8\+5)\*10;print $1,$2,$3,$4,$5,$6,$8}xx > tissue/retroTissue."$2".bed"} ' expName.tab |sed -e "s/xx/'/g" |sed -e 's/tt/"\t"/' > buildTissues.sh
tawk '{print "awk -F, xx{print $1, tt\$"$1"}xx geneTissueExtract.bed |tawk xx{\$8\=\(\$8\+5)\*10;print $1,$2,$3,$4,$5,$6,$8}xx > tissue/geneTissue."$2".bed"} ' expName.tab |sed -e "s/xx/'/g" |sed -e 's/tt/"\t"/' > buildGeneTissues.sh
chmod +x buildTissues.sh buildGeneTissues.sh
./buildTissues.sh
./buildGeneTissues.sh
mkdir -p retro.tissue
mkdir -p gene.tissue
for i in `ls tissue/retro*.bed`; do cut -f 1-3 $i | overlapSelect -selectFmt=bed stdin ../$TABLE.bed retro.$i ; done
for i in `ls tissue/gene*.bed`; do cut -f 1-3 $i | overlapSelect -selectFmt=bed stdin ../$TABLE.bed gene.$i ; done

for t in `cut -f 2 expName.tab` ; do $SCRIPT/makeHtmlRightArray.sh DEF retro.tissue/retroTissue.$t.bed array/tissue/$DB/$t all_$t;done

mkdir -p $ROOTDIR/retro/array/$DB
F=`awk -F, '{print NF}' geneTissueExtract.bed | uniq`  ; let F=F-1
echo "number of experiments plus one = $F"
cut -d , -f 2-$F geneTissueExtract.bed > geneTissue.$DB.tab
cut -d , -f 2-$F retroTissueExtract.bed > retroTissue.$DB.tab
#
#prepare for sam run
#
echo 'source("'$SCRIPT'/overHist.R")'> plotRetroVsAll.R
echo 'h<-scan("retroTissue.'$DB'.tab",sep=",")' >> plotRetroVsAll.R
echo "dim(h)<-c(79,length(h)/79)" >> plotRetroVsAll.R
echo 'm<-scan("geneTissue.'$DB'.tab",sep=",")' >> plotRetroVsAll.R
echo "dim(m)<-c(79,length(m)/79)" >> plotRetroVsAll.R
tawk '{print "pdf(\"'$ROOTDIR'/retro/array/'$DB'/gene."$2".'$DB'.pdf\")";print "overhist(h["$1"-1,],m["$1"-1,], breaks=100, xlab0=\""$2" (blue=Retros green=All Genes)\")";print "dev.off()"}' expName.tab  >> plotRetroVsAll.R
#R --slave < plotRetroVsAll.R 

#for i in `ls retroTissue.*.bed` ; do textHistogram $i -col=7 -autoScale=30 ; done > tissue.hist
wc -l retroTissueExtract.bed 
echo "run arrayCompare.sh for R analysis after all species are complete"

# post process sam data
#awk '{print $4}' $DB.siggenes.txt |sed -e 's/"//g' > $DB.siggenes.id
#$SCRIPT/selectById 1 $DB.siggenes.id 10 affy.psl > $DB.RetroSam.psl
#overlapSelect $DB.RetroSam.psl ../$TABLE.bed $DB.RetroSam.bed


