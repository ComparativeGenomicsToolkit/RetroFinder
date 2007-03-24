#
tawk '{print $1,$2,$3,$4,$5,$6}' ../$1 | overlapSelect -selectFmt=bed stdin ../retroMrnaAli.psl stdout > retroMrnaAli.psl 
#overlapSelect refGeneCds.gp retroMrnaAli.psl retroMrnaCds.out -statsOutput -overlapThreshold=0.7
#overlapSelect refGeneCds.gp retroMrnaAli.psl retroMrnaCds.psl -overlapThreshold=0.7
overlapSelect bothCds.gp retroMrnaAli.psl retroMrnaCds.out -statsOutput #-overlapSimilarity=0.7
overlapSelect bothCds.gp retroMrnaAli.psl retroMrnaCds.psl #-overlapSimilarity=0.7           

#awk -f stripversion.awk retroMrnaCds.psl  > retroMrnaAli.strip.psl
#mrnaToGene retroMrnaAli.strip.psl -cdsDb=hg18 -cdsMergeMod3 -cdsMergeSize=100 -utrMergeSize=20 -genePredExt retroParent.gp  2> retroPred.err
mrnaToGene retroMrnaAli.psl -ignoreExtension -cdsDb=hg18 -cdsMergeMod3 -cdsMergeSize=100 -utrMergeSize=20 -genePredExt retroParent.gp  2> retroPred.err > retroPred.log
ldHgGene hg18 rbRetroParent retroParent.gp -predTab -genePredExt
countFrameShift -verbose=3 -cdsDb=hg18 retroMrnaAli.psl both.gp tmp 2> err
sort -k5,5nr < tmp > type1.$1 
#sort -k5,5nr type1.$1.bed >type1.$1.sort.bed
#~markd/bin/bedToHtmlDir -dir-frame-per 25 -context-bases 500 -browser-url http://hgwdev-baertsch.cse.ucsc.edu hg18 type1.sort.bed /usr/local/apache/htdocs/retro/type1
cut -f1-6 type1.$1 > type1.short.bed
overlapSelect type1.short.bed retroMrnaAli.psl missing.psl -nonOverlapping
wc -l missing.psl

makeHtml.sh hg18 type1.$1 type1/$2 $2
#overlapSelect type1.short.bed pnasretroExp.from.hg17.bed  notKaess.bed -nonOverlapping

