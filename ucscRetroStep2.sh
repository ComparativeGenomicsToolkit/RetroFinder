#!/bin/bash 
set -beEu -o pipefail
source $1
#align human mrnas to $DB using lastz
cd $TMPMRNA
echo "working directory is $TMPMRNA."
#concatenate, sort (by name, score) and de-dup psls by chrom
mkdir -p pslFilter
for i in `awk '{print $1}' S1.len` ; do echo $i ; cat lastz/$i/*.psl | awk '{print $0, $1*3-$2}' | \
 sort -k 10,10 -k 22nr -T /scratch | awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | \
 /cluster/home/baertsch/bin/x86_64/pslFilterDups stdin pslFilter/$i.psl  ; done 

#chain blocks

ls $NIB/* > nib.lst
mkdir -p psl
for i in `awk '{print $1}' S1.len` ; do nohup $TMPMRNA/doChain $i ; done 

#convert chains to psl
#for i in `awk '{print $1}' S1.len`; do chainToPsl chainFilter/$i.chain S1.len S2.len nib.lst trim.fa psl/$i.psl; done


#reattach polyA tail and fix alignments
mkdir -p pslLift
for i in `awk '{print $1}' S1.len` ; do liftUp pslLift/$i.psl mrna.lft warn psl/$i.psl -pslQ -nohead ; done

cd $MRNABASE
pslCat -nohead $TMPMRNA/pslLift/*psl > mrnaBlastz.psl

sort -k10,10 -k14,14 -k16,16n -k12,12n mrnaBlastz.psl > mrnaBlastz.sort.psl
#pslCDnaFilter -minCover=0.05 -minId=0.65 mrnaBlastz.sort.psl mrnaBlastz.65.psl
##                       seqs    aligns
#             total:     223377  11187645
#     drop minIdent:     103075  6082119
#     drop minCover:     51568   457653
#              kept:     222565  4647873
#pslCDnaFilter -minCover=0.05 -minId=0.62 mrnaBlastz.sort.psl mrnaBlastz.psl
#                        seqs    aligns
#             total:     223377  11187645
#     drop minIdent:     74887   3622214
#     drop minCover:     52526   567156
#              kept:     223112  6998275

pslCDnaFilter -minCover=0.05 -minId=0.58 mrnaBlastz.sort.psl mrnaBlastz.psl
rm -f mrnaBlastz.sort.psl.gz
gzip mrnaBlastz.sort.psl &

#load mrna alignment track into browser
awk -f $SCRIPT/stripversion.awk mrnaBlastz.psl | hgLoadPsl $DB stdin -table=mrnaBlastz

#split for pipeline cluster run 
mkdir -p $OUTDIR/split
/cluster/home/baertsch/bin/x86_64/pslSplit nohead $OUTDIR/split mrnaBlastz.psl -chunkSize=120

cd $OUTDIR/split
for i in `ls tmp*.psl` ; do $SCRIPT/pslQueryUniq $i > temp.psl ; mv temp.psl $i ;echo $i; done
grep chr tmp* | awk '{print $1,$10}' | awk -F":" '{print $1,$2}'|awk '{print $1,$3}'|uniq  > acc.lst

#load mrna sequences into browser (with version numbers)
mkdir -p /gbdb/$DB/blastzRetro
rm -f /gbdb/$DB/blastzRetro/mrna.fa 
rm -f /gbdb/$DB/blastzRetro/refseq.fa 
ln $MRNABASE/mrna.fa /gbdb/$DB/blastzRetro/ -s
ln $MRNABASE/refseq.fa /gbdb/$DB/blastzRetro/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/refseq.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/mrna.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile

echo "ucscRetroStep2.sh mrna alignments complete"
cd $TMPMRNA
cp DEF $OUTDIR
cp S1.len $OUTDIR
cd $OUTDIR
pwd
echo "run ucscRetroStep3.sh DEF to run retro pipeline"
