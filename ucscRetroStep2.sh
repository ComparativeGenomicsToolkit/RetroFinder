#!/bin/bash 
set -beEu -o pipefail
source $1
#align human mrnas to $DB using lastz
cd $TMPMRNA
echo "working directory is $TMPMRNA."
wc -l run.0/jobList | awk '{print $1}'> jobs.cnt
find lastz -name \*.psl| wc -l | awk '{print $1}'> jobs.lst
echo Check Job Count from cluster run
diff jobs.cnt jobs.lst 
if [ $? != 0 ]; then
  echo missing jobs aborting
  wc -l jobs.cnt jobs.lst
  exit 3
fi
#concatenate, sort (by name, score) and de-dup psls by chrom
mkdir -p pslFilter
for i in `awk '{print $1}' S1.len` ; do echo $i ; cat lastz/$i/*.psl | awk '{print $0, $1*3-$2}' | \
 sort -k 10,10 -k 22nr -T /scratch | awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | \
 /cluster/home/baertsch/bin/x86_64/pslFilterDups stdin pslFilter/$i.psl  ; done 

#chain blocks

find $NIB > nib.lst
mkdir -p psl
for i in `awk '{print $1}' S1.len` ; do nohup $TMPMRNA/doChain $i ; done 

#reattach polyA tail and fix alignments
mkdir -p pslLift
for i in `awk '{print $1}' S1.len` ; do liftUp pslLift/$i.psl mrna.lft warn psl/$i.psl -pslQ -nohead ; done

cd $TMPMRNA/pslLift
pslCat -nohead *psl > $MRNABASE/mrnaBlastz.psl

cd $MRNABASE
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

#split for pipeline cluster run 
mkdir -p $OUTDIR/split
/cluster/home/baertsch/bin/x86_64/pslSplit nohead $OUTDIR/split mrnaBlastz.psl -chunkSize=120

cd $OUTDIR/split
for i in `ls tmp*.psl` ; do $SCRIPT/pslQueryUniq $i > temp.psl ; mv temp.psl $i ;echo $i; done
grep chr tmp* | awk '{print $1,$10}' | awk -F":" '{print $1,$2}'|awk '{print $1,$3}'|uniq  > acc.lst

#load mrna alignment track into browser
awk -f $SCRIPT/stripversion.awk $MRNABASE/mrnaBlastz.psl | hgLoadPsl $DB stdin -table=mrnaBlastz

echo "ucscRetroStep2.sh mrna alignments complete"
cd $MRNABASE
cp S1.len $OUTDIR
cd $OUTDIR
pwd

#load mrna sequences into browser (with version numbers)
mkdir -p /gbdb/$DB/blastzRetro
rm -f /gbdb/$DB/blastzRetro/mrna.fa 
ln $MRNABASE/mrna.fa /gbdb/$DB/blastzRetro/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/mrna.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile
rm -f /gbdb/$DB/blastzRetro/refseq.fa 
ln $MRNABASE/refseq.fa /gbdb/$DB/blastzRetro/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/refseq.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile

echo "run ucscRetroStep3.sh DEF to run retro pipeline"
