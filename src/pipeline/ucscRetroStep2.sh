#!/bin/bash 
set -bevEu -o pipefail
source $1
#align human mrnas to $DB using lastz
# cd $TMPMRNA
echo "working directory is $TMPMRNA."
wc -l $TMPMRNA/run.0/jobList | awk '{print $1}'> $TMPMRNA/jobs.cnt
# jobs.lst will then contain full path
find $TMPMRNA/lastz -name \*.psl| wc -l | awk '{print $1}'> $TMPMRNA/jobs.lst
echo Check Job Count from cluster run
diff $TMPMRNA/jobs.cnt $TMPMRNA/jobs.lst 
if [ $? != 0 ]; then
  echo missing jobs aborting
  wc -l $TMPMRNA/jobs.cnt $TMPMRNA/jobs.lst
  exit 3
fi
#concatenate, sort (by name, score) and de-dup psls by chrom
mkdir -p $TMPMRNA/pslFilter
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do echo $i ; cat $TMPMRNA/lastz/$i/*.psl | awk '{print $0, $1*3-$2}' | \
 sort -k 10,10 -k 22nr -T /scratch | awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | \
 ${BINDIR}/pslFilterDups stdin $TMPMRNA/pslFilter/$i.psl  ; done 

#chain blocks

mkdir -p $TMPMRNA/psl
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do nohup $TMPMRNA/doChain $i ; done 

#reattach polyA tail and fix alignments
mkdir -p $TMPMRNA/pslLift
for i in `awk '{print $1}' $TMPMRNA/S1.len` ; do liftUp $TMPMRNA/pslLift/$i.psl $TMPMRNA/mrna.lft warn $TMPMRNA/psl/$i.psl -pslQ -nohead ; done

# cd $TMPMRNA/pslLift
echo "MRNABASE $MRNABASE"
pslCat -nohead $TMPMRNA/plsLift/*psl > $MRNABASE/mrnaBlastz.psl

# cd $MRNABASE
sort -k10,10 -k14,14 -k16,16n -k12,12n $MRNABASE/mrnaBlastz.psl > $MRNABASE/mrnaBlastz.sort.psl
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

pslCDnaFilter -minCover=0.05 -minId=0.58 $MRNABASE/mrnaBlastz.sort.psl $MRNABASE/mrnaBlastz.psl
rm -f $MRNABASE/mrnaBlastz.sort.psl.gz
gzip $MRNABASE/mrnaBlastz.sort.psl &

#split for pipeline cluster run 
mkdir -p $OUTDIR/split
${BINDIR}/pslSplit nohead $OUTDIR/split $MRNABASE/mrnaBlastz.psl -chunkSize=120

# cd $OUTDIR/split
for i in `ls $OUTDIR/split/tmp*.psl` ; do $SCRIPT/pslQueryUniq $i > $OUTDIR/split/temp.psl ; mv $OUTDIR/split/temp.psl $i ; done
grep chr $OUTDIR/split/tmp* | awk '{print $1,$10}' | awk -F":" '{print $1,$2}'|awk '{print $1,$3}'|uniq  > $OUTDIR/split/acc.lst

export TMPDIR=/scratch/tmp
echo TMPDIR
#load mrna alignment track into browser
awk -f $SCRIPT/stripversion.awk $MRNABASE/mrnaBlastz.psl | hgLoadPsl $DB stdin -table=mrnaBlastz

echo "ucscRetroStep2.sh mrna alignments complete"
# cd $MRNABASE
cp $MRNABASE/S1.len $OUTDIR
# cd $OUTDIR
# pwd

#load mrna sequences into browser (with version numbers)
mkdir -p /gbdb/$DB/blastzRetro${VERSION}
rm -f /gbdb/$DB/blastzRetro${VERSION}/mrna.fa 
ln $MRNABASE/mrna.fa /gbdb/$DB/blastzRetro${VERSION}/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro${VERSION}/mrna.fa  -seqTbl=ucscRetroSeq${VERSION} -extFileTbl=ucscRetroExtFile${VERSION}
rm -f /gbdb/$DB/blastzRetro${VERSION}/refseq.fa 
ln $MRNABASE/refseq.fa /gbdb/$DB/blastzRetro${VERSION}/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro${VERSION}/refseq.fa  -seqTbl=ucscRetroSeq${VERSION} -extFileTbl=ucscRetroExtFile${VERSION}

echo "run ucscRetroStep3.sh DEF to run retro pipeline"
