#align human mrnas to hg18 using lastz
DB=hg18

cd /san/sanvol1/scratch/mrnaBlastz/$DB/

#concatenate, sort (by name, score) and de-dup psls by chrom

mkdir -p pslFilter
for i in `awk '{print $1}' S1.len` ; do echo $i ; cat lastz/$i/*.psl | awk '{print $0, $1*3-$2}' | \
 sort -k 10,10 -k 22nr -T /scratch | awk '{OFS=" "; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' | \
 pslFilterDups stdin pslFilter/$i.psl  ; done 

#chain blocks

mkdir -p psl
for i in `awk '{print $1}' S1.len` ; do nohup doChain $i ; done 

#convert chains to psl
ls /scratch/data/$DB/nib/*.nib > nib.lst
#for i in `awk '{print $1}' S1.len`; do chainToPsl chainFilter/$i.chain S1.len S2.len nib.lst trim.fa psl/$i.psl; done


#reattach polyA tail and fix alignments
mkdir -p pslLift
for i in `awk '{print $1}' S1.len` ; do liftUp pslLift/$i.psl mrna.lft warn psl/$i.psl -pslQ -nohead ; done

cd /cluster/data/$DB/bed/mrnaBlastz
pslCat -nohead /san/sanvol1/scratch/mrnaBlastz/$DB/pslLift/*psl > mrnaBlastz.psl

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

pslCDnaFilter -minCover=0.10 -minId=0.62 mrnaBlastz.sort.psl mrnaBlastz.psl

#load mrna alignment track into browser
awk -f ~baertsch/bin/scripts/stripversion.awk mrnaBlastz.psl | hgLoadPsl $DB stdin -table=mrnaBlastz2
#rm -f mrnaBlastz.strip.psl 

#grep -v random mrnaBlastz.psl > mrnaBlastz.norand.psl

#awk '{print $10,$14,$9}' mrnaBlastz.psl| sort|  uniq -c > qName.count &
#nice awk '$1>1{print $0}' qName.count > qName.countDups 
#select.awk qName.countDups mrnaBlastz.psl > mrnaBlastz.dups.psl 
#awk -f missingChain.awk mrnaBlastz.dups.psl > mrnaBlastz.missing.psl

#split for pipeline cluster run 
mkdir -p /hive/users/baertsch/retro/$DB/split
pslSplit nohead /hive/users/baertsch/retro/$DB/split mrnaBlastz.psl -chunkSize=120
cd /hive/users/baertsch/retro/$DB/split
for i in `ls tmp*.psl` ; do ~baertsch/bin/scripts/pslQueryUniq $i > temp.psl ; mv temp.psl $i ;echo $i; done
grep chr tmp* | awk '{print $1,$10}' | awk -F":" '{print $1,$2}'|awk '{print $1,$3}'|uniq  > acc.lst

#load mrna sequences into browser (with version numbers)
#cp mrna.fa /hive/data/genomes/$DB/bed/mrnaBlastz/mrna.fa
mkdir -p /gbdb/$DB/blastzRetro
rm -f /gbdb/hg18/blastzRetro/mrna.fa 
rm -f /gbdb/hg18/blastzRetro/refseq.fa 
ln /hive/data/genomes/$DB/bed/mrnaBlastz/mrna.fa /gbdb/$DB/blastzRetro/ -s
ln /hive/data/genomes/$DB/bed/mrnaBlastz/refseq.fa /gbdb/$DB/blastzRetro/ -s
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/refseq.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile
hgLoadSeq -replace $DB /gbdb/$DB/blastzRetro/mrna.fa  -seqTbl=ucscRetroSeq -extFileTbl=ucscRetroExtFile
