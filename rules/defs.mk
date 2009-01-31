
# must define:
#  db=mm9
#  ncbiBuild=36.3
#  retroDate=20080825
#  set=2008-10-30
#  synOrgs = bosTau2 canFam2 hg18 monDom4 panTro2 rheMac2 rn4
#  GBID = hg
#  gbRoot=/cluster/data/genbank
#  splitSize = 1000000 

ifeq (${chroms},)
   # small first
   chroms = $(shell sort -k 2,2n /cluster/data/${db}/chrom.sizes | tawk '!(/hap/||/random/||/chrM/){print $$1}')
endif

#setDir=sets/${set}
tmpDir=/san/sanvol1/scratch/mrnaBlastz/${db}
wwwDir=www/${db}/${set}

retroDb=retro

#retroData=${setDir}/${db}
genomeDir=/scratch/data/${db}/nib
mrnaAlignDir=${tmpDir}/align
retroGp = ${retroData}/retro.${db}.gp.gz
retroInfo = ${retroData}/retro.${db}.info.gz
retroIds = ${retroData}/retro.${db}.ids
mrnaSeq = mrna.fa
refseq = refseq.fa
trimSeq = trim.fa.gz
trimLen = trim.len
mrnaLen = mrna.len
mrna2bit = mrna.2bit
GBID = hg
splitSize = 1000000 

retroBed = ${hgData}/retroMrna.bed.gz

retroBin=/cluster/home/baertsch/bin/x86_64
scripts=/cluster/home/baertsch/bin/scripts
gbBin=/cluster/data/genbank/bin/x86_64
gbRoot=/cluster/data/genbank

#tmpExt = $(shell echo $$PPID).${HOST}.tmp
tmpExt = $(shell echo "1728").${HOST}.tmp

export PATH:=:${ROOT}/bin:${PATH}

paraHost = pk
#paraHost = kk
ifeq (${paraHost},pk)
   clusterRootDir=/san/sanvol1/scratch/baertsch
else
   clusterRootDir=/hive/scratch/baertsch
endif
clusterDir=${clusterRootDir}/retro/builds/${db}

wwwDir = /cluster/home/baertsch/public_html/retro/${db}
wwwDirUrl  = http://hgwdev.cse.ucsc.edu/~baertsch/retro/${db}

.SECONDARY:

