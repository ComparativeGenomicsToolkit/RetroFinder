include ../../../common.mk

PSLPSEUDO = ../../../bin/pslPseudo 
SCRIPT = ../../pipeline
DIFF = diff -u --exclude CVS -r
DB = panTro4
DBDIR = inputData/$(DB)
BED = $(DBDIR)/bed
PSEUDO = $(DBDIR)/bed/pseudo
MRNA = $(DBDIR)/mrna.2bit
BUILD=gs.19/build36
inSizes = $(DBDIR)/chrom.sizes 
RMSK = $(DBDIR)/rmsk.bed.gz 
NET1 = $(DBDIR)/netMm10.test.txt.gz 
NET2 = $(DBDIR)/netCanFam3.test.txt.gz 
TRF = $(DBDIR)/simpleRepeat.bed.gz 
MRNA_ALIGN = $(DBDIR)/all_mrna.psl.gz 
INFILES = $(inSizes) $(RMSK) $(NET1) $(NET2) $(TRF) $(MRNA_ALIGN)
OPTIONS = -minAli=0.97 -nearTop=0.005 -verbose=5 -cdsFile=$(DBDIR)/cds.tab.gz
NIBLIST = $(DBDIR)/S1.lst
GENES.ALL = $(DBDIR)/refGene.tab.gz $(DBDIR)/mgcGene.tab.gz $(DBDIR)/sortedKnownGene.tab.gz
GENES = $(DBDIR)/refGene.tab.gz $(DBDIR)/ensGene.tab.gz $(DBDIR)/refGene.tab.gz
NET3 =  $(DBDIR)/netHg19.test.txt.gz
MKDIR = mkdir -p

test:	hg18_test11 hg18_test12 hg18_test13 panTro4_test996 
	${MAKE} doPseudo name=$@ inPsl=$@.psl filtArgs='' inSizes=chrom.sizes

# nothing should be filtered with no arguments
hg17%:
	${MAKE} DB=hg17 doPseudo name=$@ inPsl=$@.psl filtArgs='' GENES="inputData/hg17/knownGene.tab.gz inputData/hg17/refGene.tab.gz inputData/hg17/mgcGene.tab.gz" NET1=inputData/hg17/netMm8.tab.gz NET2=inputData/hg17/netCanFam2.tab.gz NET3=inputData/hg17/netRheMac2.tab.gz 

hg18%:
	${MAKE} DB=hg18 PSLPSEUDO=/cluster/home/baertsch/bin/x86_64/pslPseudo oldPseudo name=$@ inPsl=$@.psl filtArgs='-skipBlatMerge -verbose=5 -minAli=0.98 -nearTop=0.005' GENES="inputData/hg18/sortedKnownGene.tab.gz inputData/hg18/refGene.tab.gz inputData/hg18/mgcGenes.tab.gz" NET1=inputData/hg18/netMm8.txt.gz NET2=inputData/hg18/netCanFam2.txt.gz NET3=inputData/hg18/netRheMac2.txt.gz 

panTro4%:
	${MAKE} DB=panTro4 doPseudo name=$@ inPsl=$@.psl filtArgs='-stripVersion -verbose=5 -minAli=0.98 -nearTop=0.005' GENES="inputData/panTro4/ensGene.tab.gz inputData/panTro4/refGene.tab.gz inputData/panTro4/refGene.tab.gz" NET1=inputData/panTro4/netMm10.test.txt.gz NET2=inputData/panTro4/netCanFam3.test.txt.gz NET3=inputData/panTro4/netHg19.test.txt.gz 

noopTest:
	${MAKE} doPseudo name=$@ inPsl=overlap.psl filtArgs=''

est5Mrna100:
	${MAKE} doPseudo name=$@ inPsl=pseudoEst5MrnaOrf100.psl filtArgs='' 

NM_%:
	${MAKE} doPseudo name=$@ inPsl=$@.psl filtArgs='' inSizes=chrom.sizes

B%:
	${MAKE} doPseudo name=$@ inPsl=$@.psl filtArgs='' inSizes=chrom.sizes

overlapTest:
	${MAKE} doPseudo name=$@ inPsl=overlap.psl inSizes=overlap.sizes \
	    filtArgs='-bestOverlap -minCover=0.15 -minId=0.96'

maxAlignTest:
	${MAKE} doPseudo name=$@ inPsl=NM_004038.3.psl inSizes=NM_004038.3.sizes \
	    filtArgs='-maxAligns=2'

# Recursive targets.  Results are sorted to allow consistent comparisons, as
# sometimes internal sorting done by pslPseudo is not stable.
#  o name - test name
#  o filtArgs - filter arguments
#  o inPsl - input psl, relative to input dir
#  o inSize - polyA sizes file, relative to input dir (optional)
#$DB $CACHE/split/tmp$1.psl $ICACHE/chrom.sizes $ICACHE/rmsk.bed.gz $ICACHE/mouseNet.txt.gz $ICACHE/dogNet.txt.gz $ICACHE/simpleRepeat.bed.gz $ICACHE/all_mrna.psl.gz $OUT/mrna$1.psl $OUT/pseudo$1.psl /tmp/pseudoMrnaLink$1.txt /tmp/pseudo$1.axt $ICACHE/S1.lst $CACHE/mrna.2bit $ICACHE/refGene.tab.gz $ICACHE/mgcGene.tab.gz $ICACHE/sortedKnownGene.tab.gz $ICACHE/rheMac2Net.txt.gz > /tmp/pseudo$1.log 2> /tmp/err$1.log
#ifneq (${inSizes},)
#	sizesOpt=-polyASizes=input/${inSizes}
#endif
outDir=output/${name}
expectDir=expected
sortPsl=sort -k 10,10 -k 12,12n -k 13,13n
doPseudo:
	@${MKDIR} ${outDir}
	${PSLPSEUDO} $(OPTIONS) ${filtArgs} $(DB) input/${inPsl} $(INFILES) ${outDir}/${name}g.psl ${outDir}/keep.psl.tmp ${outDir}/${name}.raw.bed ${outDir}/${name}.axt ${NIBLIST} $(MRNA) $(GENES) $(NET3) ${outDir}/ortho${name}.tab > ${name}.log
	tawk -f $(SCRIPT)/filterBed.awk ${outDir}/${name}.raw.bed | sort -k1,1 -k2,3n -k4,4 >${outDir}/${name}.sort.bed
	${SCRIPT}/pslPseudoPost ${outDir}/${name}.sort.bed ${outDir}/${name}.bed
	diff ${expectDir}/${name}.bed ${outDir}/${name}.bed 
oldPseudo:
	@${MKDIR} ${outDir}
	${PSLPSEUDO} $(OPTIONS) ${filtArgs} $(DB) input/${inPsl} $(INFILES) ${outDir}/${name}g.psl ${outDir}/keep.psl.tmp ${outDir}/${name}.raw.bed ${outDir}/${name}.axt ${NIBLIST} $(MRNA) $(GENES) $(NET3) 
	tawk -f $(SCRIPT)/filterBed.old.awk ${outDir}/${name}.raw.bed | sort -k1,1 -k2,3n -k4,4 >${outDir}/${name}.sort.bed
	${SCRIPT}/pslPseudoPost ${outDir}/${name}.sort.bed ${outDir}/${name}.bed
	diff ${expectDir}/${name}.bed ${outDir}/${name}.bed 
doSave:
	@${MKDIR} ${outDir}
	${PSLPSEUDO} $(OPTIONS) ${filtArgs} $(DB) input/${inPsl} $(INFILES) ${outDir}/${name}g.psl ${outDir}/keep.psl.tmp ${outDir}/${name}.raw.bed /dev/null ${NIBLIST} $(MRNA) $(GENES) >${outDir}/${name}.log 2>&1
	tawk -f $(SCRIPT)/filterBed.awk ${outDir}/${name}.raw.bed | sort -k1,1 -k2,3n -k4,4 > ${outDir}/${name}.sort.bed
	bedOverlap -noBin ${outDir}/${name}.sort.bed ${outDir}/${name}.bed
	tawk '{print $1,$2,$3,$4,$5,$6}' ${outDir}/${name}.bed | overlapSelect expected/{name}.est5.bed stdin ${outDir}/${name}.est5.bed -inFmt=bed
	@${sortPsl} ${outDir}/keep.psl.tmp > ${outDir}/keep.psl
	@rm -f ${outDir}/keep.psl.tmp
	${DIFF} expected/${name}.bed output/${name}.bed
	${DIFF} expected/${name}.est5.bed output/${name}.est5.bed
	${DIFF} expected/${name}.psl output/${name}.psl

# different parameters
doPseudoAlt:
	@${MKDIR} ${outDir}
	${PSLPSEUDO} ${filtArgs} -verbose=1 -dropped=${outDir}/drop.psl.tmp -weirdOverlapped=${outDir}/weird.psl.tmp input/${inPsl} ${outDir}/keep.psl.tmp >${outDir}/filt.out 2>&1
	@${sortPsl} ${outDir}/keep.psl.tmp > ${outDir}/keep.psl
	@${sortPsl} ${outDir}/drop.psl.tmp > ${outDir}/drop.psl 
	@${sortPsl} ${outDir}/weird.psl.tmp > ${outDir}/weird.psl
	@rm -f ${outDir}/keep.psl.tmp ${outDir}/drop.psl.tmp ${outDir}/weird.psl.tmp
	${DIFF} expected/${name} output/${name}
