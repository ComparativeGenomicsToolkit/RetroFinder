progs = pslPseudo
include ../../common.mk

# FIXME: not working, include not checked in
DB = hg18
DBDIR = /cluster/data/$(DB)
BED = $(DBDIR)/bed
PSEUDO = /hive/users/baertsch/retro/$(DB)/
MRNA = $(PSEUDO)/mrna.2bit 
OPTIONS = -minAli=0.98 -nearTop=0.005 -verbose=7 -skipBlatMerge -cdsFile=$(PSEUDO)/cds.tab.gz
VALCMD = chrom.sizes ${PSEUDO}/rmsk.bed $(PSEUDO)/mouseNet.bed $(PSEUDO)/chr13simple.bed $(PSEUDO)/chr13_mrna.psl 
TESTCMD = $(DBDIR)/chrom.sizes ${PSEUDO}/rmsk.bed.gz $(PSEUDO)/netMm8.txt.gz $(PSEUDO)/netCanFam2.txt.gz $(BED)/simpleRepeat/simpleRepeat.bed $(PSEUDO)/all_mrna.psl.gz 
ICACHE= /hive/users/baertsch/retro/$(DB)/
#TESTCMD = $(DBDIR)/chrom.sizes mtRmsk mtSyn.bed mtSimple.bed input/test48.psl
GENES = $(PSEUDO)/refGene.tab.gz $(PSEUDO)/mgcGenes.tab.gz $(PSEUDO)/knownGene.tab.gz


test::
	cd tests ; make test
#test: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test13 test14 test15 test16 test17 test18 test19 testRp131 test26 test38 test1987
t%:	input/t%.psl
	@mkdir -p output
	${bindir}/pslPseudo $(OPTIONS) $(DB) input/$@.psl $(TESTCMD)  output/$@g.psl output/$@p.psl output/$@l.txt output/$@.axt S1.$(DB).lst $(MRNA) $(GENES) $(PSEUDO)/netRheMac2.txt.gz > $@.log 
	../pipeline/pslPseudoPost output/$@l.txt output/$@.bed
	diff expected/$@.bed output/$@.bed 

clean::
	rm -fr ${O} output
