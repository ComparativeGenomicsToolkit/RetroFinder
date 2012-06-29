progs = 
include ../../common.mk

test: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test13 test14 test15 test16 test17 test18 test19 testRp131 test26 test38 test1987

clean::
	rm -fr ${O} output

t%:	${INPUT}/t%.psl
	${mkdir} output
	${bindir}/pslNetOverlap $(OPTIONS) ${INPUT}/$@.psl $(DBDIR)/chrom.sizes $(PSEUDO)/rmsk.bed.gz $(PSEUDO)/netMm8.txt.gz $(PSEUDO)/netCanFam2.txt.gz $(PSEUDO)/netRheMac2.txt.gz $(PSEUDO)/netPanTro2.txt  > output/$@l.txt 2> $@.log
	diff expected/$@.bed output/$@.bed 
