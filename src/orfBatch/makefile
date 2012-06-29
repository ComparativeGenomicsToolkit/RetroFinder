progs = orfBatch
include ../../common.mk

DB = hg18
test:	test1 test2 test3

test%:	input/test%.bed
	@mkdir -p output
	${bindir}/orfBatch $(DB) input/$@.bed output/$@p.bed output/$@.gp > output/$@b.tab 
	diff expected/$@p.bed output/$@p.bed 
	diff expected/$@p.bed output/$@p.bed | awk '{print $$15":"$$17"-"$$18}'
	diff expected/$@b.tab output/$@b.tab
	diff expected/$@.gp output/$@.gp
