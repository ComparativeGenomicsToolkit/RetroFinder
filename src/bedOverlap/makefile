progs = bedOverlap
include ../../common.mk


test:	test2 test4 test5
%:	input/%.bed
	@mkdir -p output
	${bindir}/bedOverlap -hasBin -verbose=1 input/$@.bed output/$@.bed  
	diff expected/$@.bed output/$@.bed 

%:	input/%.noBin
	@mkdir -p output
	${bindir}/bedOverlap input/$@.noBin output/$@.bed 
	diff expected/$@.bed output/$@.bed 


clean::
	rm -rf output
