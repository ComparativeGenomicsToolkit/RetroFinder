progs = countFrameShift
include ../../common.mk


test: test.5930  test.Y00345

test.%:	input/%.psl
	@mkdir -p output
	countFrameShift -verbose=3 input/$*.psl -cdsDb=hg18 input/refGene.gp output/$@.gp 2> output/$@.log 
