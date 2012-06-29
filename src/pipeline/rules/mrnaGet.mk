##
##

mrnaGet: ${refseq} ${mrnaSeq} ${mrnaLen} ${refLen} ${trimSeq} ${mrna2bit} ${mrnaAlignDir}/split.done trim.len 

${mrnaSeq}:
	@mkdir -p ${dir $@}
	rm -f ${mrnaSeq} 
	${gbBin}/gbGetSeqs -db=${GBID} -inclVersion -native -gbRoot=${gbRoot} genbank mrna stdout | tr acgt ACGT > ${mrnaSeq}.${tmpExt} 
	mv -f ${mrnaSeq}.${tmpExt} $@
${refseq}:
	@mkdir -p ${dir $@}
	rm -f ${refseq}
	${gbBin}/gbGetSeqs -db=${GBID} -inclVersion -native -gbRoot=${gbRoot} refSeq mrna stdout | tr acgt ACGT > ${refseq}.${tmpExt}
	mv -f ${refseq}.${tmpExt} $@
${mrna2bit}: ${mrnaSeq} ${refseq}
	cat ${mrnaSeq} ${refseq} | faToTwoBit stdin ${mrna2bit}
${trimSeq} : ${mrnaSeq} ${refseq} 
	@mkdir -p ${tmpDir}
	@mkdir -p ${mrnaAlignDir}
	cat ${mrnaSeq} ${refseq} | faTrimPolyA stdin stdout | gzip -c > ${trimSeq}.${tmpExt}
	mv -f ${trimSeq}.${tmpExt} ${trimSeq}
${mrnaAlignDir}/split.done: ${trimSeq} 
	faSplit about ${trimSeq} ${splitSize} ${mrnaAlignDir}/mrna
	ls ${mrnaAlignDir}/*.fa > ${mrnaAlignDir}/split.done
${mrnaLen} ${refLen} trim.len: ${trimSeq} ${mrnaSeq}
	faSize ${trimSeq} -detailed > trim.len
	faSize ${mrnaSeq} -detailed > mrna.len
	cp trim.len ${mrnaAlignDir}
	cp mrna.len ${mrnaAlignDir}

