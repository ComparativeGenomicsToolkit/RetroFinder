
DEF=${mrnaAlignDir}/DEF
SED=${mrnaAlignDir}/sed.sh

mrnaAlign: ${mrnaAlignDir}/ma.jobs.done ${mrnaAlignDir}/out2lav.jobs.done ${mrnaAlignDir}/sed.done

# run parasol jobs (don't check clusterDir dependences if done file exists)
alignParaDir = para.ma.${tmpExt}
out2lavParaDir = run.1.${tmpExt}
#runMrnaAlignJobs: ${mrnaAlignDir}/ma.jobs.done

${DEF}:
	@mkdir -p ${dir $@}
	echo "export PATH=/parasol/bin:/usr/bin:/bin:/usr/local/bin:${retroBin}:/cluster/bin/penn/x86_64:/cluster/bin/scripts:/cluster/bin/x86_64" > ${DEF}
	echo "ALIGN=blastz-run">>${DEF}
	echo "BLASTZ=blastz">>${DEF}
	echo "BLASTZ_H=2000">>${DEF}
	echo "BLASTZ_ABRIDGE_REPEATS=0">>${DEF}
	echo "SEQ1_DIR=/scratch/data/hg18/bothMaskedNibs/">>${DEF}
	echo "SEQ1_RMSK=">>${DEF}
	echo "SEQ1_SMSK=">>${DEF}
	echo "SEQ1_FLAG=-primate">>${DEF}
	echo "SEQ1_IN_CONTIGS=0">>${DEF}
	echo "SEQ1_CHUNK=10000000">>${DEF}
	echo "SEQ1_LAP=10000">>${DEF}
	echo "SEQ2_DIR=${mrnaAlignDir}">>${DEF}
	echo "SEQ2_RMSK=">>${DEF}
	echo "SEQ2_SMSK=">>${DEF}
	echo "SEQ2_FLAG=-primate">>${DEF}
	echo "SEQ2_IN_CONTIGS=1">>${DEF}
	echo "SEQ2_CHUNK=30000000">>${DEF}
	echo "SEQ2_LAP=0">>${DEF}
	echo "BASE=${mrnaAlignDir}">>${DEF}
	echo "DEF=${mrnaAlignDir}/DEF">>${DEF}
	echo "RAW=${mrnaAlignDir}/raw">>${DEF}
	echo "CDBDIR=${mrnaAlignDir}">>${DEF}
	echo "SEQ1_LEN=${mrnaAlignDir}/S1.len">>${DEF}
	echo "SEQ2_LEN=${tmpDir}/trim.len">>${DEF}
	echo "DEBUG=1">>${DEF}

${alignParaDir}/jobs.para: ${DEF}
	@mkdir -p ${dir $@}
	@rm -f $@ $@.${tmpExt}
	${retroBin}/make-joblist ${DEF} >$@.${tmpExt}
	rm -rf ${mrnaAlignDir}/raw
	bash ${mrnaAlignDir}/xdir.sh
	mv -f $@.${tmpExt} $@


${alignParaDir}/batch: ${alignParaDir}/jobs.para ${DEF}
	ssh ${paraHost} cd `pwd`/${alignParaDir} \; para create jobs.para < /dev/null

${mrnaAlignDir}/ma.jobs.done: ${alignParaDir}/batch
	@mkdir -p ${dir $@}
	${MAKE} makeMaJobs
	touch $@

makeMaJobs: ${alignParaDir}/jobs.para
	ssh ${paraHost} cd `pwd`/${alignParaDir} \; para push < /dev/null


${out2lavParaDir}/jobs.para: ${DEF}
	@mkdir -p ${dir $@}
	@rm -f $@ $@.${tmpExt}
	rm -rf ${mrnaAlignDir}/lav
	${scripts}/do.out2lav ${DEF} >$@.${tmpExt}
	mv -f $@.${tmpExt} $@


${out2lavParaDir}/batch: ${out2lavParaDir}/jobs.para ${DEF}
	ssh ${paraHost} cd `pwd`/${out2lavParaDir} \; para create jobs.para < /dev/null

${mrnaAlignDir}/out2lav.jobs.done: ${out2lavParaDir}/batch
	@mkdir -p ${dir $@}
	${MAKE} makeOut2LavJobs
	touch $@

makeOut2LavJobs: ${out2lavParaDir}/jobs.para
	ssh ${paraHost} cd `pwd`/${out2lavParaDir} \; para push < /dev/null

${mrnaAlignDir}/sed.done: 
	${scripts}/sed.sh ${mrnaAlignDir} ${db}
	touch $@

${mrnaAlignDir}/chainFilter/%.chain: ${mrnaAlignDir}/pslRaw/%.psl
	doSortFilterChain.sh $c ${mrnaAlignDir} ${DB} ${trimSeq}

# combine all chroms:
cgCombined: ${cleanGenesAnalDir}/summary.tsv ${cleanGenesAnalDir}/details.tsv

${cleanGenesAnalDir}/%.tsv: $(foreach c,${chroms},${cleanGenesAnalChromDir}/$c/all/%.tsv)
	@mkdir -p ${dir $@}
	(head -1 $<; tail -q -n +2 $^) >$@.${tmpExt}
	mv -f $@.${tmpExt} $@
${cleanGenesAnalChromDir}/%.tsv: ${cleanGenesDir}/anal.jobs.done
