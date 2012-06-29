subdirs = bedOverlap genePredFilter hgCountTissueMicroarray \
	  orfBatch pslFilterDups pslNetOverlap pslPseudo pslSplit

# not used:  countFrameShift

all: ${subdirs:%=%.all}

%.all:
	(cd $* && ${MAKE} all)

clean: ${subdirs:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)
