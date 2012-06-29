
# define variable progs with base programs names before including

root = ../..
copts = -g
MACHTYPE = x86_64
mysqllibs=$(shell mysql_config --libs)

bindir = ${root}/bin

kentsrc = ${HOME}/kent/src
kentlibdir = ${kentsrc}/lib/${MACHTYPE}

SAMDIR = /hive/data/outside/samtools/samtools-0.1.18/${MACHTYPE}
SAMLIB = ${SAMDIR}/libbam.a

TABIXDIR = /hive/data/outside/tabix/tabix-0.2.5/${MACHTYPE}
TABIXLIB = ${TABIXDIR}/libtabix.a

INCL= -I${kentsrc}/inc -I${kentsrc}/hg/inc
LIBS = ${kentlibdir}/jkhgap.a ${kentlibdir}/jkweb.a -lm ${mysqllibs} ${SAMLIB} ${TABIXLIB} -lpthread -lssl

CFLAGS = ${copts} ${INCL}

progpaths = ${progs:%=${bindir}/%}

# rules
all: ${progpaths}

clean::
	rm -f ${progpaths}

${bindir}/%: %.c
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@ $< ${LIBS}
