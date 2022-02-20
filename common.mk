
# define variable progs with base programs names before including

root = ../..
copts = -g
MACHTYPE = x86_64
mysqllibs=$(shell mysql_config --libs)

bindir = ${root}/bin

kentsrc = ${HOME}/kent/src
kentlibdir = ${kentsrc}/lib/${MACHTYPE}

# autodetect where libm is installed
ifeq (${MLIB},)
  ifneq ($(wildcard /usr/lib64/libm.a),)
      MLIB=-lm
  endif
endif
ifeq (${MLIB},)
  MLIB=-lm
endif

# autodetect UCSC installation of hal:
ifeq (${HALDIR},)
    HALDIR = /hive/groups/browser/hal/build/hal.2020-12-18
    ifneq ($(wildcard ${HALDIR}),)
        ifeq (${USE_HAL},)
          USE_HAL=1
        endif
    endif
endif
ifeq (${USE_HAL},1)
    # force static libraries to keep programs portable
    HDF5DIR=/hive/groups/browser/hal/build/hdf5-1.12.0
    HDF5LIBDIR=${HDF5DIR}/local/lib
    HDF5LIBS=${HDF5LIBDIR}/libhdf5_cpp.a ${HDF5LIBDIR}/libhdf5.a ${HDF5LIBDIR}/libhdf5_hl.a
    HALLIBS=${HALDIR}/hal/lib/libHalBlockViz.a ${HALDIR}/hal/lib/libHalMaf.a ${HALDIR}/hal/lib/libHalLiftover.a ${HALDIR}/hal/lib/libHalLod.a ${HALDIR}/hal/lib/libHal.a ${HALDIR}/sonLib/lib/sonLib.a ${HDF5LIBS} -lcurl -lstdc++
    HG_DEFS+=-DUSE_HAL
    HG_INC+=-I${HALDIR}/inc -I${HALDIR}/hal/blockViz/inc
endif

HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}
HG_INC+=-I../inc -I../../inc -I../../../inc -I../../../../inc -I../../../../../inc

MYSQLLIBS := $(shell mysql_config --libs || true)

# OK to add -lstdc++ to all MYSQLLIBS just in case it is
#    MySQL version 5.6 libraries, but no 'librt' on Mac OSX
ifeq (${IS_HGWDEV},yes)
  MYSQLLIBS += /usr/lib/gcc/x86_64-redhat-linux/4.8.5/libstdc++.a /usr/lib64/librt.a
else
  ifeq ($(UNAME_S),Darwin)
    MYSQLLIBS += -lstdc++
  else
    MYSQLLIBS += -lstdc++ -lrt
  endif
endif
FREETYPELIBS =  $(shell freetype-config --libs 2> /dev/null )


INCL= -I${kentsrc}/inc -I${kentsrc}/hg/inc
LIBS = ${kentlibdir}/jkhgapcgi.a ${kentlibdir}/jkhgap.a ${kentlibdir}/jkweb.a -lm ${MYSQLLIBS} $(kentsrc)/htslib/libhts.a ${FREETYPELIBS} ${HALLIBS} -ldl -lpthread -lssl -lcrypto -lz

CFLAGS = ${copts} ${INCL}

progpaths = ${progs:%=${bindir}/%}

# rules
all: ${progpaths}

clean::
	rm -f ${progpaths}

${bindir}/%: %.c
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@ $< ${LIBS}
