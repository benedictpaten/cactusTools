cactusRootPath=${rootPath}../cactus
include  ${cactusRootPath}/include.mk
#Location of bin and lib dirs
binPath=${rootPath}bin
libPath=${rootPath}lib
cactusLibPath=${cactusRootPath}/lib

cflags += -I ${cactusRootPath}/lib
#cflags += -I ${cactusRootPath}/lib ${cflags_prof}

# optional kent library stuff
KENTDIR=/hive/groups/recon/local/kent/src
ifneq ($(wildcard ${KENTDIR}),)
    kentInc = ${KENTDIR}/inc
    kentLib = ${KENTDIR}/lib
    kentLibWeb = ${kentLib}/${MACH}/jkweb.a
endif
