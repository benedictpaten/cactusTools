cactusLibRootPath=${rootPath}../cactus
include  ${cactusLibRootPath}/include.mk
#Location of bin and lib dirs
binPath=${rootPath}/bin
libPath=${cactusLibRootPath}/lib

# optional kent library stuff
KENTDIR=/hive/groups/recon/local/kent/src
ifneq ($(wildcard ${KENTDIR}),)
    kentInc = ${KENTDIR}/inc
    kentLib = ${KENTDIR}/lib
    kentLibWeb = ${kentLib}/${MACH}/jkweb.a
endif
