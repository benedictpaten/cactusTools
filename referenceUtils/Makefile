rootPath = ../
include ../include.mk

progs = $(notdir $(wildcard))

targets = ${binPath}/cactus_addReferenceSeq ${binPath}/cactus_getReferenceSeq

all: ${targets}

${binPath}/cactus_addReferenceSeq: cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_getReferenceSeq: cactus_getReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${targets}
