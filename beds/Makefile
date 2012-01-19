rootPath = ../
include ../include.mk
refDir = ../referenceUtils

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_bedGenerator ${binPath}/cactus_chain

all: ${targets}

${binPath}/cactus_bedGenerator : cactus_bedGenerator.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_bedGenerator cactus_bedGenerator.c ${libPath}/cactusUtils.a ${cactusLibPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_chain: cactus_chain.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_chain cactus_chain.c ${libPath}/cactusUtils.a ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${binPath}/cactus_bedGenerator ${binPath}/cactus_chain
