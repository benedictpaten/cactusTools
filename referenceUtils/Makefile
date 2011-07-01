rootPath = ../
include ../include.mk

progs = $(notdir $(wildcard cactus_addReferenceSequenceTest.py))

targets = ${binPath}/cactus_addReferenceSeq ${binPath}/cactus_getReferenceSeq

all: ${targets} ${progs:%=${binPath}/%}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_addReferenceSeq: cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}
	#${cxx} ${cflags} -g -o1 -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}
	#${cxx} ${cflags} -pg -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_getReferenceSeq: cactus_getReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${targets}
