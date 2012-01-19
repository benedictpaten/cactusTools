rootPath = ../
include ../include.mk
refDir = ../referenceUtils

progs = $(notdir $(wildcard cactus_mafToReferenceSeq.py))
targets = ${progs:%=${binPath}/%}

all :  ${targets} ${binPath}/cactus_MAFGenerator ${libPath}/cactusMafs.a ${binPath}/cactus_augmentedMaf

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_MAFGenerator : *.c *.h ${libPath}/cactusTraversal.a ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_MAFGenerator cactus_MAFGenerator.c mafs.c ${libPath}/cactusTraversal.a ${libPath}/cactusUtils.a ${cactusLibPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_augmentedMaf :  *.c *.h ${libPath}/cactusUtils.h ${libPath}/cactusUtils.a cactus_augmentedMaf.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_augmentedMaf cactus_augmentedMaf.c ${libPath}/cactusUtils.a ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${binPath}/cactus_augmentedMaf ${binPath}/cactus_MAFGenerator ${libPath}/cactusMafs.a ${libPath}/cactusMafs.h  
	rm -rf ${progs:%=${binPath}/%}
	rm -rf ${libPath}/cactusMafs.h 
	
${libPath}/cactusMafs.a :  *.c *.h  ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c mafs.c
	ar rc cactusMafs.a *.o 
	ranlib cactusMafs.a 
	rm *.o
	mv cactusMafs.a ${libPath}/
	cp cactusMafs.h ${libPath}/ 
