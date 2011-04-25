rootPath = ../
include ../include.mk

all :  ${binPath}/cactus_referenceViewer

${binPath}/cactus_referenceViewer : *.c *.h ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_referenceViewer main.c referenceViewer.c graphViz.c ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${binPath}/cactus_referenceViewer
