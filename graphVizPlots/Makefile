rootPath = ../
include ../include.mk

all :  ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 

${binPath}/cactus_adjacencyGraphViewer : *.c ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_adjacencyGraphViewer cactus_adjacencyGraphViewer.c ${cactusLibPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_treeViewer : *.c  ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeViewer cactus_treeViewer.c ${cactusLibPath}/cactusLib.a ${basicLibs}

clean :
	rm -rf *.o
	rm -rf ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 
