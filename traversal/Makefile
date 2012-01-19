rootPath = ../
include ../include.mk

libSources = cactusTraversal.c 
libHeaders = cactusTraversal.h

all : ${libPath}/cactusTraversal.a 

${libPath}/cactusTraversal.a : ${libSources} ${libHeaders} ${cactusLibPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc cactusTraversal.a *.o
	ranlib cactusTraversal.a 
	rm *.o
	mv cactusTraversal.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f ${libPath}/cactusTraversal.*

