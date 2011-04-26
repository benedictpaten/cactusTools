rootPath = ./
include ./include.mk

# order is important, libraries first
modules = genemap graphVizPlots mafs stats referenceViewer psls tuning beds referenceUtils 
.PHONY: all %.all clean %.clean

all : ${libPath}/cactusUtils.a ${modules:%=all.%}

${libPath}/cactusUtils.a : cactusUtils.h cactusUtils.c ${basicLibsDependencies}
	${cxx} ${cflags} -c cactusUtils.c  -I ${libPath}
	ar rc cactusUtils.a *.o
	ranlib cactusUtils.a
	rm *.o
	mv cactusUtils.a ${libPath}/
	cp cactusUtils.h ${libPath}/

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%} clean.cactusUtils

clean.%:
	cd $* && make clean
	rm -rf ${binPath}/*.dSYM

clean.cactusUtils:
	rm -f ${libPath}/cactusUtils.a ${libPath}/cactusUtils.h
	
test :
	python allTests.py
 