#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser

class Chain:
    def __init__(self, line):
        items = line.strip().split()
        self.chr = items[7]
        self.chrsize = int(items[8])
        self.strand = items[9]
        self.start = int(items[10])
        self.end = int(items[11])
        self.id = items[12]
        self.starts = [0]
        self.sizes = []

    def convertToPosStrand(self):
        if self.strand == "-":
            self.strand = "+"
            start = self.start
            self.start = self.chrsize - self.end
            self.end = self.chrsize - start
            #for i in range(len(self.starts)):
            #    self.starts[i] = self.chrsize - (self.starts[i] + self.sizes[i])

    def getsize(self):
        return sum(self.sizes)

    def getReverse(self):
        if self.strand == '+':
            self.strand = '-'
        else:
            self.strand = '+'

        size = self.end - self.start

        temp = self.end
        self.end = self.chrsize - self.start
        self.start = self.chrsize - temp
        for i in range( len(self.starts) ):
            self.starts[i] = size - (self.starts[i] + self.sizes[i])
        starts = []
        sizes = []
        for i in range( len(self.starts) -1, -1, -1):
            starts.append( self.starts[i] )
            sizes.append( self.sizes[i] )
        self.starts = starts
        self.sizes = sizes

    def getStr(self):
        starts = ",".join([str(s) for s in self.starts])
        sizes = ",".join([str(s) for s in self.sizes])
        return "%s %d %s %d %d %s %s %s" %(self.chr, self.chrsize, self.strand, self.start, self.end, self.id, starts, sizes)


class Bed:
    def __init__(self, line):
        items = line.strip().split()
        self.chr = items[0]
        self.start = int(items[1])
        self.end = int(items[2])
        self.name = items[3]
        self.strand = items[5]
        self.blockCount = int(items[9])
        self.blockSizes = items[10].strip(',').split(',')
        self.blockStarts = items[11].strip(',').split(',')
        for i in range(self.blockCount):
            self.blockSizes[i] = int(self.blockSizes[i])
            self.blockStarts[i] = int(self.blockStarts[i])

    def getStr(self):
        blockSizes = ",".join( [ str(s) for s in self.blockSizes] )
        blockStarts = ",".join( [ str(s) for s in self.blockStarts] )
        return "%s %d %d %s %s %d %s %s" %(self.chr, self.start, self.end, self.name, self.strand, self.blockCount, blockSizes, blockStarts) 

def readBeds(file):
    f = open(file, "r")
    beds = []
    for line in f.readlines():
        bed = Bed(line)
        beds.append(bed)

    return beds

def checkChain(chain, beds):
    if chain.strand == '-':
        chain.getReverse()

    check = False
    for bed in beds:
        if chain.chr == bed.chr and chain.start == bed.start:
            check = True
            if chain.end != bed.end:
                check = False
                sys.stdout.write("Different ends, chainEnd = %d != bedEnd = %d\n" %(chain.end, bed.end))
                continue
            if len(chain.sizes) != bed.blockCount:
                check = False
                sys.stdout.write("Different blockCount, chainBlockCount = %d, bedBlockcount = %d\n" %(len(chain.sizes), bed.blockCount))
                continue
            for i in range(bed.blockCount):
                if chain.sizes[i] != bed.blockSizes[i]:
                    check = False
                    sys.stdout.write("Different blockSize, block %d, chainSize: %d, bedSize: %d\n" %(i, chain.sizes[i], bed.blockSizes[i]))
                    continue
                if chain.starts[i] != bed.blockStarts[i]:
                    check = False
                    sys.stdout.write("Different blockStart, block %d, chainStart: %d, bedStart: %d\n" %(i, chain.starts[i], bed.blockStarts[i]))
                    sys.stdout.write("Chain: %s\n" %chain.getStr())
                    sys.stdout.write("Bed: %s\n" %bed.getStr())
                    continue
            return check
    sys.stdout.write("Could not find any bed with start %d\n" %(chain.start))
    return check

def readChains(file, beds):
    f = open(file, "r")
    totalSize = 0
    chain = Chain(f.readline())
    for line in f.readlines():
        if re.search("chain", line):
            #check prev chain:
            check = checkChain(chain, beds)
            totalSize += chain.getsize()
            #if check:
            #    sys.stdout.write("%s\n" %(chain.id))
            chain = Chain(line) 
        elif re.search("\S", line):
            items = line.strip().split()
            blocksize = int(items[0])
            chain.sizes.append(blocksize)
            if len(items) == 3:
                prevstart = chain.starts[len(chain.starts) -1]
                start = prevstart + blocksize + int(items[2])
                chain.starts.append(start)
            
    checkChain(chain, beds)
    totalSize += chain.getsize()
    f.close()
    
    sys.stdout.write("Total chain size: %d\n" %totalSize)

    return


def main():
    beds = readBeds(sys.argv[1])
    readChains(sys.argv[2], beds)
    

if __name__ == "__main__":
    main()


