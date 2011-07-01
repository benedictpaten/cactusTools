#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Wrapper functions for assisting in running the various programs of the cactus tools package.
"""

import os
import random

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import system
from sonLib.bioio import nameValue
    
#############################################
#############################################
#Runs cactus utilities.
#############################################
#############################################    
    
def runCactusTreeStats(outputFile, cactusDiskDatabaseString, flowerName='0'):
    command = "cactus_treeStats --cactusDisk '%s' --flowerName %s --outputFile %s" % (cactusDiskDatabaseString, flowerName, outputFile)
    system(command)
    logger.info("Ran the cactus tree stats command apprently okay")

def runCactusTreeStatsToLatexTables(inputFiles, regionNames, outputFile):
    assert len(regionNames) == len(inputFiles)
    k = " ".join([ "%s %s" % (i, j) for i, j in zip(inputFiles, regionNames) ])
    command = "cactus_treeStatsToLatexTables.py --outputFile %s %s" % (outputFile, k)
    system(command)
    logger.info("Ran cactus_treeStatsToLatexTables okay")
    
def runCactusTreeViewer(graphFile,
                        cactusDiskDatabaseString, 
                        flowerName="0", 
                        logLevel="DEBUG"):
    system("cactus_treeViewer --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDiskDatabaseString, flowerName, graphFile, logLevel))
    logger.info("Created a cactus tree graph")
    
def runCactusAdjacencyGraphViewer(graphFile,
                             cactusDiskDatabaseString, flowerName="0",
                             logLevel="DEBUG", includeInternalAdjacencies=False):
    includeInternalAdjacencies = nameValue("includeInternalAdjacencies", includeInternalAdjacencies, bool)
    system("cactus_adjacencyGraphViewer --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDiskDatabaseString, flowerName, graphFile, logLevel))
    logger.info("Created a break point graph of the problem")
    
def runCactusReferenceGraphViewer(graphFile,
                                  cactusDiskDatabaseString, flowerName="0",
                                  logLevel="DEBUG"):
    system("cactus_referenceViewer --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDiskDatabaseString, flowerName, graphFile, logLevel))
    logger.info("Created a cactus reference graph")
    
def runCactusMAFGenerator(mAFFile, cactusDiskDatabaseString, flowerName="0",
                          logLevel="DEBUG", orderByReference=False, referenceSequenceName=None):
    orderByReference = nameValue("orderByReference", orderByReference, bool)
    referenceSequence = nameValue("referenceSequence", referenceSequenceName, str)
    system("cactus_MAFGenerator --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s %s %s" \
            % (cactusDiskDatabaseString, flowerName, mAFFile, logLevel, orderByReference, referenceSequence))
    logger.info("Created a MAF for the given cactusDisk")
    
def runCactusAddReferenceSequence(cactusDiskDatabaseString, flowerName="0",
                          logLevel="DEBUG", referenceEventString="reference"):
    system("valgrind --leak-check=yes cactus_addReferenceSeq --cactusDisk '%s' --logLevel %s --referenceEventString %s" % (cactusDiskDatabaseString, logLevel, referenceEventString))
    #system("cactus_addReferenceSeq --cactusDisk '%s' --logLevel %s --referenceEventString %s" % (cactusDiskDatabaseString, logLevel, referenceEventString))
