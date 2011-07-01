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
from sonLib.bioio import getLogLevelString

def getLogLevelString2(logLevelString):
    """Gets the log level string for the binary
    """
    if logLevelString == None:
        return getLogLevelString()
    return logLevelString
    
#############################################
#############################################
#Runs cactus utilities.
#############################################
#############################################    
    
def runCactusTreeStats(outputFile, cactusDiskDatabaseString, flowerName='0', logLevel=None, referenceEventString=None):
    logLevel = getLogLevelString2(logLevel)
    referenceEventString = nameValue("referenceEventString", referenceEventString, str)
    command = "cactus_treeStats --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s %s" % (cactusDiskDatabaseString, flowerName, outputFile, logLevel, referenceEventString)
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
                        logLevel=None):
    logLevel = getLogLevelString2(logLevel)
    system("cactus_treeViewer --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDiskDatabaseString, flowerName, graphFile, logLevel))
    logger.info("Created a cactus tree graph")
    
def runCactusAdjacencyGraphViewer(graphFile,
                             cactusDiskDatabaseString, flowerName="0",
                             logLevel=None, includeInternalAdjacencies=False):
    logLevel = getLogLevelString2(logLevel)
    includeInternalAdjacencies = nameValue("includeInternalAdjacencies", includeInternalAdjacencies, bool)
    system("cactus_adjacencyGraphViewer --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s" \
                    % (cactusDiskDatabaseString, flowerName, graphFile, logLevel))
    logger.info("Created a break point graph of the problem")
    
def runCactusMAFGenerator(mAFFile, cactusDiskDatabaseString, flowerName="0",
                          logLevel=None, referenceEventString=None):
    logLevel = getLogLevelString2(logLevel)
    referenceEventString = nameValue("referenceEventString", referenceEventString, str)
    system("cactus_MAFGenerator --cactusDisk '%s' --flowerName %s --outputFile %s --logLevel %s %s" \
            % (cactusDiskDatabaseString, flowerName, mAFFile, logLevel, referenceEventString))
    logger.info("Created a MAF for the given cactusDisk")
