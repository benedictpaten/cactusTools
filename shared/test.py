#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
"""Common test functions used for generating inputs to run cactus workflow and running
cactus workflow and the various utilities.
"""

import random
import os
import xml.etree.ElementTree as ET

from sonLib.bioio import logger
from sonLib.bioio import getTempFile
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getRandomSequence
from sonLib.bioio import mutateSequence
from sonLib.bioio import reverseComplement
from sonLib.bioio import fastaWrite
from sonLib.bioio import printBinaryTree
from sonLib.bioio import system
from sonLib.bioio import getRandomAlphaNumericString
from sonLib.bioio import runGraphViz

from cactusTools.shared.common import runCactusTreeViewer
from cactusTools.shared.common import runCactusAdjacencyGraphViewer
from cactusTools.shared.common import runCactusTreeStats
from cactusTools.shared.common import runCactusMAFGenerator
from cactusTools.shared.common import runCactusTreeStatsToLatexTables

from sonLib.bioio import TestStatus

from sonLib.tree import makeRandomBinaryTree

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete
from jobTree.src.common import runJobTreeStats

from cactus.shared.config import checkDatabaseConf
from cactus.shared.config import CactusWorkflowExperiment
import cactus.shared.test

def runWorkflow_TestScript(sequences, newickTreeString, 
                           batchSystem="single_machine",
                           buildTrees=True, buildFaces=True, buildReference=True,
                           buildCactusPDF=False,
                           buildAdjacencyPDF=False,
                           makeCactusTreeStats=False, 
                           makeMAFs=False, 
                           configFile=None,
                           buildJobTreeStats=False):
    """Runs the workflow and various downstream utilities.
    """
    
    #Setup the temp dir
    tempDir = getTempDirectory(".")
    logger.info("Using the temp dir: %s" % tempDir)
        
    #Setup the output dir
    outputDir = getTempDirectory(tempDir)
    logger.info("Using the output dir: %s" % outputDir)
    
    experiment = cactus.shared.test.runWorkflow_TestScript(sequences, newickTreeString, 
                           outputDir=outputDir, 
                           batchSystem=batchSystem,
                           buildTrees=buildTrees, buildFaces=buildFaces, buildReference=buildReference,
                           configFile=configFile,
                           buildJobTreeStats=buildJobTreeStats)
    cactusDiskDatabaseString = experiment.getDatabaseString()
    
    #Run the cactus tree graph-viz plot
    if buildCactusPDF:
        cactusTreeDotFile = os.path.join(outputDir, "cactusTree.dot")
        cactusTreePDFFile = os.path.join(outputDir, "cactusTree.pdf")
        runCactusTreeViewer(cactusTreeDotFile, cactusDiskDatabaseString)
        runGraphViz(cactusTreeDotFile, cactusTreePDFFile)
        logger.info("Ran the cactus tree plot script")
    else:
        logger.info("Not building a cactus tree plot")
    
    #Run the cactus tree graph-viz plot
    if buildAdjacencyPDF:
        adjacencyGraphDotFile = os.path.join(outputDir, "adjacencyGraph.dot")
        adjacencyGraphPDFFile = os.path.join(outputDir, "adjacencyGraph.pdf")
        runCactusAdjacencyGraphViewer(adjacencyGraphDotFile, cactusDiskDatabaseString)
        runGraphViz(adjacencyGraphDotFile, adjacencyGraphPDFFile)
        logger.info("Ran the adjacency graph plot script")
    else:
        logger.info("Not building a adjacency graph plot")
    
    if makeCactusTreeStats:
        cactusTreeFile = os.path.join(outputDir, "cactusStats.xml")
        runCactusTreeStats(cactusTreeFile, cactusDiskDatabaseString)
        #Now run the latex script
        statsFileTEX = os.path.join(outputDir, "cactusStats.tex")
        runCactusTreeStatsToLatexTables([ cactusTreeFile ], [ "region0" ], statsFileTEX)
        logger.info("Ran the tree stats script")
    else:
        logger.info("Not running cactus tree stats")
    
    if makeMAFs:
        mAFFile = os.path.join(outputDir, "cactus.maf")
        runCactusMAFGenerator(mAFFile, cactusDiskDatabaseString)
        logger.info("Ran the MAF building script")
    else:
        logger.info("Not building the MAFs")
        
    #Now remove everything we generate
    experiment.cleanupDatabase()
    system("rm -rf %s" % tempDir)    
        
def runWorkflow_multipleExamples(inputGenFunction,
                                 testNumber=1, 
                                 testRestrictions=(TestStatus.TEST_SHORT, TestStatus.TEST_MEDIUM, \
                                                   TestStatus.TEST_LONG, TestStatus.TEST_VERY_LONG,),
                               inverseTestRestrictions=False,
                               batchSystem="single_machine",
                               buildTrees=True, buildFaces=True, buildReference=True,
                               buildReferenceSequence=False,
                               buildCactusPDF=False, buildAdjacencyPDF=False,
                               buildReferencePDF=False,
                               makeCactusTreeStats=False, makeMAFs=False,
                               configFile=None, buildJobTreeStats=False):
    """A wrapper to run a number of examples.
    """
    if (inverseTestRestrictions and TestStatus.getTestStatus() not in testRestrictions) or \
        (not inverseTestRestrictions and TestStatus.getTestStatus() in testRestrictions):
        for test in xrange(testNumber): 
            tempDir = getTempDirectory(os.getcwd())
            sequences, newickTreeString = inputGenFunction(regionNumber=test, tempDir=tempDir)
            runWorkflow_TestScript(sequences, newickTreeString,
                                   batchSystem=batchSystem,
                                   buildTrees=buildTrees, buildFaces=buildFaces, buildReference=buildReference, 
                                   buildCactusPDF=buildCactusPDF, buildAdjacencyPDF=buildAdjacencyPDF,
                                   makeCactusTreeStats=makeCactusTreeStats, makeMAFs=makeMAFs, configFile=configFile,
                                   buildJobTreeStats=buildJobTreeStats)
            system("rm -rf %s" % tempDir)
            logger.info("Finished random test %i" % test)
    
