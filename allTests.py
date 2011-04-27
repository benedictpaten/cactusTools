#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from cactusTools.graphVizPlots.cactus_graphVizTest import TestCase as graphVizTest
from cactusTools.mafs.cactus_MAFGeneratorTest import TestCase as mAFTest
from cactusTools.stats.cactus_treeStatsTest import TestCase as statsTest
from cactusTools.referenceUtils.cactus_addReferenceSequenceTest import TestCase as referenceSequenceTest
from cactusTools.referenceViewer.cactus_referenceViewerTest import TestCase as referenceViewerTest

from sonLib.bioio import parseSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(graphVizTest, 'test'),
                                   unittest.makeSuite(mAFTest, 'test'),
                                   unittest.makeSuite(statsTest, 'test'),
                                   unittest.makeSuite(referenceViewerTest, 'test'),
                                   unittest.makeSuite(referenceSequenceTest, 'test')))
    return allTests
        
def main():
    parseSuiteTestOptions()
     
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)

if __name__ == '__main__':
    import sys
    sys.exit(main())
 
