#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys
import os

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory

from cactus.shared.test import getCactusInputs_random
from cactusTools.shared.test import runWorkflow_multipleExamples
from cactusTools.shared.common import runCactusTreeStatsToLatexTables

class TestCase(unittest.TestCase):
    def testCactus_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random, 
                                     testNumber=TestStatus.getTestSetup(), 
                                     makeCactusTreeStats=True)
        
def main():
    parseCactusSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
