#!/usr/bin/env python

#Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest
import sys

from cactus.shared.test import parseCactusSuiteTestOptions
from sonLib.bioio import TestStatus

from cactus.shared.test import getCactusInputs_random
from cactus.shared.test import getCactusInputs_blanchette

from cactusTools.shared.test import runWorkflow_multipleExamples

class TestCase(unittest.TestCase):
    def testCactusReferenceSequenceAdd_Random(self):
        runWorkflow_multipleExamples(getCactusInputs_random,
                                     testNumber=TestStatus.getTestSetup(),
                                     buildTrees=False, buildFaces=False, buildReference=True,
                                     buildReferenceSequence=True)
        
    def testCactusReferenceSequenceAdd_Blanchette(self):
        runWorkflow_multipleExamples(getCactusInputs_blanchette,
                                     testRestrictions=(TestStatus.TEST_SHORT,), inverseTestRestrictions=True,
                                     buildTrees=False, buildFaces=False, buildReference=True,
                                     buildReferenceSequence=True)
    
def main():
    parseCactusSuiteTestOptions() 
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
