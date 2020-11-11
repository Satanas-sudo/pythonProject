import os
import glob


dataDir = "data"
solutionDir = "solution"

generatedCode = """#! /usr/bin/env python3

import unittest
import usefulLinesCounter

class TestUsefulLineCounter(unittest.TestCase):

"""


usecaseFilePaths = glob.glob(dataDir + os.path.sep + "*.txt")

for currentUsecaseFilePath in usecaseFilePaths:
	solutionFilePath = currentUsecaseFilePath.replace(dataDir, solutionDir)
	if os.path.exists(solutionFilePath):
		with open(solutionFilePath) as solutionFile:
			currentSolution = solutionFile.read().strip()
			solutionFileName = currentUsecaseFilePath.replace(dataDir + os.path.sep, "").replace(".txt", "")
			generatedCode += """
	def test_""" + solutionFileName + """(self):
		self.assertEqual(usefulLinesCounter.countUsefulLines('""" + currentUsecaseFilePath + """'), """ + currentSolution + """, "Should be """ + currentSolution + """")
"""


generatedCode += """


if __name__ == '__main__':
    unittest.main()
"""

with open("testUsefulLineCounterMultipleTests-generated.py", "w") as generatedFile:
	generatedFile.write(generatedCode)
