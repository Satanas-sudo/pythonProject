import unittest
from pathlib import Path
import interactomeOOP
import relatedComponents

file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")
minifile = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/toy_example.txt")


class TestInteractome(unittest.TestCase):

    #def tests_Interactome(self):
        #TO DO

    def tests_relatedComponents(self):

        relatedComponents.extractCC('A')
        relatedComponents.countCC()
        relatedComponents.writeCC()
        relatedComponents.computeCC()
        print(relatedComponents.clustering('A'))
        print(relatedComponents.density())

    def test_extractCC(self):
        t = relatedComponents.extractCC('A')
        self.assertEqual(t, ['D', 'B', 'E', 'F', 'A', 'C'])

    def test_clustering(self):
        t3 = relatedComponents.clustering('A')
        self.assertEqual(type(t3), float)


if __name__ == '__main__':
    unittest.main()
