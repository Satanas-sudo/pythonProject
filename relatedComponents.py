import pathlib
from pathlib import Path
import interactomeOOP as main

file = pathlib.Path('C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python', 'Human_HighQuality.txt')
#file = Path('C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt')
#file = Path('C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt')
minifile = Path('C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/toy_example.txt')


class Connexe(main.file):

    # Attributes
    list = []
    dict = {}
    int = []
    clean = []

    # Constructor
    def __init__(self, filepath):
        self.list = self.read_interaction_file_list(file)
        self.dict = self.read_interaction_file_dict(file)
        self.int = self.read_interaction_file(file)
        self.clean = self.clean_interactome(file)

    # Methods

    def extractCC(self, prot):
        """ This function returns every peaks of the proteins related component."""
        CC_list = [prot]
        for key in CC_list:
            for prot_cc in self.dict[key]:
                if prot_cc not in CC_list:
                        CC_list.append(prot_cc)  # Neighboor is added to the list
        return CC_list

    print(minifile.extractCC('C'))


    def countCC(self):
        """This function calculates the related components number of a graph, and gives for each of them their size (its
        number of proteins)."""

        prot_list = list(self.dict())
        CC_list = []
        CC_int = 0  # Counter initialization

        while len(prot_list) != 0:
            lcc = self.extractCC(prot_list[0])

            for CC in lcc:
                prot_list.remove(CC)  # Proteins already in a related component are removed

            CC_str = "The related component numbered" + str(CC_int)  # The components are called according to order of appearance
            CC_list.append((CC_str, len(lcc)))
            CC_int += 1
            print("The related component number", str(CC_int), "is composed of ", str(len(CC)), "proteins")

        return CC_list

    def computeCC(self):
        """This function returns a lcc list where each lcc[i] element corresponds to the related component number of the
         protein at a i location in the proteins list of the graph (which is a class attribute)."""
        lcc_list = []  # Initialization of the list containing the related component number for each protein

        for prot in self.proteins:
            lcc_list.append(self.extractCC(prot))
        return lcc_list


    def writeCC(self):
        """This function writes in a file the different related components of a graph. The format follows these rules :
        1. One line by related component.
        2. The first element of the line is the size of the related component.
        3. Then, the list of peaks that compose this related component will be added."""
        