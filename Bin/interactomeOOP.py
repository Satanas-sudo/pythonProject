from pathlib import Path

file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")

class Interactome :

    # Attributes
    int_dict = {}
    int_list = []
    proteins = []

    # Constructor
    def __init__(self, file):
        self.init_list = self.read_interaction_file_list(file)
        self.init_dict = self.read_interaction_file_dict(file)
        self.proteins = self.proteins_list()

    # Methods

    def read_interaction_file_dict(self, file):
        """This function reads an interaction graph between proteins in a tabulated file and stores it in a dictionary.

        "Peaks" proteins are considered as keys and proteins that interacts with them are considered as the values
        associated to the keys.

        :param file: interaction file
        :return: prot_graph_dict: first element as key and their neighbor(s) in second element as the key's value(s)
        :rtype: dictionary
        """
        with open(file, "r") as interaction_file:
            lines_num_int = int(interaction_file.readline())  # File first line escapement

            data_str = interaction_file.readlines()  # Storing in memory the file's data
            prot_graph_dict = {}  # Dictionary declaration
            prot_name_str = ""  # Initialization of peak protein's name

            # Loop for each line of the document in memory (this loops needs that the file's data be classified)
            for line in data_str:

                # Checking for a key that already exists for the peak protein
                if prot_name_str == (line.split()[0]):
                    # Interacting protein addition in the list of associated values to the peak protein
                    prot_graph_dict[prot_name_str].append(
                        (line.split())[1])  # Cuts a string in a list where each word is
                    # an item of the list

                else:
                    prot_name_str = (line.split()[0])  # Storing in memory of the new peak protein (as first element)
                    prot_graph_dict[prot_name_str] = list()  # Creation of the corresponding key
                    prot_graph_dict[prot_name_str].append(
                        (line.split())[1])  # Interacting protein addition in the list of
                    # associated value (in second element)

        return prot_graph_dict

    print(read_interaction_file_dict("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                     "Human_HighQuality.txt"))

    def read_interaction_file_list(self, file):
        """This function reads an interaction graph between proteins in a file and stores it in a list of couples.

        :param file: interaction file
        :return: prot_graph_list: couples of interacting proteins
        :rtype: list
        """
        with open(file, "r") as interaction_file:
            lines_num_int = int(interaction_file.readline())

            data_str = interaction_file.readlines()
            prot_graph_list = []  # List creation

            for line in data_str:  # For each line of the interaction file (and so proteins couple), replace tabulations by
                # a backslash and adds them in the list
                prot_graph_list.append(line.replace("\t", "\ ", 1).replace("\n", "", 1))

        return prot_graph_list

    print(read_interaction_file_list("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                     "Human_HighQuality.txt"))


def read_interaction_file(self, file):
    """This function returns a couple where the first element is the dictionary representing the graph and the second
    element is the interaction list representing that same graph.

    This function calls the two prior functions and stores them in a tuple.

    :param file: interaction file
    :return: couple_d_l_tuple: first element is the dictionary and the second element is the list
    :rtype: tuple
    """
    with open(file, "r") as interaction_file:
        couple_d_l_tuple = [(self.read_interaction_file_dict(file)), (self.read_interaction_file_list(file))]  # Creation of the couple with dictionary as first element and list as second
        # element
    return couple_d_l_tuple

print(read_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                            "Human_HighQuality.txt"))


def count_vertices(self):
    """This function counts the number of peaks of a graph.

    :param file: interaction file
    :return: vertices_int: number of vertices
    :rtype: integer
    """
    with open(file, "r") as interaction_file:
        vertices_int = len(self.init_dict)  # Counts the number of keys in the dictionary
        return vertices_int

print(count_vertices("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def count_edges(self):
    """This function counts the number of edges of a graph.

    The number of edges is the same than the number of interactions. In the list of interacting proteins, a line
    corresponds to an interaction.

    :param file: interaction file
    :return: edges_int: number of edges
    :rtype: integer
    """
    with open(file, "r") as interaction_file:
        edges_int = len(self.init_list)  # Counts the number of elements in the list
        return edges_int

print(count_edges("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))

def clean_interactome(self):
    """This function reads a file including an protein-protein interaction graph and removes all repeated interactions,
    and all homo-dimers from it.

    It returns a new cleaned up file."""

    with open(file, "r") as interaction_file:

        cleaned_interactome_list = []  # Initialization of the the new interactome list

        for couple_str in self.init_list:  # This loop checks for all elements in the interaction list the following condition
            couple2_str = (couple_str[1], couple_str[0])

            if couple_str not in self.init_list:
                if couple2_str not in self.init_list:  # If an interaction is not repeated
                    if couple_str[0] != couple_str[1]:  # And the first element doesn't equal to the second element (the
                        # protein isn't an homo-dimer)
                        cleaned_interactome_list.append(couple_str)  # Then, the new couple is added to the new list

        cleaned_interactome_file = open("cleaned_interaction_file.txt", "w")
        cleaned_interactome_file.write(str(len(cleaned_interactome_list)) + "\n")  # Writes the number of total interactions in the first line

        for couple_str in cleaned_interactome_list:
            cleaned_interactome_file.write(str(couple_str[0] + " " + couple_str[1] + "\n"))  # Writes the new couples into the new file
        cleaned_interactome_file.close()

    return cleaned_interactome_list

print(clean_interactome("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def get_degree(self, prot):
    """This function takes an interaction graph and a protein name inputed by the user to return the protein's degree in
     the graph.

    :param file: interaction file
    :param prot: inputed protein name
    :type prot: string
    :return: interactions_int: number of the protein's interactions
    :rtype: integer
    """
    with open(file, "r") as interaction_file:
        interactions_int = len(self.init_dict[prot])  # Takes the number of proteins in the dictionary that interacts
        # with the inputed protein

        return interactions_int

print(get_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", "ZYX_HUMAN"))


def get_max_degree(self):
    """This function gets the protein's name with the maximum of interactions and its degree.

    :param file: interaction file
    :return: prot_str, max_degree_int: protein name
    :rtype: string, integer
    """
    with open(file, "r") as interaction_file:
        max_degree_int = 0  # Initialisation of the variable
        prot_str = ""

        for vertice in self.init_dict:  # This loop checks for each element of the dictionary the following condition

            if len(self.init_dict[vertice]) > max_degree_int:  # If the number of proteins that interact with the
                # current one is superior to the current maximum degree
                max_degree_int = len(self.init_dict[vertice])  # Then, the maximum degree takes the new protein degree's value
                prot_str = vertice  # And the argument takes this protein's name

        return (prot_str, max_degree_int)

print(get_max_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def get_ave_degree(self):
    """This function calculates the proteins' average degree of the interaction graph.

    :param file: interaction file
    :return: ave_degree_float: proteins' average degree
    :rtype: float
    """
    with open(file, "r") as interaction_file:

        for vertice in self.init_dict:  # For each protein in the dictionary...
            ave_degree_float = len(self.init_dict[vertice])/len(self.init_dict)  # ... It calculates the number of interaction,
            # And the total number of interactions is divided by the total number of proteins

        return ave_degree_float

print(get_ave_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def count_degree(self, deg):
    """This function calculates the number of proteins from the interaction graph whose degree equals the parameter deg.

    :param file: interaction file
    :param deg: number of interactions inputed by the user
    :type deg: integer
    :return: prot_list: list of proteins that have a degree which equals to the inputed degree
    :rtype: list
    """
    with open(file, "r") as interaction_file:
        prot_list = []

        for vertice in self.init_dict:  # This loop will check each element of the dictionary the following condition
            if len(self.init_dict[vertice]) == deg:  # If the current protein's degree equals to the inputed one...
                prot_list.append(vertice)  # ... Then, the protein is added to the list

        return str(prot_list)

print(count_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", 185))
