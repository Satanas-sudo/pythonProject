import os
from graphReading import read_interaction_file_dict
from graphReading import read_interaction_file_list


def count_vertices(Human_HighQuality):
    """This function counts the number of peaks of a graph.

    :param Human_HighQuality: interaction file
    :return: vertices_int: number of vertices
    :rtype: integer
    """
    with open(Human_HighQuality, "r") as interaction_file:
        vertices_int = len(read_interaction_file_dict(Human_HighQuality))  # Counts the number of keys in the dictionary
        return vertices_int

print(count_vertices("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def count_edges(Human_HighQuality):
    """This function counts the number of edges of a graph.

    The number of edges is the same than the number of interactions. In the list of interacting proteins, a line
    corresponds to an interaction.

    :param Human_HighQuality: interaction file
    :return: edges_int: number of edges
    :rtype: integer
    """
    with open(Human_HighQuality, "r") as interaction_file:
        edges_int = len(read_interaction_file_list(Human_HighQuality))  # Counts the number of elements in the list
        return edges_int

print(count_edges("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))

def clean_interactome(Human_HighQuality):
    """This function reads a file including an protein-protein interaction graph and removes all repeated interactions,
    and all homo-dimers from it.

    It returns a new cleaned up file."""

    with open(Human_HighQuality, "r") as interaction_file:

        all_interactions_list = read_interaction_file_list(Human_HighQuality)
        cleaned_interactome_list = []  # Initialization of the the new interactome list

        for couple_str in all_interactions_list:  # This loop checks for all elements in the interaction list the following condition
            couple2_str = (couple_str[1], couple_str[0])

            if couple_str not in all_interactions_list:
                if couple2_str not in all_interactions_list:  # If an interaction is not repeated
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


def get_degree(Human_HighQuality, prot):
    """This function takes an interaction graph and a protein name inputed by the user to return the protein's degree in
     the graph.

    :param Human_HighQuality: interaction file
    :param prot: inputed protein name
    :type prot: string
    :return: interactions_int: number of the protein's interactions
    :rtype: integer
    """
    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        interactions_int = len(prot_graph_dict[prot])  # Takes the number of proteins in the dictionary that interacts
        # with the inputed protein

        return interactions_int

print(get_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", "ZYX_HUMAN"))


def get_max_degree(Human_HighQuality):
    """This function gets the protein's name with the maximum of interactions and its degree.

    :param Human_HighQuality: interaction file
    :return: prot_str, max_degree_int: protein name
    :rtype: string, integer
    """
    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        max_degree_int = 0  # Initialisation of the variable
        prot_str = ""

        for vertice in prot_graph_dict:  # This loop checks for each element of the dictionary the following condition

            if len(prot_graph_dict[vertice]) > max_degree_int:  # If the number of proteins that interact with the
                # current one is superior to the current maximum degree
                max_degree_int = len(prot_graph_dict[vertice])  # Then, the maximum degree takes the new protein degree's value
                prot_str = vertice  # And the argument takes this protein's name

        return (prot_str, max_degree_int)

print(get_max_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def get_ave_degree(Human_HighQuality):
    """This function calculates the proteins' average degree of the interaction graph.

    :param Human_HighQuality: interaction file
    :return: ave_degree_float: proteins' average degree
    :rtype: float
    """
    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)

        for vertice in prot_graph_dict:  # For each protein in the dictionary...
            ave_degree_float = len(prot_graph_dict[vertice])/len(prot_graph_dict)  # ... It calculates the number of interaction,
            # And the total number of interactions is divided by the total number of proteins

        return ave_degree_float

print(get_ave_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def count_degree(Human_HighQuality, deg):
    """This function calculates the number of proteins from the interaction graph whose degree equals the parameter deg.

    :param Human_HighQuality: interaction file
    :param deg: number of interactions inputed by the user
    :type deg: integer
    :return: prot_list: list of proteins that have a degree which equals to the inputed degree
    :rtype: list
    """
    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        prot_list = []

        for vertice in prot_graph_dict:  # This loop will check each element of the dictionary the following condition
            if len(prot_graph_dict[vertice]) == deg:  # If the current protein's degree equals to the inputed one...
                prot_list.append(vertice)  # ... Then, the protein is added to the list

        return str(prot_list)

print(count_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", 185))
