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

"""def clean_interactome(Human_HighQuality):
    This function reads a file including an protein-protein interaction graph and removes all repeated interactions,
    and all homo-dimers from it.

    It returns a new cleaned up file.

    with open(Human_HighQuality, "r") as interaction_file:

        all_interactions_list = read_interaction_file_list(Human_HighQuality)
        cleaned_interactome_list = []

        for couple_str in all_interactions_list:
            couple2_str = (couple_str[1], couple_str[0])

            if couple_str not in all_interactions_list:
                if couple2_str not in all_interactions_list:
                    if couple_str[0] != couple_str[1]:
                        cleaned_interactome_list.append(couple_str)

        cleaned_interactome_file = open("cleaned_interaction_file.txt", "w")
        for couple_str in cleaned_interactome_list:
            cleaned_interactome_file.write(str(couple_str[0] + " " + couple_str[1] + "\n"))

        cleaned_interactome_file.close()

    return cleaned_interactome_list

print(clean_interactome("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))"""


def get_degree(Human_HighQuality, prot):

    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        interactions_int = len(prot_graph_dict[prot])

        return interactions_int

print(get_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", "ZYX_HUMAN"))


def get_max_degree(Human_HighQuality):

    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        max_degree_int = 0
        prot_str = ""

        for vertice in prot_graph_dict:

            if len(prot_graph_dict[vertice]) > max_degree_int:
                max_degree_int = len(prot_graph_dict[vertice])
                prot_str = vertice

        return (prot_str, max_degree_int)

print(get_max_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def get_ave_degree(Human_HighQuality):

    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)

        for vertice in prot_graph_dict:
            ave_degree_float = len(prot_graph_dict[vertice])/len(prot_graph_dict)

        return ave_degree_float

print(get_ave_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def count_degree(Human_HighQuality, deg):

    with open(Human_HighQuality, "r") as interaction_file:
        prot_graph_dict = read_interaction_file_dict(Human_HighQuality)
        prot_list = []

        for vertice in prot_graph_dict:
            if len(prot_graph_dict[vertice]) == deg:
                prot_list.append(vertice)

        return str(prot_list)

print(count_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", 185))


"""def histogram_degree(Human_HighQuality, dmin, dmax):

print(count_degree("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt", 30))"""