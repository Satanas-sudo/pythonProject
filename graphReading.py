import os
import numpy as np  # needed for the array


def read_interaction_file_dict(Human_HighQuality):
    """This function reads an interaction graph between proteins in a tabulated file and stores it in a dictionary.

    "Peaks" proteins are considered as keys and proteins that interacts with them are considered as the values
    associated to the keys.

    :param Human_HighQuality: interaction file
    :return: prot_graph_dict: first element as key and their neighbor(s) in second element as the key's value(s)
    :rtype: dictionary
    """
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())  # File first line escapement

        data_str = interaction_file.readlines()  # Storing in memory the file's data
        prot_graph_dict = {}  # Dictionary declaration
        prot_name_str = ""  # Initialization of peak protein's name

        # Loop for each line of the document in memory (this loops needs that the file's data be classified)
        for line in data_str:

            # Checking for a key that already exists for the peak protein
            if prot_name_str == (line.split()[0]):
                # Interacting protein addition in the list of associated values to the peak protein
                prot_graph_dict[prot_name_str].append((line.split())[1])  # Cuts a string in a list where each word is
                # an item of the list

            else:
                prot_name_str = (line.split()[0])  # Storing in memory of the new peak protein (as first element)
                prot_graph_dict[prot_name_str] = list()  # Creation of the corresponding key
                prot_graph_dict[prot_name_str].append((line.split())[1])  # Interacting protein addition in the list of
                # associated value (in second element)

    return prot_graph_dict

print(read_interaction_file_dict("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                 "Human_HighQuality.txt"))


def read_interaction_file_list(Human_HighQuality):
    """This function reads an interaction graph between proteins in a file and stores it in a list of couples.

    :param Human_HighQuality: interaction file
    :return: prot_graph_list: couples of interacting proteins
    :rtype: list
    """
    with open(Human_HighQuality, "r") as interaction_file:

        lines_num_int = int(interaction_file.readline())

        data_str = interaction_file.readlines()
        prot_graph_list = []  # List creation

        for line in data_str:  # For each line of the interaction file (and so proteins couple), replace tabulations by
            # a backslash and adds them in the list
            prot_graph_list.append(line.replace("\t", "\ ", 1).replace("\n", "", 1))

    return prot_graph_list

print(read_interaction_file_list("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                 "Human_HighQuality.txt"))


def read_interaction_file_mat(Human_HighQuality):
    """This function reads an interaction graph between proteins in a file and stores it in an adjacency matrix.

    :param Human_HighQuality: interaction file
    :return: prot_interaction_mat: matrix of all interacting proteins that are associated with each other
    :rtype: array
    """
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())
        data_str = interaction_file.readlines()
        prot_data_list = []
        prot_graph_dict = {}

        for line in data_str:
            prot1_exist = False
            prot2_exist = False

            # This loop includes elements in a list, pass only once (when the list is empty), initialize the first
            # protein
            if prot_data_list == []:  # Checks if list is empty (not mandatory though)
                prot_data_list = [(line.split()[0])]  # Adds the value into the list
                prot_graph_dict[(line.split()[0])] = list()  # Declaration of the key without value (value is a list)
                prot_graph_dict[(line.split()[0])].append((line.split())[1])  # Adds the value as second element

            # When the function will screen the list, it will match twice only instead of passing at each loops (not
            # mandatory but it processes faster than passing every loops)
            for elt in prot_data_list:
                if line.split()[0] == elt or line.split()[1] == elt:
                    if (line.split()[0]) == elt:
                        prot1_exist = True
                        prot_graph_dict[(line.split()[0])].append((line.split())[1])
                    if (line.split()[1]) == elt:
                        prot2_exist = True
                        prot_graph_dict[(line.split()[1])].append((line.split())[0])

            if prot2_exist == False:
                prot_data_list.append((line.split()[1]))
                prot_graph_dict[(line.split()[1])] = list()
                prot_graph_dict[(line.split()[1])].append((line.split())[0])

            if prot1_exist == False:
                prot_data_list.append((line.split()[0]))
                prot_graph_dict[(line.split()[0])] = list()
                prot_graph_dict[(line.split()[0])].append((line.split())[1])

        prot_data_list.sort()  # Sorts in alphabetic order in order to find the proteins' order (the array returns the
        # number of proteins positioned as the second element only.
        lines_for_mat_dict = {}  # Dictionary initialization
        prot_interaction_mat = np.empty((len(prot_data_list), len(prot_data_list)), int)  # Creation of an empty array
        # which defines dimensions based on proteins lists length (number of columns, number of lines)

        # This loops compares the protein list from the top to the bottom. It shifts at each case of one dimension, a
        # key is declared in the dictionary for each dimension. A key is declared at each loop. At the end of the loop,
        # we retrieve each key and all values associated to a dimension.
        for i, elt in enumerate(prot_data_list):
            lines_for_mat_dict[i] = list()  # Dictionary declaration which will be used to create the array, it browses
            # all the proteins list
            print("entry in comparison", elt)  # From 0 to 27276, declares the dictionary key

            for j, el in enumerate(prot_data_list):  # For each line, browses all columns, from left to right
                compare_list = [prot_graph_dict[el]]  # Gives a comparison by list from protein in the dictionary of
                # interaction
                for k, e in enumerate(compare_list):  # compares from top to bottom, left to right by shifting case by
                    # case
                    if el == e:
                        compare_boo = True  # For each case, tests if there is a correspondence between the protein of
                        # the column and the protein of the line
                        break  # If a correspondence is found, the loop breaks (prevents from loop running until the end
                        # of the list, so if it founds a corresponding protein and it is not the last one, then it will
                        # return false and won't count it
                    else:
                        compare_boo = False  # Returns true if there is an interaction or false if there isn't
                if compare_boo == True:  # Key created before is inserted in the array
                    lines_for_mat_dict[i].append(1)  # Adds the correspondence to the right column if it returns true
                else:
                    lines_for_mat_dict[i].append(0)  # Adds the protein to the left column if it returns false

            # prot_interact_mat is the variable name, i is the dimension where i want to add it,
            # [lines_for_mat_dict[i]] is the dictionary line with value i, axis=0 is the direction where we want to add i
            np.insert(prot_interaction_mat, i, [lines_for_mat_dict[i]], axis=0)

    return prot_interaction_mat

print(read_interaction_file_mat("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                "Human_HighQuality.txt"))


def read_interaction_file(Human_HighQuality):
    """This function returns a couple where the first element is the dictionary representing the graph and the second
    element is the interaction list representing that same graph.

    This function calls the two prior functions and stores them in a tuple.

    :param Human_HighQuality: interaction file
    :return: couple_d_l_tuple: first element is the dictionary and the second element is the list
    :rtype: tuple
    """
    with open(Human_HighQuality, "r") as interaction_file:

        d_int = read_interaction_file_dict(Human_HighQuality)  # Prior function which created the dictionary is stored
        # in variable d_int
        l_int = read_interaction_file_list(Human_HighQuality)  # Prior function which created the list is stored in
        # variable l_int
        couple_d_l_tuple = [d_int, l_int]  # Creation of the couple with dictionary as first element and list as second
        # element
    return couple_d_l_tuple

print(read_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                            "Human_HighQuality.txt"))
