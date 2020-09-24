import os

# DICTIONNAIRE
def read_interaction_file_dict(Human_HighQuality):
    """résumé

    contenu détaillé

    :param Human_HighQuality:
    :type Human_HighQuality:
    :return:
    :rtype:
    """
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())  # Echappement de la première ligne du fichier

        # Mise en mémoire des données du fichier, déclaration du dictionnaire, initialisation de la variable du nom de la protéine du sommet
        data_str = interaction_file.readlines()
        prot_graph_dict = {}
        prot_name_str = ""

        # Boucle pour chaque ligne du document en mémoire (cette boucle nécessite que les données du fichier soient classées)
        for line in data_str:

            # Vérification d'une clé déja existante pour la protéine "sommet"
            if prot_name_str == (line.split()[0]):
                # Ajout de la protéine d'interaction dans la liste de valeurs associée a la protéine sommet
                prot_graph_dict[prot_name_str].append((line.split())[1])  # Découpe un string en une liste où chaque mot est un item de la liste

            else:
                # Mise en mémoire de la nouvelle protéine sommet (en premier élément), création de la clé correspondante, ajout de la protéine d'interaction dans la liste de valeur associée (en second élément)
                prot_name_str = (line.split()[0])
                prot_graph_dict[prot_name_str] = list()
                prot_graph_dict[prot_name_str].append((line.split())[1])

    return prot_graph_dict

print(read_interaction_file_dict("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))

# LISTE DE COUPLES
def read_interaction_file_list(Human_HighQuality):

    with open(Human_HighQuality, "r") as interaction_file:

        lines_num_int = int(interaction_file.readline())

        data_str = interaction_file.readlines()
        prot_graph_list = []  # Création de la liste

        for line in data_str:  # Pour chaque ligne du fichier d'interaction (et donc couple de protéines), remplace les tabulations par un backslash et les ajoute dans la liste

            prot_graph_list.append(line.replace("\t", "\ ", 1).replace("\n", "", 1))

    return prot_graph_list

print(read_interaction_file_list("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))

#MATRICE D'ADJACENCE
"""import numpy as np

def read_interaction_file_mat(Human_HighQuality):
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())
        data_str = interaction_file.readlines()
        prot_data_list = []
        prot_graph_dict = {}
        for line in data_str:
            prot1_exist = False
            prot2_exist = False
            if prot_data_list == []:
                prot_data_list = [(line.split()[0])]
                prot_graph_dict[(line.split()[0])] = list()
                prot_graph_dict[(line.split()[0])].append((line.split())[1])
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
        prot_data_list.sort()
        lines_for_mat_dict = {}
        prot_interaction_mat = np.empty((len(prot_data_list), len(prot_data_list)), int) #array
        for i, elt in enumerate(prot_data_list):
            lines_for_mat_dict[i] = list()
            print("entrée en comparaison", elt)
            for j, el in enumerate(prot_data_list):
                compare_list = [prot_graph_dict[el]]
                for k, e in enumerate(compare_list):
                    if el == e:
                        compare_boo = True
                        break
                    else:
                        compare_boo = False
                if compare_boo == True:
                    lines_for_mat_dict[i].append(1)
                else:
                    lines_for_mat_dict[i].append(0)
            np.insert(prot_interaction_mat, i, [lines_for_mat_dict[i]], axis=0)
    return prot_interaction_mat

print(read_interaction_file_mat("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))"""

# TUPLE DE COUPLES
def read_interaction_file(Human_HighQuality):

    with open(Human_HighQuality, "r") as interaction_file:

        couple_d_l_tuple = [(read_interaction_file_dict(Human_HighQuality)), (read_interaction_file_list(Human_HighQuality))]

    return couple_d_l_tuple

print(read_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))

# TEST FORMAT FICHIER
def is_interaction_file(Human_HighQuality):

    with open(Human_HighQuality, "r") as interaction_file:

        if

read_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")