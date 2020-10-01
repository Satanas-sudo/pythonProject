import os
import numpy as np  # needed for the array


def read_interaction_file_dict(Human_HighQuality):
    """This function reads an interaction graph between proteins in a tabulated file and stores it in a dictionnary.

    "Peaks" proteins are considered as keys and proteins that interacts with them are considered as the values
    associated to the keys.

    :param Human_HighQuality: interaction file
    :return: prot_graph_dict
    :rtype: dictionnary
    """
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())  # File first line escapement

        data_str = interaction_file.readlines()  # Storing in memory the file's data
        prot_graph_dict = {}  # Dictionnary declaration
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
                prot_graph_dict[prot_name_str] = list()  # Corresponding key creation
                prot_graph_dict[prot_name_str].append((line.split())[1])  # Interacting protein addition in the list of
                # associated value (in second element)

    return prot_graph_dict

print(read_interaction_file_dict("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                 "Human_HighQuality.txt"))


def read_interaction_file_list(Human_HighQuality):
    """This function reads an interaction graph between proteins in a file and stores it in a list of couples.

    :param Human_HighQuality: interaction file
    :return: prot_graph_list
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
    with open(Human_HighQuality, "r") as interaction_file:
        lines_num_int = int(interaction_file.readline())
        data_str = interaction_file.readlines()
        prot_data_list = []
        prot_graph_dict = {}
        for line in data_str:
            prot1_exist = False
            prot2_exist = False
            # Cette boucle insère des élements dans une liste, passe qu'une fois (quand la liste est vide), initialise la première prot
            if prot_data_list == []:  # Vérifie que la liste est vide mais on pourrait s'en passer puisque la ligne suivant met direct la valeur, les ajoute
                prot_data_list = [(line.split()[0])]  #Piste d'amélioration de la boucle -> liste de compréhension
                prot_graph_dict[(line.split()[0])] = list() # Déclaratin de la clé sans valeur. Déf d'une val du dict, cette val associée à clé est une liste
                prot_graph_dict[(line.split()[0])].append((line.split())[1])  # Ajout de la valeur en 2e élément
            # Quand la fct fera le tour de la liste, elle va matcher que 2 fois au lieu passer à chaque toutes les boucles
            #
            for elt in prot_data_list:  # piste d'amélioration -> remplacer les for par compréhension de listes
                if line.split()[0] == elt or line.split()[1] == elt:  # pas obligatoire, ça va juste plus vite plutôt que de passer toutes les boucles
                    if (line.split()[0]) == elt:  #6'
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
        prot_data_list.sort()  # trie dans l'ordre alphabétique, permet de retrouver l'ordre des prot (la matrice retourne que le num des prot de droite)
        lines_for_mat_dict = {}
        prot_interaction_mat = np.empty((len(prot_data_list), len(prot_data_list)), int)  # création d'un array vide, déf des dim basée sur la lg des listes des prot (nb colo, nb lignes)
        for i, elt in enumerate(prot_data_list):  # compare liste de prot du haut vers le bas (décale à chq fois d'une dim, pr chq dim on déclare une clé ds dico, à la fin je récupe chaque clé et ttes les val asso corresp à une dim) (décla d'une clé à chq tor de boucle for)
            lines_for_mat_dict[i] = list() # déclaration du dico qui va servir à créa de la matrice, parcourt tte la liste de prot
            print("entrée en comparaison", elt) # de 0 à 27276, déclare la clé '103" du dico (corresp à 1 dim)
            for j, el in enumerate(prot_data_list):  # pour chaque ligne, parcourt ttes les colo de dim/ttes les listes/la matrice de gauche à droite (corresp aussi à liste)
                compare_list = [prot_graph_dict[el]] # sort une comp par liste à partir de la prot dans le dico d'inter # le plus grd nb d'inter est 10 ou 11 donc y'a juste ça à parcourir pr les comp
                for k, e in enumerate(compare_list):  # compare à chaque fois, on se déplace de case en case, gauche à d, h en b
                    if el == e:
                        compare_boo = True  # pour chaque case, test si la case où y s'trouve y'a une corresp entre prot de la colo et prot de la ligne
                        break  # si corresp trouvée, on casse la boucle sinon ça continue de tourner jusqu'à fin de liste (dc si prot corresp mais c'est pas la dernière il sort false et il la compte pas)
                    else:
                        compare_boo = False  # t ou f si y'a interaction ou pas, c'est pcq parcourt 2 fois tte la liste avec un test
                if compare_boo == True: # la clé créée ligne 118 est insérée dans l'array
                    lines_for_mat_dict[i].append(1) # si c'est true, y'a corresp donc ça l'ajoute à 1
                else:
                    lines_for_mat_dict[i].append(0) # sinon à 0
            np.insert(prot_interaction_mat, i, [lines_for_mat_dict[i]], axis=0)  # nom variable/tab, la dim où j'veux l'ajouter i, [valeur que j'veux ajouter donc ici la ligne du dico avc val i pcq c'est tjrs la mm ligne, dans quel sens on veut l'ajouter
    return prot_interaction_mat

print(read_interaction_file_mat("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt"))


def read_interaction_file(Human_HighQuality):
    """This function returns a couple where the first element is the dictionnary representing the graph and the second
    element is the interaction list representing that same graph.

    This function calls the two prior functions and stores them in a tuple.

    :param Human_HighQuality: interaction file
    :return: couple_d_l_tuple
    :rtype: tuple
    """
    with open(Human_HighQuality, "r") as interaction_file:

        couple_d_l_tuple = [(read_interaction_file_dict(Human_HighQuality)), (read_interaction_file_list(Human_HighQuality))]

    return couple_d_l_tuple

print(read_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                            "Human_HighQuality.txt"))
