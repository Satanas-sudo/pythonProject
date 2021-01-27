from pathlib import Path
import proteinDomains
import pickle
from tqdm import tqdm
import tempfile
import re
import urllib
from collections import Counter


file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")
uniprotFile = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality_uniprot.tab")  # Only proteins from Human_HighQuality.txt in Uniprot
proteome_uniprotFile = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_proteome_uniprot.tab")  # Entire human proteome in Uniprot
uniprotDict_file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/uniprotDic.pkl")

class Interactome:

    # Attributes
    prot_dict = {}
    prot_list = []
    proteins = []
    dict_graph = {}
    dict_graphe = {}
    domains = []

    # Constructor
    def __init__(self, file, uniprotFile, uniprotDict_file):
        self.file = str(file)  # Global variable is stocked with the file name
        self.uniprotFile = str(uniprotFile)
        self.prot_list = self.read_interaction_file_list(file)
        self.prot_dict = self.read_interaction_file_dict(file)
        self.proteins = list(self.prot_dict.keys())
        self.dict_graph = self.xlink_Uniprot(uniprotFile)
        self.dict_graph2 = self.x_link_domains(uniprotFile)
        self.dict_graphe = self.load_uniprotDict(uniprotDict_file)
        self.domains = self.ls_domains()

    # Methods

    def read_interaction_file_dict(self, file):
        """Reads an interaction graph between proteins in a tabulated file and stores it in a dictionary.
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
                    prot_graph_dict[prot_name_str].append((line.split())[1])  # Cuts a string in a list where each word
                    # is an item of the list
                else:
                    prot_name_str = (line.split()[0])  # Storing in memory of the new peak protein (as first element)
                    prot_graph_dict[prot_name_str] = list()  # Creation of the corresponding key
                    prot_graph_dict[prot_name_str].append((line.split())[1])  # Interacting protein addition in the list of
                    # associated value (in second element)
        return prot_graph_dict
    #print(read_interaction_file_dict(file))


    def read_interaction_file_list(self, file):
        """Reads an interaction graph between proteins in a file and stores it in a list of couples.

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
                prot_graph_list.append(line.replace("\t", " \ ", 1).replace("\n", "", 1))
        return prot_graph_list
    #print(read_interaction_file_list(file))


    def read_interaction_file(self, file):
        """Returns a couple where the first element is the dictionary representing the graph and the second
        element is the interaction list representing that same graph.
        This function calls the two prior functions and stores them in a tuple.

        :param file: interaction file
        :return: couple_d_l_tuple: first element is the dictionary and the second element is the list
        :rtype: tuple
        """
        with open(file, "r"):
            couple_d_l_tuple = [(self.read_interaction_file_dict(file)), (self.read_interaction_file_list(file))]
            # Creation of the couple with dictionary as first element and list as second element
        return couple_d_l_tuple
    #print(read_interaction_file(file))


    def count_vertices(self):
        """Counts the number of peaks of a graph.

        :param file: interaction file
        :return: vertices_int: number of vertices
        :rtype: integer
        """
        with open(file, "r"):
            vertices_int = len(self.prot_dict)  # Counts the number of keys in the dictionary
            return vertices_int
    #print(count_vertices(file))


    def count_edges(self):
        """Counts the number of edges of a graph. The number of edges is the same than the number of interactions.
        In the list of interacting proteins, a line corresponds to an interaction.

        :param file: interaction file
        :return: edges_int: number of edges
        :rtype: integer
        """
        with open(file, "r"):
            edges_int = len(self.prot_list)  # Counts the number of elements in the list
            return edges_int


    def clean_interactome(self, file):
        """Reads a file including an protein-protein interaction graph and removes all repeated interactions,
        and all homo-dimers from it. It returns a new cleaned up file.

        :param file: interaction file to clean
        :return: cleaned_interactome_list: all non-redundant interacting proteins
        :rtype: list
        """
        with open(file, "r") as file_to_clean:
            cleaned_interactome_list = []  # Initialization of the the new interactome list
            nb_interaction_int = 0
            file_to_clean.readline()

            for line in file_to_clean:
                couple_list = line.split()
                if tuple(couple_list) and tuple(couple_list[::-1]) not in cleaned_interactome_list:  # Checks if an interaction is not repeated
                    if len(couple_list) == 2:
                        if couple_list[0] != couple_list[1]:  # And if the first element doesn't equal to the second element (the
                            # protein isn't an homo-dimer)
                            nb_interaction_int +=1
                            cleaned_interactome_list.append(tuple(couple_list))  # Then, the new couple is added to the new list
        file_to_clean.close()

        with open("cleaned_interaction_file.txt", "w") as cleaned_interactome_file:
            cleaned_interactome_file.write(str(len(cleaned_interactome_list)) + "\n")  # Writes the number of total interactions in the first line
            for couple_list in cleaned_interactome_list:
                cleaned_interactome_file.write(str(couple_list[0] + " " + couple_list[1] + "\n"))  # Writes the new couples into the new file
        cleaned_interactome_file.close()

        return cleaned_interactome_list


    def get_degree(self, prot):
        """Takes an interaction graph and a protein name inputed by the user to return the protein's degree in the graph.

        :param file: interaction file
        :param prot: inputed protein name
        :type prot: string
        :return: nb_interactions_int: number of the protein's interactions
        :rtype: integer
        """
        with open(file, "r"):
            nb_interactions_int = len(self.prot_dict[prot])  # Takes the number of proteins in the dictionary that interacts
            # with the inputed protein
            return nb_interactions_int
    #get_degree(file, "ZYX_HUMAN")


    def get_max_degree(self):
        """Gets the protein's name with the maximum of interactions and its degree.

        :param file: interaction file
        :return: prot_str, max_degree_int: protein name
        :rtype: string, integer
        """
        with open(file, "r"):
            max_degree_int = 0  # Initialization of the variable
            prot_str = ""

            for vertice in self.prot_dict:  # This loop checks for each element of the dictionary the following condition
                if len(self.prot_dict[vertice]) > max_degree_int:  # If the number of proteins that interact with the
                    # current one is superior to the current maximum degree
                    max_degree_int = len(self.prot_dict[vertice])  # Then, the maximum degree takes the new protein degree's value
                    prot_str = vertice  # And the argument takes this protein's name
            return (prot_str, max_degree_int)


    def get_ave_degree(self):
        """Computes the proteins' average degree of the interaction graph.

        :param file: interaction file
        :return: ave_degree_float: proteins' average degree
        :rtype: float
        """
        with open(file, "r"):
            for vertice in self.prot_dict:  # For each protein in the dictionary...
                ave_degree_float = len(self.prot_dict[vertice])/len(self.prot_dict)  # ... It calculates the number of interaction,
                # And the total number of interactions is divided by the total number of proteins
            return ave_degree_float


    def count_degree(self, deg):
        """Computes the number of proteins from the interaction graph whose degree equals the parameter deg.

        :param file: interaction file
        :param deg: number of interactions inputed by the user
        :type deg: integer
        :return: prot_list: list of proteins that have a degree which equals to the inputed degree
        :rtype: list
        """
        with open(file, "r"):
            prot_list = []
            for vertice in self.prot_dict:  # This loop will check each element of the dictionary the following condition
                if len(self.prot_dict[vertice]) == deg:  # If the current protein's degree equals to the inputed one...
                    prot_list.append(vertice)  # ... Then, the protein is added to the list
            return str(prot_list)
    #print(count_degree(file, 185))


    def histogramm_degree(self, dmin, dmax):
        """Computes, for all degrees between dmin and dmax, the number of proteins having a degree d. It
        returns the result by an histogramm.

        :param dmin: minimum degree of the interval
        :param dmax: maximum degree of the interval
        :return: nb_prot_int: number of proteins within the interval
        """
        with open(file, "r"):
            nb_prot_int = 0

            for vertice_str in self.prot_dict:
                if dmin <= len(self.prot_dict[vertice_str]) <= dmax:  # Verifies that each vertice's degree is within the interval
                    nb_prot_int += 1  # If so, the number of protein is incremented
            deg_int = dmin

            while deg_int <= dmax:  # Histogram created between dmin (here 'degree') and dmax
                deg_int += 1  # While the different degrees in the interval are explored,
                hist_str = str(deg_int) + " "  # Each line is a string composed of the degree + 'space' + 'n*'
                for vertice_str in self.prot_dict:  # Each vertice in the dictionnary is checked
                    if len(self.prot_dict[vertice_str]) == deg_int:  # Each time a vertice has the same degree as the current degree,
                        hist_str += "*"  # the line of the current degree gains one '*'
                print(hist_str)

            return nb_prot_int
    # We observe that most proteins interact with few proteins

########## RELATED COMPONENTS ##########

    def extractCC(self, prot):
        """Computes the related components number of a graph, and gives for each of them their size (its
        number of proteins).

        :param prot: inputed protein
        :return: CC_list: list of proteins in the same connected component than the protein inputed
        """
        CC_list = [prot]
        dict_list = self.prot_dict.get(prot)

        for prot in dict_list:
            if prot not in CC_list:
                CC_list.append(prot)  # Protein neighbor is added to the list
        return CC_list
    #print(file.extractCC('C'))


    def countCC(self):
        """Computes the related components number of a graph, and gives for each of them their size (its
        number of proteins).

        :return: CC_list: list of connected componants and their size
        """
        prot_list = list(self.prot_dict.keys())
        CC_list = []
        number_CC_int = 1  # Counter initialization

        while len(prot_list) != 0:
            lcc = self.extractCC(prot_list[0])
            for element in lcc[:]:
                prot_list.remove(prot_list[0])  # Proteins already in a related component are removed,
                # port_list[0] is the first element (first protein, str type) implicated in the connected componant list
            CC_list.append((str(number_CC_int), len(lcc)))
            number_CC_int += 1
            print("The related component number ", str(number_CC_int), "is composed of ", str(len(element)), "proteins")
        return CC_list


    def computeCC(self):
        """Returns a lcc list where each lcc[i] element corresponds to the related component number of the
         protein at a i location in the proteins list of the graph (which is a class attribute).

        :return: lcc_list: list of related componants
        """
        lcc_list = []  # Initialization of the list containing the related component number for each protein
        for prot in self.proteins:
            lcc_list.append(self.extractCC(prot))
        return lcc_list


    def writeCC(self):
        """Writes in a file the different related components of a graph. The format follows these rules :
        1. One line by related component.
        2. The first element of the line is the size of the related component.
        3. Then, the list of peaks that compose this related component will be added.

        :return: None
        """
        lcc = self.countCC()
        prot_list = list(self.prot_dict.keys())
        CC_file = open('CC_file.txt', 'w')

        for tuple_iterator in lcc:
            CC_name_str = tuple_iterator[0]
            len_int = tuple_iterator[1]
            CC = self.extractCC(prot_list[0])
            CC_file.write(str(CC_name_str) + "\t" + str(len_int) + "\t" + str(CC) + "\n")
        CC_file.close()
        return


    def density(self):
        """Returns the density of the (non-oriented) graph. The density is the number of present edges
        divided by the number of total edges which could theoretically be.

        :return: D_float: density of the graph
        :rtype: float
        """
        D_float = 2*int(self.count_edges())/(int(self.count_vertices())*(int(self.count_vertices())-1))
        return D_float


    def clustering(self, prot):
        """Returns the clustering coefficient of the "prot" summit in the graph. The clustering coefficient
        is a measure of the nodes grouping in a network. This coefficient is the probability that two nodes are connected
        knowing they have a common neighbor. The local clustering coefficient of a node i is the fraction of its connected
        neighbors pairs, equal to 0 is di ≤ 1 by convention.

        :param prot: protein inputed
        :return: clustering coefficient
        :rtype: float
        """
        triangle_int = 0
        #if len(self.prot_dict[prot]) != None:
        for first_iter in self.prot_dict[prot]:
            for sec_iter in self.prot_dict[first_iter]:
                if first_iter in self.prot_dict[sec_iter]:
                    triangle_int += 1
            pair_neighbor_float = (len(self.prot_dict[prot]) * (len(self.prot_dict[prot]) - 1)) / 2
        return float(triangle_int/pair_neighbor_float)

########## PROTEIN DOMAINS ##########

    def xlink_Uniprot(self, uniprotFile):
        """Extends, for each protein of the interactome, the dictionary stocking the interactome.

        :param uniprotFile: interaction file
        :return: graph_dict: dictionary with UniprotID linked with the protein file path
        """
        graph_dict = {}
        Uniprot_dict = proteinDomains.Proteome("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality_uniprot.txt")

        for prot in self.prot_dict.keys():
            if prot in Uniprot_dict.ID_name_dict.keys():
                graph_dict[prot] = {}
                graph_dict[prot]["UniprotID"] = Uniprot_dict.ID_name_dict[prot]
                graph_dict[prot]["Neighbor"] = self.prot_dict[prot]
        return graph_dict


    def get_protein_domains_list(self, prot):
        """Allows to search for the Uniprot page corresponding to a protein entry and returns the domain(s)'s
        name(s) contained by that protein. urllib, tempfile (creates temporary files) and re are needed.

        :param prot: inputed protein
        :return: domain_names_list: list of domains associated with the inputed protein
        """
        url_str = "https://www.uniprot.org/uniprot/{}.txt".format(prot)

        with urllib.request.urlopen(url_str) as rep_request:
            with tempfile.TemporaryFile() as tmp_file:  # Creation of the temporary file which contains the request result
                tmp_file.write(rep_request.read())
                tmp_file.seek(0)  # All documents are found (argument 0)
                # We search for the lines that have this regular expression
                # Use of re to find the lines beginning by "DR   Pfam" and stocks them into a list
                result_domains = re.findall(r'DR\s{3}Pfam;\s\w+;\s\w+', str(tmp_file.read()))
                # 4th elements of each lines (domain name) are extracted and returned as a list
                domain_names_list = []

                for element_str in result_domains:
                    domain_names_list.append(element_str.split()[3])  # The domain is the 3rd element of each line found
                return domain_names_list


    def x_link_domains(self, uniprotFile):
        """Allows for each protein of the interactome to extend our dictionary data structure to add the proteins' domains.

        :param uniprotFile: interaction file
        :return: dictionary
        """
        with open("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/uniprotDict.pkl", 'wb') as file:  # Creation of the uniprot file and dict_graph is stored into it
            for prot in tqdm(self.dict_graph.keys()):
                # Search the protein domain for each protein in Pfam and append it in the dictionary
                self.dict_graph[prot]["domains"] = self.get_protein_domains_list(prot)
                #self.dict_graph[prot]["domains"] = proteinDomains.Proteome.get_protein_domains_list(prot)
            pickle.dump(self.dict_graph, file)
        file.close()
        return self.dict_graph

########## DOMAINS COMPOSITION ANALYSIS ##########

    def load_uniprotDict(self, uniprotDict_file):
        """Initializes the dictionary with the Uniprot domains and PFAM keys.

        :param uniprotDict_file: pickle file
        :return: dict
        """
        with open(uniprotDict_file, 'rb') as f:
            dict = pickle.load(f)
        f.close()
        return dict


    def ls_proteins(self):
        """Returns the non-redundant list of all proteins found in the protein-protein interaction graph.
        :return: list of non-redundant proteins
        """
        return self.proteins


    def ls_domains(self):
        """Returns the non-redundant list of all domains found in the proteins composing the interaction graph.

        :return: domain_list: list of all the non-redundant domains
        """
        with open(uniprotDict_file, 'rb') as domain_dict:
             domain_list = []
             proteome_dict = pickle.load(domain_dict)  # Loads the global dictionnary with name_prot, UniprotID, neighbors and domains

             for prot in proteome_dict.keys():  # Searches for each protein
                 domains = proteome_dict[prot]["domains"]
                 for dom in domains:
                     if dom not in domain_list:
                         domain_list.append(dom)  # Adds the domain to the domains' list
                     if not domains:
                        del domain_list[domain_list.index(domains)]  # Removing of empty domains
        return domain_list


    def ls_domains_n(self, n):
        """Returns the non-redundant list of all domains found at least n times in the proteins composing
        the interaction graph.

        :param n: inputed number of least times a domain appears
        :return: nDomain_list: list of domains found
        """
        dict_graphe_values_list = list(self.dict_graphe.values())
        domain_list = []
        nDomain_list = []
        # Recovery of all domains from the nested interaction dictionary
        for i in dict_graphe_values_list:
            for dom in i.get('domains'):
                domain_list.append(dom)

        # Counts all proteins and sort in a sub-class type dictionary
        count_dict = Counter(sorted(domain_list))
        for dom in count_dict.keys():
            if count_dict[dom] >= n:
                nDomain_list.append(dom)
        return nDomain_list


    def co_occurence(self, dom_x, dom_y):
        """Computes and returns the co-occurrences of x and y domains number in the proteins of the interactome.

        :param dom_x: inputed domain
        :param dom_y: inputed domain
        :return: count_occ_int: number of co-occurrences where the domains x and y appear
        """
        count_occ_int = int(0)
        for i in self.dict_graphe.values():
            if dom_x in i.get('domains') and dom_y in i.get('domains'):
                count_occ_int += 1
        return count_occ_int


### SCRIPT TEST ###

#print(Interactome(file, uniprotFile, uniprotDict_file).read_interaction_file_dict(file))
#print(Interactome(file, uniprotFile, uniprotDict_file).read_interaction_file_list(file))
#print(Interactome(file, uniprotFile, uniprotDict_file).read_interaction_file(file))
#print(Interactome(file, uniprotFile, uniprotDict_file).count_vertices())  # 6660
#print(Interactome(file, uniprotFile, uniprotDict_file).count_edges())  # 27276
#print(Interactome(file, uniprotFile, uniprotDict_file).clean_interactome(file))
#print(Interactome(file, uniprotFile, uniprotDict_file).get_degree("ZYX_HUMAN"))  # 2
#print(Interactome(file, uniprotFile, uniprotDict_file).get_max_degree())  # ('ATX1_HUMAN', 185)
#print(Interactome(file, uniprotFile, uniprotDict_file).get_ave_degree())  # 0.0003003003003003003
#print(Interactome(file, uniprotFile, uniprotDict_file).count_degree(185))  # ['ATX1_HUMAN']
#print(Interactome(file, uniprotFile, uniprotDict_file).histogramm_degree(1, 100))

#print(Interactome(file, uniprotFile, uniprotDict_file).extractCC('ZYX_HUMAN'))
#print(Interactome(file, uniprotFile, uniprotDict_file).countCC())
#print(Interactome(file, uniprotFile, uniprotDict_file).computeCC())
#print(Interactome(file, uniprotFile, uniprotDict_file).writeCC())
#print(Interactome(file, uniprotFile, uniprotDict_file).density())  # 0.0012300632213532049
#print(Interactome(file, uniprotFile, uniprotDict_file).clustering('ATX1_HUMAN'))

#print(Interactome(file, uniprotFile, uniprotDict_file).xlink_Uniprot(uniprotFile))
#print(Interactome(file, uniprotFile, uniprotDict_file).xlink_domains(uniprotFile))
#print(Interactome(file, uniprotFile, uniprotDict_file).load_uniprotDict(uniprotDict_file))
#print(Interactome(file, uniprotFile, uniprotDict_file).ls_proteins())
#print(Interactome(file, uniprotFile, uniprotDict_file).ls_domains())
#print(Interactome(file, uniprotFile, uniprotDict_file).ls_domains_n(3))
#print(Interactome(file, uniprotFile, uniprotDict_file).co_occurence('NEBL_HUMAN','ZYX_HUMAN'))