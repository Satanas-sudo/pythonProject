import numpy as np
import proteinDomains
import pickle
#import interactomeExplore
from collections import Counter
from matplotlib import pyplot as plt
from tqdm import tqdm
import networkx as nx
import csv
from pathlib import Path


uniprotDict_file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/uniprotDic.pkl")
cooccurenceDict_file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/cooccurenceDic.pkl")
cooccurenceDict_file_withoutNb = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/coocurrenceDic_without_nb_occurence.pickle")

class DomainGraph:

    # Attributes
    dict_graph = {}
    co_occurrence_graph = {}
    edges = int
    domains = []

    # Constructor
    def __init__(self, uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb):
        self.dict_graph = self.load_uniprotDict(uniprotDict_file)
        self.co_occurrence_graph = self.load_co_occurenceDict(cooccurenceDict_file)
        #self.edges = self.count_edges()
        self.nb_vertices = len(self.domains)
        self.domains = self.ls_domains()

    # Methods
    def load_uniprotDict(self, uniprotDict_file):
        """Initializes the dictionary with the Uniprot domains and PFAM keys.

        :param uniprotDict_file: pickle file
        :return: dict
        """
        with open(uniprotDict_file, 'rb') as file:
            dict = pickle.load(file)
        file.close()
        return dict


    def load_co_occurenceDict(self, cooccurenceDict_file):
        """Loads the co-occurrence domain dictionary created with generate_weighted_cooccurrence_graph() where keys are
        the domain's name and the value associated is every domain present with the key in a minimum of one protein.
        """
        with open(cooccurenceDict_file, 'rb') as file:
            dict = pickle.load(file)
        file.close()
        return dict


    def load_cooccurenceDict_withoutNb(self, cooccurenceDict_file_withoutNb):
        with open(cooccurenceDict_file_withoutNb, 'rb') as file:
            dict = pickle.load(file)
        return dict


    def generate_cooccurence_graph(self):
        """Creates an instance and returns a new DomainGraph object of which the peaks are the domains. There must be an interaction
        between the peaks x and y only when these two domains are co-occurrent in at least one protein of the graph.
        """
        cooccurrence_dict = {}
        for prot in self.dict_graph.values():
            domains_list = list(prot.get('domains'))  # The domains are stored in a list
            if len(domains_list) > 1:  # Proteins that have more than one domain
                for domain1 in domains_list:
                    for domain2 in domains_list:
                        # Verifies if different domain and protein have more 1 domain
                        if domain2 != domain1 and len(domain1) > 0 and len(domain2) > 0:
                            if domain1 not in cooccurrence_dict: cooccurrence_dict[domain1] = {domain2: 0}
                            if domain2 not in cooccurrence_dict: cooccurrence_dict[domain2] = {domain1: 0}
                            if domain2 not in cooccurrence_dict[domain1]: cooccurrence_dict[domain1] = {domain2: 0}
                            if domain1 not in cooccurrence_dict[domain2]: cooccurrence_dict[domain2] = {domain1: 0}
                            # If the two domains are present, the counter increases
                            cooccurrence_dict[domain1][domain2] += 1
                            cooccurrence_dict[domain2][domain1] += 1
        # Returns the dictionary created in a file
        print(cooccurrence_dict)
        pickle.dump(cooccurrence_dict, open("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/cooccurenceDic.pkl", 'wb'))
        #return("Number of co-occurences between domains : %s" %(sum(len(v) for v in iter(cooccurrence_dict.values()))))
        # Number of co-occurences between domains : 2459


    def density(self):
        """Returns the density of the (non-oriented) graph. The density is the number of present edges
        divided by the number of total edges which could theoretically be."""
        #D_float = 2*int(self.count_edges())/(int(self.nb_vertices())*(int(self.nb_vertices())-1))
        #return D_float

        count = 0
        for dom in self.co_occurrence_graph.keys():
            for other_domain in self.co_occurrence_graph[str(dom)].keys():
                count += 1
        return count/(len(self.co_occurrence_graph.keys())*(len(self.co_occurrence_graph.keys()))-1)


    def top_domains_max_neighbors(self):
        """Returns the ten domains with the highest number of neighbours.
        """
        dict_neighbors = {}
        dict_domain = {}
        with open(cooccurenceDict_file_withoutNb, 'rb') as file:
            coocc_dict = pickle.load(file)
            for dom in coocc_dict.keys():
                if dom not in dict_neighbors.keys():
                    dict_neighbors[dom] = []
            dict_neighbors[dom] = len(coocc_dict[dom])

            domains = Counter(coocc_dict)
            top = domains.most_common(10)  # Finds the 10 domains which have the more neighbors
            top = dict(top)  # Converts top in a dictionary to retrieve the key values
            dict_domain = top
            return dict_domain.keys()
            # ['zf-CHY', 'IMS', 'IMS_C', 'zf-Di19', 'Mcm10', 'ClpS', 'PRT6_C', 'NosD', 'BRAP2', 'zf-UBP_var']

    def low_domains_min_neighbors(self):
        """Returns the ten domains with the lowest number of neighbours.
        """
        dict_neighbors = {}
        dict_domain = {}
        with open(cooccurenceDict_file_withoutNb, 'rb') as file:
            coocc_dict = pickle.load(file)
            for dom in coocc_dict.keys():
                if dom not in dict_neighbors.keys():
                    dict_neighbors[dom] = []
            dict_neighbors[dom] = len(coocc_dict[dom])
            domains = Counter(coocc_dict)
            n = 10
            less = domains.most_common()[:-n - 1:-1]  # Finding the ten domains with the least number of neighbors
            less = dict(less)  # Converts top in a dictionary to retrieve the key values
            dict_domain = less
            return dict_domain.keys()
        # ['AhpC-TSA', '2-Hacid_dh_C', '2-Hacid_dh', 'P4Ha_N', 'Myotub-related', 'Ribosomal_S2', 'SNF', 'His_Phos_1', 'NAD_binding_2', 'Rhodopsin_N'])


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


    def compare_most_neighbors_most_frequent(self):
        """Compares two dictionaries of interaction to check if the domains having the most neighbors are the most frequent in proteins.
        """
        most_neighbor = self.top_domains_max_neighbors()

        # Computes the domains' frequencies
        freq_dom = {}
        for dom in self.domains:
            if dom not in freq_dom:
                freq_dom[dom] = 0
        for prot in self.dict_graph.keys():
            for value in self.dict_graph[prot]["domains"]:
                for sub_value in value:
                    freq_dom[value] += 1

        # Extracts the most frequent domains
        freq_dom_10 = {}
        list_n = []
        k = Counter(freq_dom)
        high = k.most_common(10)  # Finds the 10 highest values
        for i in high:
            list_n.append(i[0])  # Adds only the domain to the list (i[1] = number of neighbours)

        # Compare the most frequent domains with those with the most neighbours
        print(most_neighbor)
        print(list_n)
        if len(most_neighbor) == len(list_n):
            for i in set(most_neighbor):
                for j in set(list_n):
                    if i == j:
                        return ('TRUE')
                    else:
                        return ('FALSE : the domains with the most neighbors are not the most frequent')


    def co_occurrence(self, dom_x, dom_y):
        """Computes and returns the co-occurrences of x and y domains number in the proteins of the interactome.

        :param dom_x: inputed domain
        :param dom_y: inputed domain
        :return: count_occ_int: number of co-occurrences where the domains x and y appear
        """
        count_occ_int = int(0)
        for i in self.dict_graph.values():
            if dom_x in i.get('domains') and dom_y in i.get('domains'):
                count_occ_int += 1
        return count_occ_int


    def self_cooccurrence(self):
        """Computes how many domains are co-occurrent using the function co_occurrence(dom_x, dom_y).
        """
        for dom in self.domains:
            number_int = self.co_occurrence(dom, dom)
        return number_int
        # Only 1 self-occurrence


    def co_occurrence_graph_to_look(self, number):
        dict_return = {}
        count_co_occurrence = 0
        count_co_occurrence_domain = 0
        for domain in self.co_occurrence_graph.keys():
            for co_occurrence_domain in self.co_occurrence_graph[str(domain)].keys():
                co_occurrence_number = self.co_occurrence_graph[str(domain)][str(co_occurrence_domain)]

                if co_occurrence_number >= number:
                    dict_return[domain] = {co_occurrence_domain: co_occurrence_number}
                    count_co_occurrence += dict_return.get(domain, {}).get(co_occurrence_domain)
                    count_co_occurrence_domain += 1
        print(self.co_occurrence_graph)
        print("Dict created with interaction with a minimum of %s contains %s domains with a total of %s interactions" %
              (number, count_co_occurrence_domain, count_co_occurrence))
        return dict_return


    def matplolib_visualization(self,number):
        #Collection creation
        domain_list=list(self.co_occurrence_graph.keys())

        G = nx.Graph()

        for domain in self.co_occurrence_graph.keys():
            for other_domain in self.co_occurrence_graph[str(domain)].keys():
                co_occurrence_number = self.co_occurrence_graph[str(domain)][str(other_domain)]
                if co_occurrence_number >= number :
                    G.add_node(domain)
                    G.add_node(other_domain)
                    G.add_edge(domain,other_domain, weight=co_occurrence_number)

        # Create positions of all nodes and save them
        pos = nx.spring_layout(G)
        # Draw the graph according to node positions
        nx.draw(G, pos, with_labels=True)
        # Create edge labels
        labels = nx.get_edge_attributes(G,'weight')
        # Draw edge labels according to node positions
        nx.draw_networkx_edge_labels(G, pos,edge_labels=labels)

        plt.show()
        return(G)


### SCRIPT TEST ###

#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).load_uniprotDict(uniprotDict_file))
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).load_co_occurenceDict(cooccurenceDict_file))
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).generate_cooccurence_graph())
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).density())  # 0.00040666944505083784
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).top_domains_max_neighbors())
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).low_domains_min_neighbors())
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).ls_domains())
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).compare_most_neighbors_most_frequent())
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).self_cooccurrence())  # 1
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).co_occurence_graph_to_look(3))
#print(DomainGraph(uniprotDict_file, cooccurenceDict_file, cooccurenceDict_file_withoutNb).matplolib_visualization(40))
