from pathlib import Path
import csv
import urllib.request
import tempfile
import re


### PRELIMINARY QUESTIONS ###

# What kind of informations contains Pfam (take INSR_HUMAN for example) ?
# It countains: Entry, Protein names, Gene names, Organism, Length, Entry name
# As there is a great heterogeneity of gene or protein names (same official nomenclature, which contains synonyms, and of
# which the syntaxis is not always followed correctly), it is more cautious to use a identifier. What is the Uniprot
# identifier of INSR_HUMAN 2 ?
# P06213

# In Pfam, what are the domains of INSR_HUMAN ? Nevertheless, the information in the Pfam datanase is not very handy to
# get by a program. Do you retrieve this information in the Uniprot page ?
# I found 4 domains in the Uniprot page: Fibronectin type-III 1, Fibronectin type-III 2, Fibronectin type-III 3 and
# Protein kinase.

########### PROTEINS IDENTIFIERS RETRIEVAL ###########

file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")
uniprotFile = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality_uniprot.tab")  # Only proteins from Human_HighQuality.txt in Uniprot
proteome_uniprotFile = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_proteome_uniprot.tab")  # Entire human proteome in Uniprot

class Proteome:

    # Attributes
    prot_list = []  # List of proteins of the proteome
    name_ID_dict = {}  # This dictionary contains the proteins' names as keys and the Uniprot ID identifier as value
    ID_name_dict = {}  # This dictionary countains the Uniprot ID identifier as key and the associated protein's name as value

    # Constructor
    def __init__(self, uniprotFile):
        self.uniprotFile = str(uniprotFile)
        self.prot_list = self.get_protein_name_list(uniprotFile)
        self.name_ID_dict = self.get_name_uniprotID_dict(uniprotFile)
        self.ID_name_dict = self.get_uniprotID_name_dict(uniprotFile)

    # Methods

########### PROTEINS DOMAINS RETRIEVAL ###########

### PRELIMINARY QUESTIONS ###
# To retrieve information from the Uniprot page of a protein, it is easier to go through textual version of Uniprot than
# HTML pages.
# 1. From HTML page of "P06213 3", find Pfam domains.
# ['Furin', 'Insulin_TMD', 'PK_Tyr_Ser', 'Recep_L_domain']
# 2. By adding .txt to the url a the protein's page, we get a version manageable more easily.
# Where are the Pfam domains that we identified previously ? How do we know how many times each domain is present in the protein ?

    def get_protein_name_list(self, uniprotFile):
        """Reads the Uniprot file and return its elements in a list.

        :param uniprotFile: interaction file
        :return: prot_list: list of proteins' names
        """
        uniprotFile = csv.reader(open(uniprotFile), delimiter="\t")
        prot_list = []
        for line in uniprotFile:
            if line:
                prot_list.append(line[0])
        return prot_list

    def get_name_uniprotID_dict(self, uniprotFile):
        """Reads the Uniprot file and creates a dictionary with the protein's name associated with the Uniprot ID.

        :param uniprotFile: interaction file
        :return: nameID_dict: dictionary of proteins' names with their corresponding Uniprot ID.
        """
        uniprotFile = csv.reader(open(uniprotFile), delimiter="\t")
        nameID_dict = {}
        for line in uniprotFile:
            if line:
                nameID_dict.update({line[1]:line[0]})
        return nameID_dict

    def get_uniprotID_name_dict(self, uniprotFile):
        """Reads the Uniprot file and creates a dictionary with the UniprotID associated with the protein's name.

        :param uniprotFile: interaction file
        :return: nameID_dict: dictionary of the proteins' UniprotID associated with their names.
        """
        uniprotFile = csv.reader(open(uniprotFile), delimiter="\t")
        IDname_dict = {}
        for line in uniprotFile:
            if line:
                IDname_dict.update({line[0]:line[1]})
        return IDname_dict

    def get_protein_domains_list(self, prot):
        """This function allows to search for the Uniprot page corresponding to a protein entry and returns the domain(s)'s
        name(s) contained by that protein. urllib, tempfile (creates temporary files) and re are needed."""

        url_str = "https://www.uniprot.org/uniprot/{}.txt".format(prot)

        with urllib.request.urlopen(url_str) as rep_request:
            with tempfile.TemporaryFile() as tmp_file:  # Creation of the temporary file which contains the request result
                tmp_file.write(rep_request.read())
                tmp_file.seek(0)  # All document are found (argument 0)
                # We search for the lines that have this regular expression
                # Use of re to find the lines beginning by "DR   Pfam" and stocks them into a list
                result_domains = re.findall(r'DR\s{3}Pfam;\s\w+;\s\w+', str(tmp_file.read()))
                # 4th elements of each lines (domain name) are extracted and returned as a list
                domain_names_list = []

                for element_str in result_domains:
                    domain_names_list.append(element_str.split()[3])  # The domain is the 3rd element of each line found
                return domain_names_list

### SCRIPT TEST ###

#print(Proteome(uniprotFile).get_protein_name_list(uniprotFile))
#print(Proteome(uniprotFile).get_name_uniprotID_dict(uniprotFile))
#print(Proteome(uniprotFile).get_uniprotID_name_dict(uniprotFile))
#print(Proteome(uniprotFile).get_protein_domains_list("P06213"))