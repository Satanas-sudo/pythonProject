import unittest
import main
from pathlib import Path


class TestFunction(unittest.TestCase):

    def test_read_interaction_file(self):
        """This function applies several tests to check several features. """
        interaction_file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/"
                                "Human_HighQuality.txt")

        """Tests if the output is a dictionary"""
        dict_type = type(main.read_interaction_file_dict(interaction_file))
        self.assertEqual(type(dict_type), dict)  # assertEqual() tests the equality between two parameters

        """Tests if the output is a list"""
        list_type = type(main.read_interaction_file_list(interaction_file))
        self.assertEqual(type(list_type), list)

        """Tests if the output is a list of tuples"""
        tuple_type = type(main.read_interaction_file(interaction_file))
        self.assertEqual(type(tuple_type), tuple)

    unittest.main()


"""This function checks if the interaction file is in proper format to be read correctly.

Returns True if correct format.
1st case : le fichier ne comporte pas la première ligne qui compte le nombre d'interactions
2nd case : empty file
3rd case : fichier dont la 1ere ligne contient un nb qui n'est pas le nb d'interactions
4th case : fichier contenant une ligne qui ne comporte pas le bon nombre de colonne
5th case : fichier ne contenant que des interlignes"""
def is_interaction_file(File_Empty):

    with open(File_Empty, "r") as interaction_file:
        first_char = interaction_file.read(1)  # get the first character
        if not first_char:
            print("file is empty")  # first character is the empty string
        else:
            interaction_file.seek(0)  # first character wasn't empty, returns to the beginning of the file

is_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/File_Empty.txt")

def is_interaction_file(File_Empty):
    return os.stat(File_Empty).st_size==0
    #return os.stat(File_Empty).getsize > 0
    print("File is empty")
is_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/File_Empty.txt")

"""def is_interaction_file2(Whitespaces):

    with open(Whitespaces, "r") as interaction_file:
        content = interaction_file.read()
        if interaction_file.search(r'^\s*$', content):
            #^ vérifie position au début de la ligne
            #\s match avec un espace, revient à [\ r\ n\ t\ f\ v]
            #* Quantifier — Matches between zero and unlimited times, as many times as possible, giving back as needed (greedy)
            #$ vérifie la position à la fin de la ligne
            return True
    print("file contains only blank lines")
is_interaction_file2("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Whitespaces.txt")"""

def is_interaction_file(FirstLine_Wrong):

    with open(FirstLine_Wrong, "r") as interaction_file:
        first_line = int(interaction_file.readline(1))  # Gets the first line
        lines_num = count.interaction_file()
        if first_line == (lines_num-1):
            print("number of interactions is correct")
        else:
            print("number of interactions is false")

is_interaction_file("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/FirstLine_Wrong.txt")
