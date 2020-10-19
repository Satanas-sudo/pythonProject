import unittest
import graphReading
from pathlib import Path
import time

interaction_file = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Human_HighQuality.txt")
file_empty = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/File_Empty.txt")
first_line_wrong = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/FirstLine_Wrong.txt")
wrong_num_column = Path("C:/Users/a/Documents/Cours_rennes1_master_bioinfo/S9_ADG/projet_python/Line_WrongColumn.txt")


class TimeFunction:

    def time_read_interaction_file(self):
        """This function measures the time in seconds to execute the functions read_interaction_file_dict and
        read_interaction_file_list in order to compare their performances."""

        dict_start = time.time()
        graphReading.read_interaction_file_dict(interaction_file)
        dict_finished = time.time() - dict_start

        print("Time to execute the function read_interaction_file_dict is {} seconds".format(dict_finished))


        list_start = time.time()
        graphReading.read_interaction_file_list(interaction_file)
        list_finished = time.time() - list_start

        print("Time to execute the function read_interaction_file_list is {} seconds".format(list_finished))

# Making a list is faster than making a dictionary


class TestFunction(unittest.TestCase):

    def test_read_interaction_file(self):
        """This function applies several tests to check several features. """

        """Tests if the output is a dictionary"""
        dict_test = graphReading.read_interaction_file_dict(interaction_file)
        dict_type = type(dict_test)
        self.assertEqual(dict_type, dict)  # assertEqual() tests the equality between two parameters

        """Tests if the output is a list"""
        list_test = graphReading.read_interaction_file_list(interaction_file)
        list_type = type(list_test)
        self.assertEqual(list_type, list)

        """Tests if the output is a list of tuples"""
        tuple_test = graphReading.read_interaction_file(interaction_file)
        tuple_type = type(tuple_test)
        self.assertEqual(tuple_type, tuple)

        """Tests if a file is empty"""
        first_char = int(file_empty.read(1))  # Gets the first character (if there is one, then it equals to 1)
        #print(type(first_char))
        if not first_char:
            no_char = 0  # If there is not a character in the file, then it equals to 0
        self.assertEqual(first_char, no_char)

        """Tests if a file whose first line has a wrong number which counts the number of interactions"""
        first_line = int(first_line_wrong.readline(1))  # Gets the first line
        lines_num = len(first_line_wrong.readlines())  # Counts the number of lines
        interaction_int = lines_num - 1
        self.assertEqual(first_line, interaction_int)

        """Tests if a file has not a first line with the number of interactions"""
        # TO DO

        """Tests if the number of columns is correct"""
        # TO DO

    unittest.main()
