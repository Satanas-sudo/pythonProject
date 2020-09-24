import unittest
import main

class TestFunction(unittest.TestCase):
    def test_interaction_file_dict(self):
        dict_type=type(main.read_interaction_file_dict("Human_HighQuality.txt"))
        self.assertEqual(type(dict_type), dict)  # assertEqual() teste l'égalité entre les deux paramètres, ici que le type du dictionnaire soit bien un dictionnaire)

    def test_interaction_file_list(self):
        list_type=type(main.read_interaction_file_list("Human_HighQuality.txt"))
        self.assertEqual(type(list_type), list)

    def test_interaction_file(self):
        tuple_type=type(main.read_interaction_file("Human_HighQuality.txt"))
        self.assertEqual(type(tuple_type), tuple)

if __name__ == '__main__':
    unittest.main()
