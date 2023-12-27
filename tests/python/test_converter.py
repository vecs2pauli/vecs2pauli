import unittest
from vecs2pauli.converter import Generator, GeneratorList, PauliX, PauliY, PauliI, PauliZ
import numpy as np


class TestGenerator(unittest.TestCase):

    def test(self):
        string = "+X"
        check_row = np.array([True, False, False])
        for inp in [string, check_row]:
            gen = Generator(data=inp)
            assert(gen._phase == 1)
            assert(len(gen._paulis) == 1)
            assert(isinstance(gen._paulis[0], PauliX))
            assert(str(gen) == "+X")
    
        string = "-I"
        check_row = np.array([False, False, True])
        for inp in [string, check_row]:
            gen = Generator(data=inp)
            assert(gen._phase == -1)
            assert(len(gen._paulis) == 1)
            assert(isinstance(gen._paulis[0], PauliI))
            assert(str(gen) == "-I")
    
        string = "-YYY"
        check_row = np.array([True, True, True, True, True, True, True])
        for inp in [string, check_row]:
            gen = Generator(data=inp)
            assert(gen._phase == -1)
            assert(len(gen._paulis) == 3)
            assert(isinstance(gen._paulis[0], PauliY))
            assert(isinstance(gen._paulis[1], PauliY))
            assert(isinstance(gen._paulis[2], PauliY))
            assert(str(gen) == "-YYY")

class TestGeneratorList(unittest.TestCase):

    def test(self):

        check_matrix = np.array([[True, True, False, False, False], [False, False, True, True, False]])
        string_list = ["+XX", "+ZZ"]
        for inp in [string_list, check_matrix]:
            genlist = GeneratorList(data=inp)
            assert(genlist.to_check_matrix().shape == (2, 5))
            outputted_cm = genlist.to_check_matrix()
            for row_ix in range(2):
                for column_ix in range(5):
                    assert(outputted_cm[row_ix, column_ix] == check_matrix[row_ix, column_ix])
            assert(genlist.to_string_list() == string_list)
    
    
        check_matrix = np.array([[True, True, False, False, False], [True, True, True, True, True]])
        string_list = ["+XX", "-YY"]
        for inp in [string_list, check_matrix]:
            genlist = GeneratorList(data=inp)
            assert(genlist.to_check_matrix().shape == (2, 5))
            outputted_cm = genlist.to_check_matrix()
            for row_ix in range(2):
                for column_ix in range(5):
                    assert(outputted_cm[row_ix, column_ix] == check_matrix[row_ix, column_ix])
            assert(genlist.to_string_list() == string_list)
    
    
        check_matrix = np.array([[False, False, False, False, False], [True, False, True, False, True]])
        string_list = ["+II", "-YI"]
        for inp in [string_list, check_matrix]:
            genlist = GeneratorList(data=inp)
            assert(genlist.to_check_matrix().shape == (2, 5))
            outputted_cm = genlist.to_check_matrix()
            for row_ix in range(2):
                for column_ix in range(5):
                    assert(outputted_cm[row_ix, column_ix] == check_matrix[row_ix, column_ix])
            assert(genlist.to_string_list() == string_list)


if __name__ == "__main__":
    unittest.main()
