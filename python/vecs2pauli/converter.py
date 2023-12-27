import numpy as np
from vecs2pauli._vecs2pauli import _expand_generators

def _add_plus_in_front_of_string_if_necessary(string):
    if string[0] != "-":
        return "+" + string
    else:
        return string

def _remove_plus_in_front_of_string_if_present(string):
    if string[0] == "+":
        return string[1:]
    else:
        return string




class Pauli:

    string = None
    check_pair = None
    matrix = None

    def __str__(self):
        return self.string

    def commutes_with(self, other):
        if not isinstance(other, Pauli):
            raise TypeError


class PauliI(Pauli):

    string = "I"
    check_pair = [False, False]
    matrix = np.array([[1, 0], [0, 1]], dtype=complex)

    def commutes_with(self, other):
        super().commutes_with(other=other)
        return True


class PauliX(Pauli):

    string = "X"
    check_pair = [True, False]
    matrix = np.array([[0, 1], [1, 0]], dtype=complex)

    def commutes_with(self, other):
        super().commutes_with(other=other)
        return isinstance(other, PauliI) or isinstance(other, PauliX)


class PauliY(Pauli):

    string = "Y"
    check_pair = [True, True]
    matrix = np.array([[0, -1j], [1j, 0]], dtype=complex) 

    def commutes_with(self, other):
        super().commutes_with(other=other)
        return isinstance(other, PauliI) or isinstance(other, PauliY)


class PauliZ(Pauli):

    string = "Z"
    check_pair = [False, True]
    matrix = np.array([[1, 0], [0, -1]], dtype=complex)

    def commutes_with(self, other):
        super().commutes_with(other=other)
        return isinstance(other, PauliI) or isinstance(other, PauliZ)


class Generator:
    # TODO rename to LocalPauliTransformation? and then assert that the factor is +-1?
    """
    Represents a generators of a stabiliser group.

    Parameters
    ----------
    data : :obj:`~numpy.array` of Booleans or str

    Example
    -------
    >>> string = "-X"
    >>> gen = Generator(data=string)
    >>> print(gen)
    >>> # prints "-X"
    >>>
    >>> check_row = np.array([True, False, True])  # -X
    >>> gen = Generator(data=check_row)
    >>> # prints "-X"
    """

    def __init__(self, data):
        self._phase = None
        self._paulis = None
        self._num_qubits = None
        try:
            self.from_string(string=data)
        except TypeError:
            try:
                self.from_check_row(check_row=data)
            except TypeError:
                raise TypeError

    @property
    def num_qubits(self):
        return self._num_qubits

    def from_check_row(self, check_row):
        if not type(check_row).__module__ == np.__name__:
            raise TypeError("check_row is not a numpy array")
        length = check_row.shape[0]
        if not length % 2 == 1:
            raise TypeError("number of columns is check_row is not twice an integer plus one")
        if length < 3:
            raise TypeError("number of columns of check_row should be at least 3")
        self._num_qubits = int((length - 1) / 2)


        if check_row.dtype != bool:
            if check_row.dtype == int or check_row.dtype == float:
                # convert to bools
                num_columns = check_row.shape[0]

                check_row_boolean = np.array([False for __ in range(num_columns)])
                for column_ix in range(num_columns):
                    if int(check_row[column_ix]) == 1:
                        check_row_boolean[column_ix] = True
                    elif int(check_row[column_ix]) == 0:
                        check_row_boolean[column_ix] = False 
                    else:
                        raise TypeError("check_row is not an array of booleans or 0/1 ints")
                check_row = check_row_boolean
            else:
                raise TypeError("check_row is not an array of booleans or 0/1 ints")


        if check_row[-1]:
            self._phase = -1
        else:
            self._phase = 1
        self._paulis = []

        # set individual paulis
        for qubit_ix in range(self.num_qubits):
            if check_row[qubit_ix]:
                if check_row[qubit_ix + self.num_qubits]:
                    self._paulis.append(PauliY())
                else:
                    self._paulis.append(PauliX())
            else:
                if check_row[qubit_ix + self.num_qubits]:
                    self._paulis.append(PauliZ())
                else:
                    self._paulis.append(PauliI())

    @property
    def check_row(self):
        return Generator(data=self.to_check_row()).to_check_row()

    def to_check_row(self):
        return np.array([pauli.check_pair[0] for pauli in self._paulis] +\
                        [pauli.check_pair[1] for pauli in self._paulis] +\
                        [self._phase == -1])

    def __str__(self):
        return self.string()

    def from_string(self, string):
        if not isinstance(string, str):
            raise TypeError("string is not of string type")
        if len(string) < 2:
            raise TypeError("stabilizer as a string should have at least two characters")
        self._paulis = []

        phase = string[0]
        if phase == "+":
            self._phase = 1
        elif phase == "-":
            self._phase = -1
        else:
            raise TypeError

        self._num_qubits = len(string) - 1
        for ix in range(self._num_qubits):
            for pauli_cls in [PauliI, PauliX, PauliY, PauliZ]:
                if pauli_cls.string == string[ix + 1]:
                    self._paulis.append(pauli_cls())

    def string(self):
        if self._phase == 1:
            s = "+"
        else:
            s = "-"
        for pauli in self._paulis:
            s += str(pauli)
        return s

    def commutes_with(self, other):
        if not isinstance(other, Generator):
            raise TypeError
        else:
            if self.num_qubits != other.num_qubits:
                raise TypeError
            ret = True
            for ix in range(self.num_qubits):
                if not self._paulis[ix].commutes_with(other._paulis[ix]):
                    ret = not ret
            return ret

    def to_matrix(self):
        res = np.identity(1, complex)
        for qubit_ix in range(self.num_qubits):
            res = np.kron(res, self._paulis[qubit_ix].matrix)
        if(self._phase == -1):
            res *= -1
        return res

    def to_projector(self, positive_subspace=True):
        for row_ix in range(self.num_qubits):
            pauli = self.to_matrix()
            phase = 1 if positive_subspace else -1
            return np.identity(2 ** self.num_qubits, complex) + phase * pauli



class GeneratorList:
    # TODO make subscriptable, i.e. so glist[0] becomes possible. Also add this to a test
    # TODO add function reduce() which reduces to a list of only independent elements
    # TODO overwrite __repr__
    # TODO do we want a property 'group_size', which is nontrivial to compute because GeneratorList might contain '+I' or not be reduced to only contain independent elements. Also add this to a test
    """
    Represents a set of generators of a stabiliser group.

    Parameters
    ----------
    data : list of :obj:`~vecs2pauli.converter.Generator` or :obj:`~numpy.array` of Booleans or list of str
    """
    
    def __init__(self, data=None):
        self._reduced = False
        if data is None:
            self.generators = None
        else:
            self.generators = []
            try:
                self.from_check_matrix(check_matrix=data)
            except:
                try:
                    self.from_string_list(strings=data)
                except:
                    raise TypeError("Input to GeneratorList not in known format")

    @property
    def num_generators(self):
        # TODO rename to 'size'?
        # TODO what if the generators are not independent?
        return len(self.generators)

    @property
    def group_size(self):
        # TODO first call 'reduce()'
        # TODO remove the following line
        new_generators = []
        for g in self.generators:
            if str(g) != "+" + "I" * self.num_qubits:
                new_generators.append(g)
        self.generators = new_generators
        return 2 ** self.num_generators

    def reduce(self):
        if not self._reduced:
            # TODO fill
            pass
            self._reduced = True

    def add(self, generator, reduce=False):
        if reduce:
            # TODO how to set 'self._reduced'? Should make a call
            # TODO
            self._reduced = True
        else:
            self._reduced = False

    @property
    def num_qubits(self):
        if self.num_generators > 0:
            return self.generators[0].num_qubits
        else:
            return None

    def _check_same_num_qubits(self):
        if self.num_generators > 0:
            for generator in self.generators:
                assert(generator.num_qubits == self.num_qubits)

    def from_check_matrix(self, check_matrix):
        if not type(check_matrix).__module__ == np.__name__:
            raise TypeError("check_matrix is not a numpy array of booleans")

        num_rows = check_matrix.shape[0]
        self.generators = []
        for row_ix in range(num_rows):
            check_row = check_matrix[row_ix, 0:]
            generator = Generator(data=check_row)
            self.generators.append(generator)
        self._check_same_num_qubits() 
        self._check_all_commute()

    def to_check_matrix(self):
        return np.array([list(generator.to_check_row()) for generator in self.generators])

    def from_string_list(self, strings):
        if not isinstance(strings, list):
            raise TypeError("strings is not of type list")
        for ix in range(len(strings)):
            if not isinstance(strings[ix], str):
                raise TypeError("element {} in list strings is not a string", strings[ix])
        self.generators = [Generator(data=string) for string in strings]
        self._check_same_num_qubits() 
        self._check_all_commute()

    def to_string_list(self):
        return [str(generator) for generator in self.generators]

    def _check_all_commute(self):
        for gen_a in self.generators:
            for gen_b in self.generators:
                assert(gen_a.commutes_with(gen_b))

    def to_projector(self, positive_subspace=None):
        """
        Parameters
        ----------
        positive_subspace : list of bool
        """
        if positive_subspace is None:
            positive_subspace = [True for __ in range(self.num_qubits)]
        else:
            assert(len(positive_subspace) == self.num_qubits)
            for el in positive_subspace:
                assert(isinstance(el, bool))

        ret = np.identity(2 ** self.num_qubits, complex)
        for generator_ix in range(self.num_generators):
            generator = self.generators[generator_ix]
            proj = generator.to_projector(positive_subspace=positive_subspace[generator_ix])
            ret = np.matmul(ret, proj)
        return ret

    def includes(self, generator):
        assert(isinstance(generator, Generator))
        # TODO currently slow, make faster!
        # TODO remove need for removing the 'plus'
        group = _expand_generators([_remove_plus_in_front_of_string_if_present(str(g)) for g in self.generators])
        return (_remove_plus_in_front_of_string_if_present(str(generator)) in group)




def check_row2stabilizer_generator(check_row):
    """
    Converts a row of a check matrix to a strings.

    Example
    -------
    >>> check_row = np.array([True, True, False, False, False])
    >>> generator = check_row2stabilizer_generator(check_row)
    >>> print(generator)
    >>> # prints "+XX"
    """
    return str(Generator(data=check_row))

def check_matrix2stabilizer_generators(check_matrix):
    """
    Converts a check matrix to a list of strings.

    Example
    -------
    >>> check_matrix = np.array([[True, True, False, False, False], [False, False, True, True, True]])
    >>> generators = check_matrix2stabilizer_generators(check_matrix)
    >>> print(generators)
    >>> # prints ["+XX", "-ZZ"]
    """
    return str(Generator(data=check_matrix))
