from vecs2pauli.converter import *
from vecs2pauli._vecs2pauli import _extend_check_matrix


def produce_comp_basis_vec(num_qubits, comp_basis_ix):
    basis_vec = [0 for __ in range(2 ** num_qubits)]
    basis_vec[comp_basis_ix] = 1
    return basis_vec


def full_check_matrix_to_vector(check_matrix):
    num_rows = check_matrix.shape[0]
    num_columns = check_matrix.shape[1]
    num_variables = num_columns - 1
    num_qubits = int(num_variables / 2)

    # produce projector
    common_projector = GeneratorList(data=check_matrix).to_projector()
    
    # produce vector
    comp_basis_ix = 0
    projected_vec = None
    while((projected_vec is None) or (not np.any(projected_vec))):
        basis_vec = produce_comp_basis_vec(num_qubits, comp_basis_ix)
        projected_vec = np.dot(common_projector, basis_vec)
        comp_basis_ix += 1

    return projected_vec


def _find_basis_for_stabilised_subspace(check_matrix, output_has_unit_norm=True):
    """
    Example usage:
    >>> check_matrix = np.array([[True, True, False, False, False]])
    >>> vecs = find_basis_for_stabiliser_subspace(check_matrix=check_matrix, output_has_unit_norm=False)
    >>> print(vecs)
    >>> # [array([2.+0.j, 0.+0.j, 0.+0.j, 2.+0.j]), array([0.+0.j, 2.+0.j, 2.+0.j, 0.+0.j])]
    """
    num_rows = check_matrix.shape[0]
    num_columns = check_matrix.shape[1]
    num_variables = num_columns - 1
    num_qubits = int(num_variables / 2)

    extended_check_matrix = _extend_check_matrix(check_matrix)

    ret = []
    for vec_ix in range(2 ** (num_qubits - num_rows)):
        extended_cm_with_phases_changed = np.copy(extended_check_matrix) 
        ix_in_binary = ("{0:0" + str(num_qubits - num_rows) + "b}").format(vec_ix)
       
        # change phases
        for j in range(len(ix_in_binary)):
            if ix_in_binary[j] == "1":
                extended_cm_with_phases_changed[j + num_rows][num_columns - 1] = True
        # get vector
        vec = full_check_matrix_to_vector(check_matrix=extended_cm_with_phases_changed)
        if(output_has_unit_norm):
            vec = vec / np.linalg.norm(vec)
        ret.append(vec)
    return ret


def find_basis_for_stabilised_subspace(stabilizer_generators, output_has_unit_norm=True):
    # check if input is a check matrix
    generator_list = GeneratorList(data=stabilizer_generators)
    check_matrix = generator_list.to_check_matrix()
    return _find_basis_for_stabilised_subspace(check_matrix=check_matrix, output_has_unit_norm=output_has_unit_norm)
