from vecs2pauli._vecs2pauli import _get_local_pauli_transformations, _get_stabilizers, _intersect_stabilizer_groups, _intersect_cosets
from vecs2pauli.converter import _add_plus_in_front_of_string_if_necessary, _remove_plus_in_front_of_string_if_present




def get_local_pauli_transformations(source_vec, target_vec):
    """
    Output convention:
    - 'None' if no transformations possible
    - string "all" if all transformations possible
    - otherwise: (factor, Pauli string, group generators)
    """
    # TODO parse output as done in 'get_stabilizers'
    (factor, pauli, group, all_possible, no_possible) = _get_local_pauli_transformations(source_vec, target_vec)
    if all_possible:
        return "all"
    elif no_possible:
        return None
    else:
        return (factor, pauli, group)

def get_stabilizers(vector):
    generator_strings = _get_stabilizers(vector)
    return [_add_plus_in_front_of_string_if_necessary(string) for string in generator_strings]

def intersect_stabilizer_groups(generators_a, generators_b):
    gens_a = [_remove_plus_in_front_of_string_if_present(string) for string in generators_a]
    gens_b = [_remove_plus_in_front_of_string_if_present(string) for string in generators_b]
    gens = _intersect_stabilizer_groups(gens_a, gens_b)
    return [_add_plus_in_front_of_string_if_necessary(string) for string in gens]


def intersect_cosets(factor_a, coset_repr_a, generators_a, factor_b, coset_repr_b, generators_b):
    """
    Parameters
    ----------
    factor_a : complex number, None or "any"
    factor_b : complex number, None or "any"
    """
    if factor_a == "any":
        return (factor_b, coset_repr_b, generators_b)
    if factor_b == "any":
        return (factor_a, coset_repr_a, generators_a)
    if factor_a is None or factor_b is None:
        return (None, coset_repr_a, generators_a)
    if factor_a == 0.0:
        if factor_b == 0.0:
            return (factor_a, coset_repr_a, generators_a)
        else:
            return (None, coset_repr_a, generators_a)
    if factor_b == 0.0:
        return (None, coset_repr_a, generators_a)
    if factor_a != factor_b:
        return (None, coset_repr_a, generators_a)
    else:
        gens_a = [_remove_plus_in_front_of_string_if_present(string) for string in generators_a]
        gens_b = [_remove_plus_in_front_of_string_if_present(string) for string in generators_b]
        coset_repr, gens = _intersect_cosets(factor_a, coset_repr_a, gens_a, factor_b, coset_repr_b, gens_b)
        return (factor_a, coset_repr, [_add_plus_in_front_of_string_if_necessary(string) for string in gens])
