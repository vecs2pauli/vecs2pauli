
import vecs2pauli as vtp

# the usual logical basis states |0>_L and |1>_L of the 5-qubit Shor code
logical_zero = [1, 0, 0, -1, 0, 1, -1, 0, 0, 1, 1, 0, -1, 0, 0, -1, 0, -1, 1, 0, 1, 0, 0, -1, -1, 0, 0, -1, 0, -1, -1, 0]
logical_one = [0, -1, -1, 0, -1, 0, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0, -1, 0, 0, -1, 0, 1, 1, 0, 0, -1, 1, 0, -1, 0, 0, 1]
minus_logical_one = [-1 * x for x in logical_one]

# finding a logical X operator
print(vtp.get_local_pauli_transformations(logical_zero, logical_one))
# prints ((1-0j), 'iIZZIY', ['ZYYZI', 'XIXZZ', '-IZIXX', 'IXZZX', '-IIYZY'])
# i.e. all possible logical X operators are of the form iIZZIY * g
# where g is either IIIII or can be written as product of a subset of
# ['ZYYZI', 'XIXZZ', '-IZIXX', 'IXZZX', '-IIYZY']


# finding a logical Z operator
transf_satisfying_zero = vtp.get_local_pauli_transformations(logical_zero, logical_zero)
transf_satisfying_one = vtp.get_local_pauli_transformations(logical_one, minus_logical_one)
print(vtp.intersect_cosets(
    transf_satisfying_zero[0], transf_satisfying_zero[1], transf_satisfying_zero[2],
    transf_satisfying_one[0], transf_satisfying_one[1], transf_satisfying_one[2])
      )
# prints ((1+0j), '-IZIXX', ['+ZYYZI', '+XIXZZ', '+IXZZX', '+IZYYZ'])
# TODO can we learn anything from the fact that the stabiliser group is not maximal?
# TODO I think we want Coset.includes(lp_tarnsformation) and StabGroup.includes(lp_transformation) functions, that return a boolean, to check if the resulting coset contains the usual logical Z operator ZZZZZ
