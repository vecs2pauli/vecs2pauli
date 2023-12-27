import vecs2pauli as vtp
import numpy as np


num_qubits = 3
state = [3, -1, -1, -1, -1, -1, -1, 3]

# finding the generators of all stabilizers
stabilizer_generators = vtp.get_stabilizers(state)
print("\nGenerators of the stabilizer group of the state", state, ":")
print(stabilizer_generators)
# prints ["+XXX"]

print("\nAn upper bound to the stabilizer rank of the state is therefore:", 2 ** (num_qubits - len(stabilizer_generators)))

# find a basis of stabilizer states in which 'state' can be decomposed
basis = vtp.find_basis_for_stabilised_subspace(stabilizer_generators)
normalized_basis = [vector / np.linalg.norm(vector) for vector in basis]

print("\nA basis of stabilizer states in which the state", state, "can be decomposed:\n")
print(normalized_basis)
print("...with coefficients:")
print([np.vdot(state, vector) for vector in normalized_basis])
