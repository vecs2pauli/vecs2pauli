import vecs2pauli as vtp
import numpy as np

# generators of the stabilizer code
stabilizer_generators = ["+ZZI", "+IZZ"]

unnormalized_basis_vectors = vtp.find_basis_for_stabilised_subspace(stabilizer_generators)

normalized_basis_vectors = [vector / np.linalg.norm(vector) for vector in unnormalized_basis_vectors]

print("An orthogonal basis of stabilizer states which span the code space of the code", stabilizer_generators, ":")
print(normalized_basis_vectors)
# prints a list of [1, 0, 0, 0, 0, 0, 0, 0] and [0, 0, 0, 0, 0, 0, 0, 1],
# i.e. |000> and |111>
