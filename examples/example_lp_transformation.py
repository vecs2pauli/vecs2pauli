import vecs2pauli as vtp


source = ([1, 2, 3, 4, 5, 6, 7, 8])
target = ([7j, -8j, -5j, 6j, 3j, -4j, -1j, 2j])

transformations = vtp.get_local_pauli_transformations(source, target)
print(transformations)
