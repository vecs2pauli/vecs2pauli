import vecs2pauli as vtp

# the logical basis states |0 > _L and |1 > _L of the
# 3-qubit repetition code
logical_zero = [1 , 0 , 0 , 0 , 0 , 0 , 0 , 0] # |000 >
logical_one = [0 , 0 , 0 , 0 , 0 , 0 , 0 , 1] # |111 >

# finding the stabilizers of the code
print("Stabilizers of the logical |0> state of the 3-qubit repetition code:")
stabilizer_group_zero = vtp.get_stabilizers(logical_zero) #{ZII,IZI,IIZ}
print(stabilizer_group_zero)

print("\nStabilizers of the logical |1> state of the 3-qubit repetition code:")
stabilizer_group_one = vtp.get_stabilizers(logical_one) #{-ZII,-IZI,-IIZ}
print(stabilizer_group_one)


three_qubit_code_in_stabilizer_representation = vtp.intersect_stabilizer_groups(stabilizer_group_zero,stabilizer_group_one)

print("\nStabilizers of the 3-qubit repetition code:")
print(three_qubit_code_in_stabilizer_representation)
# prints ["ZZI", "IZZ"]
