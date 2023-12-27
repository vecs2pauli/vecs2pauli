import unittest
import vecs2pauli as vtp


class TestGetStabilizers(unittest.TestCase):

    def test_get_stabilizers_single_qubit_stabilizer_states(self):
        # |0>
        vec = [1.0, 0.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "+Z")

        # |1>
        vec = [0.0, 1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "-Z")

        # |+>
        vec = [1.0, 1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "+X")

        # |->
        vec = [1.0, -1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "-X")

        # |+i>
        vec = [1.0, 1.0j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "+Y")

        # |+i>
        vec = [1.0, -1.0j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "-Y")

        # 2|+i>
        vec = [2.0, -2.0j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(str(gens.generators[0]) == "-Y")


    def test_get_stabilizers_single_qubit_nonstabilizer_states(self):

        vec = [0.0, 0.0]
        print(vtp.get_stabilizers(vec))
        # TODO now what?
#        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
#        print(gens)
#        self.assertTrue(gens.group_size == 1)
#        self.assertTrue(str(gens.generators[0]) == "-Y")


        vec = [2.0, 1.0j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        #TODO we want this to equal 0, and that 'III' is always a generator
        self.assertTrue(gens.group_size == 1)


        vec = [2.0, -0.1 + 1.0j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        #TODO we want this to equal 0, and that 'III' is always a generator
        self.assertTrue(gens.group_size == 1)


        vec = [-0.5 + 0.3j, 0.7j]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        #TODO we want this to equal 0, and that 'III' is always a generator
        self.assertTrue(gens.group_size == 1)


    def test_get_stabilizers_multi_qubit_stabilizer_states(self):

        # |00> + |11>
        vec = [1.0] + [0.0] * 2 + [1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 4)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+ZZ")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+XX")))

        # |00> + |11>
        vec = [1.0] + [0.0] * 2 + [1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 4)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+ZZ")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+XX")))

        # |000> + |111>
        vec = [1.0] + [0.0] * 6 + [1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 8)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+IZZ")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+ZZI")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+XXX")))

        # |00000> + |11111>
        vec = [1.0] + [0.0] * 30 + [1.0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 32)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+IZZII")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+ZZIII")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+IIZZI")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+IIIZZ")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="+XXXXX")))


    def test_get_stabilizers_multi_qubit_nonstabilizer_states(self):

        # |001> + |010> + |100>
        vec = [0, 1, 1, 0, 1, 0, 0, 0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-ZZZ")))

        # 0.3|001> + 0.5|010> - 0.2i|100>
        vec = [0, 0.3, 0.5, 0, -0.2j, 0, 0, 0]
        print("here:")
        print([str(g) for g in gens.generators])
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-ZZZ")))

        # 0.3|001> + 0.5|010> - 0.2i|100> + 0.7|111>
        vec = [0, 0.3, 0.5, 0, -0.2j, 0, 0, 0.7]
        print("here:")
        print([str(g) for g in gens.generators])
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 2)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-ZZZ")))

        # 0.3|001> + 0.5|010> - 0.2i|100> + 0.7|110>
        vec = [0, 0.3, 0.5, 0, -0.2j, 0, 0.7, 0]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 1)


        vec = [2, 3, -2j, 3j, 2j, 3j, -2, 3]
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        self.assertTrue(gens.group_size == 4)
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-ZYZ")))
        self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-XXZ")))

        # zero vector
        vec = [0] * 16
        gens = vtp.converter.GeneratorList(data=vtp.get_stabilizers(vec))
        # TODO fix
        #self.assertTrue(gens.group_size == 4)
        #self.assertTrue(gens.includes(generator=vtp.converter.Generator(data="-IIII")))


if __name__ == "__main__":
    unittest.main()
