#pragma once

#ifndef VECS2PAULI_H_
#define VECS2PAULI_H_

#include "../extern/dd_package/include/dd/LimTable.hpp"
#include "../extern/dd_package/include/dd/PauliAlgebra.hpp"
#include "../extern/dd_package/include/dd/PauliUtilities.hpp"
#include <utility>
#include <complex>



namespace dd {

	bool powerOfTwo(int number);

    typedef std::pair<LimEntry<>, StabilizerGroupValue> Coset;

    struct Alpha
    {
        std::complex<double> value;
        bool allValues;
        bool noValues;

        Alpha(){
            std::complex<double> complexZero(0.0, 0.0);
            allValues = false;
            noValues = false;
            value = complexZero;
        }
    };


    struct PauliLIMCoset
    {
        Alpha alpha;
        LimEntry<> LIM;
        StabilizerGroupValue stab;
    };
    

    std::string combinePauliStrings(std::string pauliop, std::string paulistring);

    int findNumQubits(std::vector<std::complex<double>> vec);

    void printAlpha(Alpha alpha);

    bool alphaValueEqual(Alpha alpha1, Alpha alpha2);

    void resetAlpha(Alpha alpha);

    bool vecIsZero(std::vector<std::complex<double>> vec);

    void printVec(std::vector<std::complex<double>> vec);

    bool complexApproxZero(std::complex<double> complex1);

    bool complexApproxEqual(std::complex<double> complex1, std::complex<double> complex2);

    bool complexApproxMinusEqual(std::complex<double> complex1, std::complex<double> complex2);

    bool complexApproxiEqual(std::complex<double> complex1, std::complex<double> complex2);

    bool complexApproxMinusiEqual(std::complex<double> complex1, std::complex<double> complex2);

    std::complex<double> multiplyByMinusOne(std::complex<double> orig);

    std::complex<double> multiplyByi(std::complex<double> orig);

    std::complex<double> multiplyByMinusi(std::complex<double> orig);

    Coset stab2Coset(StabilizerGroupValue stab);

    StabilizerGroup toStabilizerGroup(StabilizerGroupValue groupvalue);

    StabilizerGroupValue findGeneratingSet(StabilizerGroupValue stab);

    std::pair<LimEntry<NUM_QUBITS>, bool> findCosetIntersectionElement(Coset coset1, Coset coset2, Qubit numQubits);

    Coset findCosetIntersection(Coset coset1, Coset coset2, Qubit numQubits);

    PauliLIMCoset findSolution(Coset Map1, Coset Map2, std::string pauli, StabilizerGroupValue Stab);

    StabilizerGroupValue findStabilizerGroup(std::vector<std::complex<double>> vec);

    std::pair<LimEntry<>, Alpha> handleAlpha(Alpha alpha, std::string pauliOp);

    std::pair<LimEntry<>, Alpha> vecs2PauliBase(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2, bool &foundsol);

    PauliLIMCoset vecs2Pauli(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2, bool &foundSol);
}

#endif // VECS2PAULI_H_
