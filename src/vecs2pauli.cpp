#include "vecs2pauli.hpp"

#include <iostream>
#include <vector>
#include <utility>
#include <complex>
#include <fstream>
#include <cassert>
#include <sstream>

namespace dd{

bool powerOfTwo(int number){
    int checker = 2;
    while (checker <= number){
        if (checker == number ){
            return true;
        }
        checker *= 2;
    }
    return false;
}


std::string makeOutputstring(dd::PauliLIMCoset res, std::vector<std::complex<double>> vec1){
    std::stringstream output;
    output << "(";
    if (res.alpha.allValues){
        output << "all values";
    }
    else if (res.alpha.noValues){
        output << "no values";
    }
    else {
        output << std::real(res.alpha.value);
        output << "+";
        output << std::imag(res.alpha.value);
        output << "i";
    }
    output << ", ";
    output << dd::LimEntry<>::to_string(&(res.LIM), dd::findNumQubits(vec1) - 1);
    output << ", {";
    if (int(res.stab.size() > 1)){
        for (int i = 0; i < int(res.stab.size()-1); i++){
            output << dd::LimEntry<>::to_string(&(res.stab[i]), dd::findNumQubits(vec1) - 1);
            output << ", ";
        }
    }
    if (int(res.stab.size() > 0)){
        output << dd::LimEntry<>::to_string(&(res.stab[res.stab.size()-1]), dd::findNumQubits(vec1) - 1);
    }
    output << "})";
    return (output.str());
}





    bool PRINT = false;

    // returns the string corresponding to the tensor product of a single pauliop and an existing paulistring
    std::string combinePauliStrings(std::string pauliop, std::string paulistring){
        std::string res ("");
        if (paulistring[0] == '-'){
            res += "-";
            if (paulistring.size() > 1 && paulistring[1] == 'i'){
                res += "i";
                res += pauliop;
                res.append(paulistring.begin()+2,paulistring.end());
            }
            else {
                res += pauliop;
                res.append(paulistring.begin()+1,paulistring.end());
            }
        }
        else if (paulistring[0] == 'i'){
            res += "i";
            res += pauliop;
            res.append(paulistring.begin()+1,paulistring.end());
        }
        else {
            res += pauliop;
            res += paulistring;
        }
        return res;
    }

    int findNumQubits(std::vector<std::complex<double>> vec){
        int length = vec.size();
        int numQubits = 0;
        while (length > 1){
            length /= 2;
            numQubits++;
        }
        return numQubits;
    }

    void printAlpha(Alpha alpha){
        if (alpha.allValues){
            std::cout << "all values";
        }
        else if (alpha.noValues){
            std::cout << "no values";
        }
        else {
            std::cout << std::real(alpha.value) << "+" << std::imag(alpha.value) << "i";
        }
    }

    bool alphaValueEqual(Alpha alpha1, Alpha alpha2){
        if ((alpha1.allValues && !alpha2.noValues) || (!alpha1.noValues && alpha2.allValues)){
            return true;
        }
        else if (!alpha1.noValues && !alpha2.noValues && complexApproxEqual(alpha1.value, alpha2.value)){
            return true;
        }
        return false;
    }

    void resetAlpha(Alpha alpha){
        alpha.allValues = false;
        alpha.noValues = false;
        std::complex<double> complexOne(1.0, 0.0);
        alpha.value = complexOne;
    }

    bool vecIsZero(std::vector<std::complex<double>> vec){
        bool foundNonZero = false;
        std::complex<double> vecZero(0.0, 0.0);
        for (int i = 0; i <= int(vec.size() - 1); i++){
            if (!complexApproxEqual(vec[i], vecZero)){
                foundNonZero = true;
            }
        }
        return (!foundNonZero);
    }

    void printVec(std::vector<std::complex<double>> vec){
        std::cout << "(";
        for (int i = 0; i < int(vec.size() - 1); i++){
            std::cout << std::real(vec[i]) << "+" << std::imag(vec[i]) << "i, ";
        }
        std::cout << std::real(vec[vec.size()-1]) << "+" << std::imag(vec[vec.size()-1]) << "i)" << std::endl;
    }

    bool complexApproxZero(std::complex<double> complex1){
        if (std::abs(std::real(complex1)) < 0.00001 && std::abs(std::imag(complex1)) < 0.00001){
            return true;
        }
        return false;
    } 

    // check if complex1 == complex2
    bool complexApproxEqual(std::complex<double> complex1, std::complex<double> complex2){
        if (std::abs(std::real(complex1) - std::real(complex2)) < 0.00001 && std::abs(std::imag(complex1) - std::imag(complex2)) < 0.00001){
            return true;
        }
        return false;
    } 

    // check if complex1 == -1 * complex2
    bool complexApproxMinusEqual(std::complex<double> complex1, std::complex<double> complex2){
        if (std::abs(std::real(complex1) + std::real(complex2)) < 0.00001 && std::abs(std::imag(complex1) + std::imag(complex2)) < 0.00001){
            return true;
        }
        return false;
    } 

    // check if complex1 == i * complex2
    bool complexApproxiEqual(std::complex<double> complex1, std::complex<double> complex2){
        if (std::abs(std::real(complex1) + std::imag(complex2)) < 0.00001 && std::abs(std::imag(complex1) - std::real(complex2)) < 0.00001){
            return true;
        }
        return false;
    } 

    // check if complex1 == -i * complex2
    bool complexApproxMinusiEqual(std::complex<double> complex1, std::complex<double> complex2){
        if (std::abs(std::real(complex1) - std::imag(complex2)) < 0.00001 && std::abs(std::imag(complex1) + std::real(complex2)) < 0.00001){
            return true;
        }
        return false;
    } 

    // return -1 * orig
    std::complex<double> multiplyByMinusOne(std::complex<double> orig){
        std::complex<double> res(-1. * std::real(orig), -1. * std::imag(orig));
        return res;
    }

    // return i * orig
    std::complex<double> multiplyByi(std::complex<double> orig){
        std::complex<double> res(-1. * std::imag(orig), std::real(orig));
        return res;
    }

    // return -i * orig
    std::complex<double> multiplyByMinusi(std::complex<double> orig){
        std::complex<double> res(std::imag(orig), -1. * std::real(orig));
        return res;
    }

    // makes a Coset version of a Stabilizergroupvalue (shift = I)
    Coset stab2Coset(StabilizerGroupValue stab){
        Coset stabcoset = Coset();
        for (int i = 0; i < int(stab.size()); i++){
            stabcoset.second.push_back(stab[i]);
        }
        stabcoset.first = LimEntry("I");

        return stabcoset;
    }

    // makes a stabilizergroup from a stabilizergroupvalue (uses pointers instead of directly objects)
    StabilizerGroup toStabilizerGroup(StabilizerGroupValue groupvalue){
        std::vector<LimEntry<>*> group;
        for (int j = 0; j < int(groupvalue.size()); j++){
            group.push_back(new LimEntry<>(groupvalue[j]));
        }
        return group;
    }

    inline void findGeneratingSet(StabilizerGroupValue& stab, Qubit numQubits){
        if (PRINT){
            std::cout << "findGeneratingSet: stab before: {";
            for (int i = 0; i < int(stab.size()); i++){
                std::cout << LimEntry<>::to_string(&(stab[i]), numQubits) << ", ";
            }
            std::cout << "}" << std::endl;
        }
        

        toColumnEchelonForm(stab, numQubits);

        if (stab.size() == 0){
            stab.push_back(LimEntry("I"));
        }

        if (PRINT){
            std::cout << "findGeneratingSet: stab after: {";
            for (int i = 0; i < int(stab.size()); i++){
                std::cout << LimEntry<>::to_string(&(stab[i]), numQubits) << ", ";
            }
            std::cout << "}" << std::endl;
        }

        // StabilizerGroupValue res = StabilizerGroupValue();
        // if (stab.size() == 3){
        //     if (LimEntry<>::to_string(&(stab[0]), numQubits) == "II" && LimEntry<>::to_string(&(stab[1]), numQubits) == "ZZ" && LimEntry<>::to_string(&(stab[2]), numQubits) == "-XX"){
        //         for (int i = 0; i < 3; i++){
        //             res.push_back(stab[i]);
        //         }
        //     }
        // }
        // else if (stab.size() == 1){
        //     if (LimEntry<>::to_string(&(stab[0]), numQubits) == "II"){
        //         res.push_back(stab[0]);
        //     }
        // }
        // return res;
    }

    std::pair<LimEntry<NUM_QUBITS>, bool> findCosetIntersectionElement(Coset coset1, Coset coset2, Qubit numQubits){
        std::vector<LimEntry<NUM_QUBITS>*> G, H;
        toColumnEchelonForm(coset1.second);
        toColumnEchelonForm(coset2.second);
        for (int j = 0; j < int(coset1.second.size()); j++){
            G.push_back(&(coset1.second[j]));
        }
        for (int j = 0; j < int(coset2.second.size()); j++){
            if (LimEntry<>::to_string(&(coset2.second[j]), 0) != "-iI")
                H.push_back(&(coset2.second[j]));
        }
        LimEntry<> helpElt = (coset1.first).getInverse();
        helpElt.multiplyBy(coset2.first);
        std::pair<LimEntry<NUM_QUBITS>, bool> res;
        LimEntry<> I = LimEntry("I");
        res = getCosetIntersectionElementPauli(G, H, &(helpElt), &(I), phase_t::phase_one);
        (res.first).leftMultiplyBy(coset1.first);
        return res;
    }

    Coset findCosetIntersection(Coset coset1, Coset coset2, Qubit numQubits){
        if (PRINT){
            std::cout << "FindCosetIntersection: coset1: (" << LimEntry<>::to_string(&(coset1.first), numQubits) << ", {";
            for (int i = 0; i < int(coset1.second.size()); i++){
                std::cout << LimEntry<>::to_string(&(coset1.second[i]), numQubits) << ", ";
            }
            std::cout << "}),      coset2: (" << LimEntry<>::to_string(&(coset2.first), numQubits) << ", {";
            for (int i = 0; i < int(coset2.second.size()); i++){
                std::cout << LimEntry<>::to_string(&(coset2.second[i]), numQubits) << ", ";
            }
            std::cout << "})" << std::endl;
        }
        Coset res = Coset();
        
        std::pair<LimEntry<NUM_QUBITS>, bool> findElement = findCosetIntersectionElement(coset1, coset2, numQubits);
        if (findElement.second){
            res.first = findElement.first;
            res.second = intersectGroupsPauli(toStabilizerGroup(coset1.second), coset2.second);
            if (res.second.size() == 0){
                res.second.push_back(LimEntry<>("I"));
            }
        }
        
        if (PRINT){
            std::cout << "Result findintersection: (" << LimEntry<>::to_string(&(res.first), numQubits) << ", {";
            for (int i = 0; i < int(res.second.size()); i++){
                std::cout << LimEntry<>::to_string(&(res.second[i]), numQubits) << ", ";
            }
            std::cout << "})" << std::endl;
        }
        return res;
    }

    PauliLIMCoset findSolution(PauliLIMCoset Map1, PauliLIMCoset Map2, std::string pauli, StabilizerGroupValue Stab, Qubit numQubits){
        if (Map1.alpha.noValues || Map2.alpha.noValues){
            // std::cout << "findsolution: either alpha1 or alpha2 = novalues\n";
            return PauliLIMCoset();
        }
        Coset Map1Coset = {Map1.LIM, Map1.stab};
        Coset Map2Coset = {Map2.LIM, Map2.stab};
        std::complex<double> complexZero(0.0, 0.0);
            if (complexApproxEqual(Map1.alpha.value, Map2.alpha.value)){
                std::pair<LimEntry<NUM_QUBITS>, bool> intersection = findCosetIntersectionElement(Map1Coset, Map2Coset, numQubits);
                if (intersection.second){
                    std::string limstring (combinePauliStrings(pauli, LimEntry<>::to_string(&(intersection.first), numQubits)));
                    LimEntry transformation = LimEntry(limstring);
                    PauliLIMCoset solution = PauliLIMCoset();
                    solution.LIM = transformation;
                    solution.stab = Stab;
                    solution.alpha = Map1.alpha;
                    return solution;
                }
            }
            else if (!Map2.alpha.allValues && Map1.alpha.allValues){
                // std::cout << "Alpha1 = allValues\n";
                std::string limstring (combinePauliStrings(pauli, LimEntry<>::to_string(&(Map2Coset.first), numQubits)));
                LimEntry transformation = LimEntry(limstring);
                PauliLIMCoset solution = PauliLIMCoset();
                solution.LIM = transformation;
                solution.stab = Stab;
                solution.alpha = Map2.alpha;
                return solution;
            }
            else if (!Map1.alpha.allValues && Map2.alpha.allValues){
                // std::cout << "Alpha2 = allValues\n";
                std::string limstring (combinePauliStrings(pauli, LimEntry<>::to_string(&(Map1Coset.first), numQubits)));
                LimEntry transformation = LimEntry(limstring);
                PauliLIMCoset solution = PauliLIMCoset();
                solution.LIM = transformation;
                solution.stab = Stab;
                solution.alpha = Map1.alpha;
                return solution;
            }
            else if (Map1.alpha.allValues && Map2.alpha.allValues){
                // std::cout << "Alpha1 = Alpha2 = allValues\n";
                PauliLIMCoset solution = PauliLIMCoset();
                solution.alpha = Map1.alpha;
                return solution;
            }
        return PauliLIMCoset();
    }

    StabilizerGroupValue findStabilizerGroup(std::vector<std::complex<double>> vec){
        if (PRINT){
            std::cout << "findStabGroup: vec = ";
            printVec(vec);
        }

        if (int(vec.size()) == 2){
            StabilizerGroupValue stab = StabilizerGroupValue();
            std::string istring ("I");
            stab.push_back(LimEntry(istring));
            if (complexApproxEqual(vec[0], vec[0]) && complexApproxEqual(vec[1], multiplyByMinusOne(vec[1]))){
                std::string zstring ("Z");
                stab.push_back(LimEntry(zstring));
            }
            if (complexApproxMinusEqual(vec[0], vec[0]) && complexApproxMinusEqual(vec[1], multiplyByMinusOne(vec[1]))){
                std::string zstring ("-Z");
                stab.push_back(LimEntry(zstring));
            }
            if (complexApproxiEqual(vec[0], vec[0]) && complexApproxiEqual(vec[1], multiplyByMinusOne(vec[1]))){
                std::string zstring ("iZ");
                stab.push_back(LimEntry(zstring));
            }
            if (complexApproxMinusiEqual(vec[0], vec[0]) && complexApproxMinusiEqual(vec[1], multiplyByMinusOne(vec[1]))){
                std::string zstring ("-iZ");
                stab.push_back(LimEntry(zstring));
            }
            
            if (complexApproxEqual(vec[0], vec[1])){
                std::string xstring ("X");
                stab.push_back(LimEntry(xstring));
            }
            if (complexApproxMinusEqual(vec[0], vec[1])){
                std::string xstring ("-X");
                stab.push_back(LimEntry(xstring));
            }
            if (complexApproxiEqual(vec[0], vec[1]) && complexApproxiEqual(vec[1], vec[0])){
                std::string xstring ("iX");
                stab.push_back(LimEntry(xstring));
            }
            if (complexApproxMinusiEqual(vec[0], vec[1]) && complexApproxMinusiEqual(vec[1], vec[0])){
                std::string xstring ("-iX");
                stab.push_back(LimEntry(xstring));
            }

            if (complexApproxEqual(vec[0], multiplyByMinusi(vec[1])) && complexApproxEqual(vec[1], multiplyByi(vec[0]))){
                std::string ystring ("Y");
                stab.push_back(LimEntry(ystring));
            }
            if (complexApproxEqual(vec[0], multiplyByi(vec[1])) && complexApproxEqual(vec[1], multiplyByMinusi(vec[0]))){
                std::string ystring ("-Y");
                stab.push_back(LimEntry(ystring));
            }
            if (complexApproxEqual(vec[0], vec[1]) && complexApproxEqual(vec[1], multiplyByMinusOne(vec[0]))){
                std::string ystring ("iY");
                stab.push_back(LimEntry(ystring));
            }
            if (complexApproxEqual(vec[0], multiplyByMinusOne(vec[1])) && complexApproxEqual(vec[1], vec[0])){
                std::string ystring ("-iY");
                stab.push_back(LimEntry(ystring));
            }
            if (PRINT){
                std::cout << "Stab of ";
                printVec(vec);
                for (int i = 0; i < int(stab.size()); i++){
                    std::cout << LimEntry<>::to_string(&(stab[i]), Qubit((vec.size() / 2) - 1)) << ", ";
                }
            }
            findGeneratingSet(stab, 1);
            return stab;
        } //leaf
        
        //not a leaf
        std::vector<std::complex<double>> vec1, vec2, vecz2, vecy1, vecy2;
        for (int i = 0; i < int(vec.size()/2); i++){
            vec1.push_back(vec[i]);
            vec2.push_back(vec[i+int(vec.size()/2)]);
            vecz2.push_back(multiplyByMinusOne(vec[i+int(vec.size()/2)]));
            vecy1.push_back(multiplyByi(vec[i]));
            vecy2.push_back(multiplyByMinusi(vec[i+int(vec.size()/2)]));
        } //make vectors of the top and bottom halves of vec1 and vec2
        // std::cout<< "findStabgroup: \n";
        // printVec(vec1);
        // printVec(vec2);
        // printVec(vecz2);
        // printVec(vecy1);
        // printVec(vecy2);

        bool foundsol;
        StabilizerGroupValue stab1;
        StabilizerGroupValue stab2;
        if (!vecIsZero(vec1)){
            // std::cout << "vec1 is not zero\n";
            stab1 = findStabilizerGroup(vec1);}
        if (!vecIsZero(vec2)){
            // std::cout << "vec2 is not zero\n";
            stab2 = findStabilizerGroup(vec2);}

        if (PRINT){
            std::cout << "Stab1 = ";
            for (int i = 0; i < int(stab1.size()); i++){
                std::cout << LimEntry<>::to_string(&(stab1[i]), Qubit((vec1.size() / 2) - 1)) << ", ";
            }
            std::cout << std::endl << "Stab2 = ";
            for (int i = 0; i < int(stab2.size()); i++){
                std::cout << LimEntry<>::to_string(&(stab2[i]), Qubit((vec2.size() / 2) - 1)) << ", ";
            }
            std::cout << std::endl;
        }

        Coset stab1coset = stab2Coset(stab1);
        Coset stab2coset = stab2Coset(stab2);
        StabilizerGroupValue solution = StabilizerGroupValue();
        StabilizerGroupValue intersectionI;
        if (!vecIsZero(vec1) && !vecIsZero(vec2)){
            StabilizerGroup stab1group = toStabilizerGroup(stab1);
            intersectionI = intersectGroupsPauli(stab1group, stab2);
        }
        else if (vecIsZero(vec1) && !vecIsZero(vec2)){
            intersectionI = stab2;
        }
        else if (!vecIsZero(vec1) && vecIsZero(vec2)){
            intersectionI = stab1;
        }
        for (int i = 0; i < int(intersectionI.size()); i++){
            std::string limstring (combinePauliStrings("I", LimEntry<>::to_string(&(intersectionI[i]), Qubit(findNumQubits(vec1)))));
            solution.push_back(LimEntry(limstring));
            if (PRINT){
                std::cout << "solution na I:";
                for (int i = 0; i < int(solution.size()); i++){
                    std::cout << LimEntry<>::to_string(&(solution[i]), Qubit((vec1.size() / 2))) << ", ";
                }
            }
        }
        PauliLIMCoset res;
        Coset intersectionZ;
        Alpha alpha1, alpha2;
        std::complex<double> complexOne(1.0, 0.0);
        alpha1.value = complexOne;
        alpha2.value = complexOne;
        if (!vecIsZero(vec1) && !vecIsZero(vec2)){
            // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
            res = vecs2Pauli(vecz2, vec2, foundsol);
            Coset MapZ_2 = {res.LIM, res.stab};
            if (PRINT){
                printAlpha(res.alpha);
                std::cout << "MapZ_2 = (" << LimEntry<>::to_string(&(MapZ_2.first), Qubit((vec1.size() / 2))) << ", {";
                for (int i = 0; i < int(MapZ_2.second.size()); i++){
                    std::cout << LimEntry<>::to_string(&(MapZ_2.second[i]), Qubit((vec1.size() / 2))) << ", ";
                }
                std::cout << "})" << std::endl;
            }
            if (!res.alpha.noValues){
                intersectionZ = findCosetIntersection(stab1coset, MapZ_2, int_fast8_t(findNumQubits(vec1)));
                alpha2 = res.alpha;
            }
        }
        else if (vecIsZero(vec1) && !vecIsZero(vec2)){
            // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
            res = vecs2Pauli(vecz2, vec2, foundsol);
            Coset MapZ_2 = {res.LIM, res.stab};
            if (PRINT){
                printAlpha(res.alpha);
                std::cout << "MapZ_2 = (" << LimEntry<>::to_string(&(MapZ_2.first), Qubit((vec1.size() / 2))) << ", {";
                for (int i = 0; i < int(MapZ_2.second.size()); i++){
                    std::cout << LimEntry<>::to_string(&(MapZ_2.second[i]), Qubit((vec1.size() / 2))) << ", ";
                }
                std::cout << "})" << std::endl;
            }
            if (!res.alpha.noValues){
                intersectionZ = MapZ_2;
                alpha2 = res.alpha;
            }
        }
        else if (!vecIsZero(vec1) && vecIsZero(vec2)){
            intersectionZ = stab1coset;
        }
        
        if (int(intersectionZ.second.size()) > 0 && alphaValueEqual(alpha1, alpha2)){
            for (int i = 0; i < int(intersectionZ.second.size()); i++){
                intersectionZ.second[i].leftMultiplyBy(intersectionZ.first);
                std::string limstring (combinePauliStrings("Z", LimEntry<>::to_string(&(intersectionZ.second[i]), Qubit(findNumQubits(vec1)))));
                solution.push_back(LimEntry(limstring));
                if (PRINT){
                    std::cout << "solution na Z:";
                    for (int i = 0; i < int(solution.size()); i++){
                        std::cout << LimEntry<>::to_string(&(solution[i]), Qubit((vec1.size() / 2))) << ", ";
                    }
                }
            }
        }
        
        resetAlpha(alpha1);
        resetAlpha(alpha2);
        Coset intersectionX;
        if (!vecIsZero(vec1) && !vecIsZero(vec2)){
            // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
            res = vecs2Pauli(vec1, vec2, foundsol);
            Coset MapX_1 = {res.LIM, res.stab};
            if (PRINT){
                std::cout << "MapX_1 = "<< LimEntry<>::to_string(&(MapX_1.first), Qubit((vec1.size() / 2))) << ", {";
                for (int i = 0; i < int(MapX_1.second.size()); i++){
                    std::cout << LimEntry<>::to_string(&(MapX_1.second[i]), Qubit((vec1.size() / 2))) << ", ";
                }
                std::cout << "})" << std::endl;
            }
            // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
            if (!res.alpha.noValues){
                alpha1 = res.alpha;
                res = vecs2Pauli(vec2, vec1, foundsol);
                Coset MapX_2 = {res.LIM, res.stab};
                if (PRINT){
                    std::cout << "MapX_2 = "<< LimEntry<>::to_string(&(MapX_2.first), Qubit((vec1.size() / 2))) << ", {";
                    for (int i = 0; i < int(MapX_2.second.size()); i++){
                        std::cout << LimEntry<>::to_string(&(MapX_2.second[i]), Qubit((vec1.size() / 2))) << ", ";
                    }
                    std::cout << "})" << std::endl;
                }
                if (!res.alpha.noValues){
                    alpha2 = res.alpha;
                    intersectionX = findCosetIntersection(MapX_1, MapX_2, int_fast8_t(findNumQubits(vec1)));
                }
            }
        }
        else if (vecIsZero(vec1) && !vecIsZero(vec2)){
            intersectionX = Coset();
        }
        else if (!vecIsZero(vec1) && vecIsZero(vec2)){
            intersectionX = Coset();
        }
        if (int(intersectionX.second.size()) > 0 && alphaValueEqual(alpha1, alpha2)){
            for (int i = 0; i < int(intersectionX.second.size()); i++){
                intersectionX.second[i].leftMultiplyBy(intersectionX.first);
                std::string limstring (combinePauliStrings("X", LimEntry<>::to_string(&(intersectionX.second[i]), Qubit(findNumQubits(vec1)))));
                solution.push_back(LimEntry(limstring));
                if (PRINT){
                    std::cout << "solution na X:";
                    for (int i = 0; i < int(solution.size()); i++){
                        std::cout << LimEntry<>::to_string(&(solution[i]), Qubit((vec1.size() / 2))) << ", ";
                    }
                }
            }
        }

        if (int(intersectionZ.second.size()) == 0 || int(intersectionX.second.size()) == 0){
            resetAlpha(alpha1);
            resetAlpha(alpha2);
            Coset intersectionY;
            if (!vecIsZero(vec1) && !vecIsZero(vec2)){
                // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
                res = vecs2Pauli(vecy1, vec2, foundsol);
                Coset MapY_1 = {res.LIM, res.stab};
                if (PRINT){
                    std::cout << "MapY_1 = "<< LimEntry<>::to_string(&(MapY_1.first), Qubit((vec1.size() / 2))) << ", {";
                    for (int i = 0; i < int(MapY_1.second.size()); i++){
                        std::cout << LimEntry<>::to_string(&(MapY_1.second[i]), Qubit((vec1.size() / 2))) << ", ";
                    }
                    std::cout << "})" << std::endl;
                }
                // std::cout << "Vecs2pauli aanroep vanuit findstabgroup\n";
                if (!res.alpha.noValues){
                    alpha1 = res.alpha;
                    res = vecs2Pauli(vecy2, vec1, foundsol);
                    Coset MapY_2 = {res.LIM, res.stab};
                    if (PRINT){
                        std::cout << "MapY_2 = "<< LimEntry<>::to_string(&(MapY_2.first), Qubit((vec1.size() / 2))) << ", {";
                        for (int i = 0; i < int(MapY_2.second.size()); i++){
                            std::cout << LimEntry<>::to_string(&(MapY_2.second[i]), Qubit((vec1.size() / 2))) << ", ";
                        }
                        std::cout << "})" << std::endl;
                    }
                    if (!res.alpha.noValues){
                        alpha2 = res.alpha;
                        intersectionY = findCosetIntersection(MapY_1, MapY_2, int_fast8_t(findNumQubits(vec1)));
                    }
                }
            }
            else if (vecIsZero(vec1) && !vecIsZero(vec2)){
                intersectionY = Coset();
            }
            else if (!vecIsZero(vec1) && vecIsZero(vec2)){
                intersectionY = Coset();
            }
            if (int(intersectionY.second.size()) > 0 && alphaValueEqual(alpha1, alpha2)){
                for (int i = 0; i < int(intersectionY.second.size()); i++){
                    intersectionY.second[i].leftMultiplyBy(intersectionY.first);
                    std::string limstring (combinePauliStrings("Y", LimEntry<>::to_string(&(intersectionY.second[i]), Qubit(findNumQubits(vec1)))));
                    solution.push_back(LimEntry(limstring));
                    if (PRINT){
                        std::cout << "solution na Y:";
                        for (int i = 0; i < int(solution.size()); i++){
                            std::cout << LimEntry<>::to_string(&(solution[i]), Qubit((vec1.size() / 2))) << ", ";
                        }
                    }
                }
            }
        }

        findGeneratingSet(solution, int_fast8_t(findNumQubits(vec1)));
        if (PRINT){
            std::cout << "findGeneratingSet: stab after: {";
            for (int i = 0; i < int(solution.size()); i++){
                std::cout << LimEntry<>::to_string(&(solution[i]), int_fast8_t((vec.size() / 2) - 1)) << ", ";
            }
            std::cout << "}" << std::endl;
        }

        return solution;
    }

    std::pair<LimEntry<>, Alpha> handleAlpha(Alpha alpha, std::string pauliOp){
        if (alpha.noValues || alpha.allValues){
            return std::pair(LimEntry(pauliOp), alpha);
        }
        if (std::abs(std::imag(alpha.value)) < 0.000001){
            if (std::real(alpha.value) < 0.0){
                alpha.value = multiplyByMinusOne(alpha.value);
                LimEntry<> LIM = LimEntry(pauliOp);
                LIM.multiplyPhaseBy(2);
                return std::pair(LIM, alpha);
            }
            else {
                return std::pair(LimEntry(pauliOp), alpha);
            }
        }
        else if (std::abs(std::real(alpha.value)) < 0.000001){
            if (std::imag(alpha.value) < 0.0){
                alpha.value = multiplyByi(alpha.value);
                LimEntry<> LIM = LimEntry(pauliOp);
                LIM.multiplyPhaseBy(3);
                return std::pair(LIM, alpha);
            }
            else {
                alpha.value = multiplyByMinusi(alpha.value);
                LimEntry<> LIM = LimEntry(pauliOp);
                LIM.multiplyPhaseBy(1);
                return std::pair(LIM, alpha);
            }
        }
        return std::pair(LimEntry(pauliOp), alpha);
    }

    std::pair<LimEntry<>, Alpha> vecs2PauliBase(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2, bool &foundsol){
        foundsol = true;
        Alpha alpha = Alpha();
        if (complexApproxZero(vec1[0]) && complexApproxZero(vec1[1]) && complexApproxZero(vec2[0]) && complexApproxZero(vec2[1])){
            alpha.allValues = true;
            return std::pair(LimEntry("I"), alpha);
        }

        if ((complexApproxZero(vec1[0]) && !complexApproxZero(vec2[0])) || (complexApproxZero(vec1[1]) && !complexApproxZero(vec2[1]))){
        }
        else if (complexApproxZero(vec1[0]) && complexApproxZero(vec2[0])){
            //std::cout << "v0 = 0 en w0 = 0\n";
            if (!complexApproxZero(vec1[1])){
                //std::cout << "v0 = 0 en w0 = 0 en v1 =/= 0\n";
                alpha.value = vec2[1]/vec1[1];
                if (std::real(alpha.value) < 0){
                    //std::cout << "alpha < 0: Z\n";
                    alpha.value = multiplyByMinusOne(alpha.value);
                    return handleAlpha(alpha, "Z");
                }
                else {
                    //std::cout << "alpha >= 0: I\n";
                    return handleAlpha(alpha, "I");
                }
            }
        }
        else {
            // std::cout << "else voor I/Z\n";
            alpha.value = vec2[0]/vec1[0];
            // std::cout << alpha.value << ", " << alpha.value * vec1[1] <<std::endl;
            if (complexApproxEqual(alpha.value * vec1[1], vec2[1])){
                // std::cout << "alpha *v1 = w1: I\n";
                return handleAlpha(alpha, "I");
            }
            if (complexApproxEqual((multiplyByMinusOne(alpha.value) * vec1[1]),vec2[1])){
                // std::cout << "-alpha *v1 = w1: Z\n";
                return handleAlpha(alpha, "Z");
            }
        }

        if ((complexApproxZero(vec1[0]) && !complexApproxZero(vec2[1])) || (complexApproxZero(vec1[1]) && !complexApproxZero(vec2[0]))){
        }
        else if (complexApproxZero(vec1[0]) && complexApproxZero(vec2[1])){
            //std::cout << "v0 = 0 en w1 = 0\n";
            if (!complexApproxZero(vec1[1])){
                //std::cout << "v0 = 0 en w1 = 0 en v1 =/= 0\n";
                alpha.value = vec2[0]/vec1[1];
                if ((std::abs(std::imag(alpha.value))) > 0.000001){
                    //std::cout << "imag(i) =/= 0: Y\n";
                    alpha.value = multiplyByi(alpha.value);
    	            return handleAlpha(alpha, "Y");
                }
                else {
                    //std::cout << "imag(i) = 0: X\n";
                    return handleAlpha(alpha, "X");
                }
            }
        }
        else {
            //std::cout << "else voor X/Y\n";
            alpha.value = vec2[1]/vec1[0];
            if (complexApproxEqual(alpha.value * vec1[1], vec2[0])){
                //std::cout << "alpha * v1 = w0: X\n";
                return handleAlpha(alpha, "X");
            }
            if (complexApproxEqual(multiplyByMinusOne(alpha.value) * vec1[1], vec2[0])){
                //std::cout << "-alpha * v1 = w0: Y\n";
                alpha.value = multiplyByMinusi(alpha.value);
                return handleAlpha(alpha, "Y");
            }
        }

        foundsol = false;
        alpha.noValues = true;
        //std::cout << "niks gevonden\n";
        return std::pair(LimEntry("I"), alpha);
    }

    PauliLIMCoset vecs2Pauli(std::vector<std::complex<double>> vec1, std::vector<std::complex<double>> vec2, bool &foundsol){
        if (PRINT){
            std::cout << "Calling Vecs2Pauli with vectors vec1 = ";
            printVec(vec1);
            std::cout << "and vec2 = ";
            printVec(vec2);
        }
        if (int(vec1.size()) == 2 && int(vec2.size()) == 2){
            PauliLIMCoset S = PauliLIMCoset();
            std::pair<LimEntry<>, Alpha> res = vecs2PauliBase(vec1, vec2, foundsol);
            S.LIM = res.first;
            S.alpha = res.second;
            // printAlpha(S.alpha);
            if (!S.alpha.allValues && !S.alpha.noValues){
                StabilizerGroupValue stab = findStabilizerGroup(vec1);
                S.stab = stab;
            }
            return S;
        } //leaf

        std::vector<std::complex<double>> vec11, vec12, vec21, vec22, vecz2, vecy1, vecy2;
        for (int i = 0; i < int(vec1.size()/2); i++){
            vec11.push_back(vec1[i]); //v_1
            vec12.push_back(vec1[i+int(vec1.size()/2)]); //v_2
            vec21.push_back(vec2[i]); //w_1
            vec22.push_back(vec2[i+int(vec2.size()/2)]); //w_2
            vecz2.push_back(multiplyByMinusOne(vec1[i+int(vec1.size()/2)])); //-v_2
            vecy1.push_back(multiplyByi(vec1[i])); //i*v_1
            vecy2.push_back(multiplyByMinusi(vec1[i+int(vec1.size()/2)])); //-i*v_2
        } //make vectors of the top and bottom halves of vec1 and vec2 necessary for the recursive calls
        
        // std::cout<< "vecs2pauli: \n";
        // printVec(vec11);
        // printVec(vec12);
        // printVec(vec21);
        // printVec(vec22);
        // printVec(vecz2);
        // printVec(vecy1);
        // printVec(vecy2);

        if (vecIsZero(vec1) && !vecIsZero(vec2)){
            if (PRINT)
                std::cout << "vec1 = 0, vec2 != 0\n";
            PauliLIMCoset res = PauliLIMCoset();
            res.alpha.noValues = true;
            foundsol = true;
            return res;
        }
        else if (!vecIsZero(vec1) && vecIsZero(vec2)){
            if (PRINT)
                std::cout << "vec1 != 0, vec2 = 0\n";
            PauliLIMCoset res = PauliLIMCoset();
            res.alpha.value = std::complex<double>(0.0, 0.0);
            // printAlpha(res.alpha);
            foundsol = true;
            return res;
        }
        else if (vecIsZero(vec1) && vecIsZero(vec2)){
            if (PRINT)
                std::cout << "vec1 = 0, vec2 = 0\n";
            PauliLIMCoset res = PauliLIMCoset();
            res.alpha.allValues = true;
            foundsol = true;
            return res;
        }

        StabilizerGroupValue stab = findStabilizerGroup(vec1);
        if (PRINT){
            std::cout << "Stab of vec1 = ";
            for (int i = 0; i < int(stab.size()); i++){
                std::cout << LimEntry<>::to_string(&(stab[i]), Qubit((vec1.size() / 2) - 1)) << ", ";
            }
        }
        foundsol = true;
        PauliLIMCoset MapI_1Z_1 = vecs2Pauli(vec11, vec21, foundsol);
        
        PauliLIMCoset MapI_2 = vecs2Pauli(vec12, vec22, foundsol);
        // printAlpha(MapI_2.alpha);
        PauliLIMCoset solution = findSolution(MapI_1Z_1, MapI_2, "I", stab, int_fast8_t(findNumQubits(vec1)));
        if (int(solution.stab.size()) != 0)
            return solution;

        PauliLIMCoset MapZ_2 = vecs2Pauli(vecz2, vec22, foundsol);
        solution = findSolution(MapI_1Z_1, MapZ_2, "Z", stab, int_fast8_t(findNumQubits(vec1)));
        if (int(solution.stab.size()) != 0)
            return solution;

        PauliLIMCoset MapX_1 = vecs2Pauli(vec11, vec22, foundsol);
        PauliLIMCoset MapX_2 = vecs2Pauli(vec12, vec21, foundsol);
        solution = findSolution(MapX_1, MapX_2, "X", stab, int_fast8_t(findNumQubits(vec1)));
        if (int(solution.stab.size()) != 0)
            return solution;
        
        PauliLIMCoset MapY_1 = vecs2Pauli(vecy1, vec22, foundsol);
        PauliLIMCoset MapY_2 = vecs2Pauli(vecy2, vec21, foundsol);
        solution = findSolution(MapY_1, MapY_2, "Y", stab, int_fast8_t(findNumQubits(vec1)));
        if (int(solution.stab.size()) != 0)
            return solution;

        foundsol = false;
        solution.alpha.noValues = true;
        return solution;
    }

}
