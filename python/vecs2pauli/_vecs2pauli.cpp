#include "vecs2pauli.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>
#include <utility>
#include <complex>
#include <fstream>
#ifdef __cplusplus
extern "C" {
#endif
	#include "binarylinalg/include/extend_check_matrix.h"
#ifdef __cplusplus
}
#endif

namespace py = pybind11;

namespace dd {

void printOutput(dd::PauliLIMCoset res, std::vector<std::complex<double>> vec1){
    std::cout << "Result: PauliLIMCoset = (";
    dd::printAlpha(res.alpha);
    std::cout << ", " << dd::LimEntry<>::to_string(&(res.LIM), dd::findNumQubits(vec1) - 1) << ", {";
    if (int(res.stab.size() > 1)){
        for (int i = 0; i < int(res.stab.size()-1); i++){
            std::cout << dd::LimEntry<>::to_string(&(res.stab[i]), dd::findNumQubits(vec1) - 1) << ", ";
        }
    }
    if (int(res.stab.size() > 0)){
        std::cout << dd::LimEntry<>::to_string(&(res.stab[res.stab.size()-1]), dd::findNumQubits(vec1) - 1);
    }
    std::cout << "})" << std::endl << std::endl;
}

void print_vector(const std::vector<std::complex<double>> &v) {
    dd::printVec(v);
}

std::vector<std::string> intersectGroups(std::vector<std::string> vec1, std::vector<std::string> vec2){
    dd::StabilizerGroupValue group1, group2;
    int numqubits = vec1[0].length();
    if (vec1[0][0] == '-' || vec1[0][0] == 'i'){
        if (vec1[0][1] == 'i'){
            numqubits--;
        }
        numqubits--;
    }
    for (int i = 0; i < vec1.size(); i++){
        dd::LimEntry<> elt(vec1[i]);
        group1.push_back(elt);
    }
    for (int i = 0; i < vec2.size(); i++){
        dd::LimEntry<> elt(vec2[i]);
        group2.push_back(elt);
    }
    dd::StabilizerGroupValue res = dd::intersectGroupsPauli(dd::toStabilizerGroup(group1), group2);
    std::vector<std::string> intersectedGroup;
    for (int i = 0; i < int(res.size()); i++){
        intersectedGroup.push_back(dd::LimEntry<>::to_string(&(res[i]), numqubits - 1));
    }
    return intersectedGroup;
}

//void intersectCosets(std::complex<double> factor_a, std::string coset_repr_a, std::vector<std::string> group_a, std::complex<double> factor_b, std::string coset_repr_b, std::vector<std::string> group_b){
std::tuple<std::string, std::vector<std::string>> intersectCosets(std::complex<double> factor_a, std::string coset_repr_a, std::vector<std::string> group_a, std::complex<double> factor_b, std::string coset_repr_b, std::vector<std::string> group_b){
//std::tuple<std::complex<double>, std::string, std::vector<std::string>, bool, bool> intersectCosets(std::complex<double> factor_a, std::string coset_repr_a, std::vector<std::string> group_a, std::complex<double> factor_b, std::string coset_repr_b, std::vector<std::string> group_b){

    int num_qubits = coset_repr_a.length();

    //make a Coset object for the first coset
    dd::Coset coset_a = dd::Coset();
    for(int i = 0; i < int(group_a.size()); i++){
    	coset_a.second.push_back(dd::LimEntry(group_a[i]));
	}
    coset_a.first = dd::LimEntry(coset_repr_a);

    //make a Coset object for the second coset
    dd::Coset coset_b = dd::Coset();
    for(int i = 0; i < int(group_b.size()); i++){
    	coset_b.second.push_back(dd::LimEntry(group_b[i]));
	}
    coset_b.first = dd::LimEntry(coset_repr_b);


    dd::Coset res = dd::findCosetIntersection(coset_a, coset_b, num_qubits);

    //bool foundsol = true;
    //dd::PauliLIMCoset res = dd::vecs2Pauli(v, w, foundsol);
    // printOutput(res, v);

    // create output
    //std::tuple<std::complex<double>, std::string, std::vector<std::string>, bool, bool> output;
    //std::tuple<std::string, std::vector<std::string>> output;
    std::tuple<std::string, std::vector<std::string>> output;
    std::vector<std::string> stabgens;
    
    for (int i = 0; i < int(res.second.size()); i++){
        stabgens.push_back(dd::LimEntry<>::to_string(&res.second[i], num_qubits)); //dd::findNumQubits(v) - 1
        //stabgens.push_back(dd::LimEntry<>::to_string(&(res.second[i])), num_qubits - 1); //dd::findNumQubits(v) - 1));
    }
    //output = {res.alpha.value, dd::LimEntry<>::to_string(&(res.LIM), num_qubits), stabgens, res.alpha.allValues, res.alpha.noValues};
    //output = {dd::LimEntry<>::to_string(&res.LIM), stabgens};
    output = {dd::LimEntry<>::to_string(&res.first, num_qubits), stabgens};

    /*
    dd::PauliLLIMCoset pcoset_a;
    pcoset_a.LIM = coset_a.first;
    pcoset_a.stab = coset_a.second;
    pcoset.alpha = factor_a;
    */

    return output;
}


std::vector<std::string> makeFullGroup(std::vector<std::string> gens){
    std::vector<std::string> fullGroup;
    if (gens.size() == 0){
        return gens;
    }
    int numqubits = gens[0].length();
    if (gens[0][0] == '-' || gens[0][0] == 'i'){
        if (gens[0][1] == 'i'){
            numqubits--;
        }
        numqubits--;
    }
    for (int i = 0; i < gens.size(); i++){
        fullGroup.push_back(gens[i]);
    }
    bool newFound = true;
    while (newFound){
        newFound = false;
        int oldSize = fullGroup.size();
        for (int i = 0; i < oldSize; i++){
            dd::LimEntry<> elt(fullGroup[i]);
            dd::LimEntry<> elt2(fullGroup[i]);
            int newOldSize = fullGroup.size();
            for (int j = 0; j < newOldSize; j++){
                dd::LimEntry<> other(fullGroup[j]);
                elt.multiplyBy(other);
                elt2.leftMultiplyBy(other);
                std::string new1 = dd::LimEntry<>::to_string(&(elt), numqubits - 1);
                std::string new2 = dd::LimEntry<>::to_string(&(elt2), numqubits - 1);
                if(std::find(fullGroup.begin(), fullGroup.end(), new1) == fullGroup.end()) {
                    fullGroup.push_back(new1);
                    newFound = true;
                }
                if(std::find(fullGroup.begin(), fullGroup.end(), new2) == fullGroup.end()) {
                    fullGroup.push_back(new2);
                    newFound = true;
                }
            }
        }
    }
    return fullGroup;

}


std::tuple<std::complex<double>, std::string, std::vector<std::string>, bool, bool> runvecs2pauli(const std::vector<std::complex<double>> &v, const std::vector<std::complex<double>> &w){
    if (v.size() != w.size()){
        std::cout << "Please input vectors of equal length\n";
        return {};
    }
    if (!dd::powerOfTwo(int(v.size()))){
        std::cout << "Please input vectors of length 2^n\n";
        return {};
    }
    bool foundsol = true;
    dd::PauliLIMCoset res = dd::vecs2Pauli(v, w, foundsol);
    // printOutput(res, v);
    std::tuple<std::complex<double>, std::string, std::vector<std::string>, bool, bool> output;
    std::vector<std::string> stabgens;
    for (int i = 0; i < int(res.stab.size()); i++){
        stabgens.push_back(dd::LimEntry<>::to_string(&(res.stab[i]), dd::findNumQubits(v) - 1));
    }
    output = {res.alpha.value, dd::LimEntry<>::to_string(&(res.LIM), dd::findNumQubits(v) - 1), stabgens, res.alpha.allValues, res.alpha.noValues};
    return output;
}


std::vector<std::string> runstabgroup(const std::vector<std::complex<double>> &v){
    std::vector<std::string> stabgens;
    if (!dd::powerOfTwo(int(v.size()))){
        std::cout << "Please input vector of length 2^n\n";
        return stabgens; 
    }
    dd::StabilizerGroupValue stab = dd::findStabilizerGroup(v);
    for (int i = 0; i < int(stab.size()); i++){
        stabgens.push_back(dd::LimEntry<>::to_string(&(stab[i]), dd::findNumQubits(v) - 1));
    }
    return stabgens;
}

py::array_t<double> extend_check_matrix(py::array_t<bool> check_matrix){
	  // check input dimensions
	  if ( check_matrix.ndim()     != 2 )
	  {
	    throw std::runtime_error("Input should be 2-D NumPy array");	
	  }
	  size_t num_rows = check_matrix.shape()[0];
	  size_t num_columns = check_matrix.shape()[1];
	  size_t num_variables = num_columns - 1;
	  assert(num_variables % 2 == 0);
	  size_t num_qubits = (size_t) (num_variables / 2);

	  py::buffer_info buf = check_matrix.request();
	  bool* ptr = (bool*) buf.ptr;

	
	  //allocate
        bool** extended_check_matrix = (bool**) malloc(sizeof(bool*) * num_qubits);
        bool* row;
    for(size_t i=0;i<num_qubits;i++){
    	row = (bool*) malloc(sizeof(bool) * (num_variables));
	*(extended_check_matrix + i) = &(row[0]);
    }

    //initialise the top part of extended_check_matrix to equal check_matrix; the rest of the rows only contain zeroes
    for(size_t i=0;i<num_rows;i++){
    	for(size_t j=0;j<num_variables;j++){
	    *(*(extended_check_matrix + i) + j) = ptr[i * num_columns + j];
	}
    }
    for(size_t i=num_rows;i<num_qubits;i++){
    	for(size_t j=0;j<num_variables;j++){
	    *(*(extended_check_matrix + i) + j) = false;
	}
    }

    //extend the check matrix
    extend_rref_check_matrix(num_rows, num_qubits, extended_check_matrix);

        constexpr size_t elsize = sizeof(bool);
        size_t shape[2]{num_qubits, num_columns};
        size_t strides[2]{num_columns * elsize,elsize};
        auto a = py::array_t<bool>(shape, strides);
        auto view = a.mutable_unchecked<2>();

        for(size_t i = 0; i < a.shape(0); i++)
        {
          for(size_t j = 0; j < a.shape(1); j++)
          {
              view(i,j) = *(*(extended_check_matrix + i) + j);
          }
	  if(i < num_rows){
	  	view(i, num_columns - 1) = ptr[i * num_columns + num_columns - 1];
	  }
	  else{
	  	view(i, num_columns - 1) = false;
	  }
        }

	//clean up
    for(size_t i=0;i<num_qubits;i++){
	    free(extended_check_matrix[i]);
    }


        return a;
}

}

namespace vecs2pauli{

PYBIND11_MODULE(_vecs2pauli, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("printOutput", &dd::printOutput, "Printing function for vecs2Pauli result");
    //m.def("print_vector", &dd::print_vector, "prints a vector");
    m.def("_get_local_pauli_transformations", &dd::runvecs2pauli, "run vecs2pauli");
    m.def("_get_stabilizers", &dd::runstabgroup, "run findStabilizerGroup");
    m.def("_expand_generators", &dd::makeFullGroup, "generate full group from generators");
    m.def("_intersect_stabilizer_groups", &dd::intersectGroups, "intersect 2 groups of Pauli strings");
    m.def("_extend_check_matrix", &dd::extend_check_matrix, "extend a nonmaximal check matrix to a maximal one");
    m.def("_intersect_cosets", &dd::intersectCosets, "intersect cosets");
}

} // namespace vecs2pauli
