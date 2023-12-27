# Vecs2pauli: library for converting quantum state vectors to local Pauli transformations and back

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Vecs2pauli is a Python/C++ software library for obtaining all local Pauli symmetries of a vector which represents a quantum state of qubits.
It also contains a bunch of related useful functionalities, such as finding all local Pauli operators mapping one vector to another, intersecting stabiliser groups, and more.

Vecs2pauli is useful for various reasons:
- **fast classical simulation of quantum circuits**: converting from the slow vector-based simulation to the fast stabilizer-based simulation, if possible
- **educating quantum error correction**: students can create quantum error-correction codes in the state-vector picture instead of the less beginner-friendly stabilizer-picture. Specifically, they can start by defining the logical states as vectors, and subsequently use Vecs2Pauli to ...
    - ...find the stabilizers of the code space (the usual representation of a code) 
	- ...find errors which the code cannot correct, i.e. Pauli strings which map one vector in the span to another
- **testing hypotheses** on quantum states for only a few qubits
- (more to be added)

For more information, see the [project website](https://vecs2pauli.github.io/).

# Examples

Vecs2pauli can find a generating set for the Pauli stabiliser group of the two-qubit state `|00> + |11>` (unnormalised for the sake of simplicity):

```python
import vecs2pauli as vtp
myvector = [1, 0, 0, 1]  # the |00> + |11> state
print(vtp.get_stabilizers(myvector))
# prints ["+ZZ", "+ZI"]
```

It can also be used for finding all single-qubit Pauli operations between two vectors:

```python
import vecs2pauli as vtp
source = [1, 0, 0, 1]  # the |00> + |11> state
target = [0, -1j, -1j, 0]  # the |01> + |10> state, multiplied by a global factor -i
print(vtp.get_local_pauli_transformations(source, target))
# prints ((1-0j), '-iIX', ['ZZ', 'XX']), implying all local-Pauli transformations
# from source to target can be written as -i * IX * g where g is II or a product of
# a subset of ["+ZZ", "+XX"]
```

Another example is finding the stabiliser formulation of a [quantum error-correction code](https://en.wikipedia.org/wiki/Quantum_error_correction):

```python
import vecs2pauli as vtp

# the logical basis states |0>_L and |1>_L of the 3-qubit repetition code
logical_zero = [1, 0, 0, 0, 0, 0, 0, 0]  # |000>
logical_one = [0, 0, 0, 0, 0, 0, 0, 1]  # |111>

# finding the stabilizers of the code
stabilizer_group_zero = vtp.get_stabilizers(logical_zero)  # {ZII, IZI, IIZ}
stabilizer_group_one = vtp.get_stabilizers(logical_one)  # {-ZII, -IZI, -IIZ}

print(vtp.intersect_stabilizer_groups(stabilizer_group_zero, stabilizer_group_one))
# prints ["+ZZI", "+IZZ"]
```

For more examples, see [the vecs2pauli website](https://vecs2pauli.github.io/examples.html).


# Installation

Building vecs2pauli requires that you have installed:

* Python 3
* A C++17-compliant compiler, such as `g++`
* CMake (version `3.9` or higher)

Before installation, you might want to make a [Python virtual environment](https://docs.python.org/3/tutorial/venv.html) in order for the vecs2pauli installation not to interfere with other installed python packages.

First, clone the repository:
```
git clone --recursive [url-to-clone-this-repository]
```

The following instructions are for Linux. For Windows, replace `python3` by `py`.

First, ensure that your current working directory is the top-level directory of the freshly cloned repository (the directory is probably named `vecs2pauli`).
The recommended means of installation is through the package manager [pip](https://pip.pypa.io/en/stable/):

- For Linux: `python3 -m pip install .`
- For Windows: `py -m pip install .`

Now, you can use vecs2pauli as part of a python script.
For example, create the text file `example.py` with the content

```
import vecs2pauli as vtp
myvector = [1, 0, 0, 0]  # the |00> state
print(vtp.get_stabilizers(myvector))
```
and run it (`python3 example.py`) which should print `["ZZ", "ZI"]`.


**Alternative route using `setup.py`**

Alternatively to the installation instructions above, one can directly install using `setup.py`.
This first requires one to install all packages mentioned in the file `requirements.txt`: `build`, `scikit-build` and `pybind11`. This is most easily done using `pip` by `python3 -m pip install -r requirements.txt`.
Next, use:

```
python3 setup.py build
python3 setup.py install
```

(As alternative to `python3 setup.py install`: the command `python3 setup.py build` results in a wheel file in the folder `dist` which can be directly installed using `python3 -m pip vecs2pauli-xxx.whl` where `xxx` is a string that is specific to your operating system.)




# Testing vecs2pauli

The Python test suite can be run by first `pip`-installing the Python package
and then running `pytest` from the top-level directory:

```
python -m pip install .
pytest
```


# Advanced: obtaining the C++ build artifacts

The installation procedure above suffices to use vecs2pauli as Python package.
To directly access the C++ build, first ensure that you have C++17-compliant compiler installed and CMake installed (one way to do the latter, is through `pip`: `python3 -m pip cmake`).

Again, ensure that your current working directory is the top-level directory of the freshly cloned repository (the directory which is probably named `vecs2pauli`).
Within this folder, create a directory named `build` (in Linux, this is done by `mkdir build`) and then run

```
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

If you get the error

```
Could not find a package configuration file provided by "pybind11" with any of the following names:

    pybind11Config.cmake
    pybind11-config.cmake
```
then install pybind11 using pip and add 

- for Linux: `export CMAKE_PREFIX_PATH=mypythondir/lib/python3.11/site-packages/pybind11/share/cmake/pybind11/`
- for Windows: `set CMAKE_PREFIX_PATH=mypythondir\Lib\site-packages\pybind11\share\cmake\pybind11\`

where `mypythondir` is the full path to where your python is installed, for example your virtual environment directory.

The resulting binary can be found in `build/app/Debug/` and can be called in Windows as `\vecs2pauli_app.exe 1 0 0 0 0 0 1 0` to obtain for example all Pauli operators mapping the vector (1 + 0j, 0 + 0j) to (0 + 0j, 1 + 0j).

The build process can be customized with the following CMake variables, which can be set by adding `-D<var>={ON, OFF}` to the `cmake` call:

* `BUILD_TESTING`: Enable building of the test suite (default: `ON`)
* `BUILD_PYTHON`: Enable building the Python bindings (default: `ON`)

For testing, (i.e. when build with `-DBUILD_TESTING=ON`), the C++ test suite of `vecs2pauli` can be run using `ctest` from the build directory:

```
cd build
ctest
```


# Documentation

To get started, see 'Installation' and 'Examples' above. More examples can be found in the 'examples' folder.

Building a full documentation of the API is work in progress. To build and see the current status of the documentation, go to the `docs` folder, install all packages in `requirements.txt` using pip (in Linux, the command is `pip3 install -r requirements.txt`), then run `make html`. Finally, open the file `index.html` in the `_build/html` folder using your web browser.

# Contributors

- Lieke Vertegaal (main researcher)
- Tim Coopmans (project lead and [point of contact](https://www.universiteitleiden.nl/en/staffmembers/tim-coopmans))

The Pauli coset intersection functionality is taken from [previous work with part of the authors involved](https://quantum-journal.org/papers/q-2023-09-11-1108/): `vecs2pauli` calls (a fork of) [its implementation](https://github.com/cda-tum/dd_package/tree/limdd) .
