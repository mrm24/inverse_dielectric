# IDieL

**IDieL**(**I**nverse **Diel**ectric **L**ibrary) is a modern Fortran library providing a common framework for GW codes, specifically for the most common operations involving the inverse dielectric function. The code is parallelized via multithreading (OPENMP) and GPUs (via MAGMA library). 

**Supported operations:**
* _Anisotropic averaging of ε<sup>-1</sup>(**q**→**0**) (3D)_: the library contains functions for the computation of the inverse dielectric function at **q**=**Γ** for 3D systems. This last by numerically computing the anisotropic average in a small region (i.e. a small version of the reciprocal cell) centered around  **Γ** using a 131<sup>st</sup>-order Lebedev grid.
*_Anisotropic averaging for **q**=**Γ** (2D)_*:
* _Body inversion_: the code provides functions for inverting the body of the dielectric matrix either using CPUs or GPUs.

**Experimental:**
* _Interfaces_: several major methods are exposed to C++, and to Python (via pybind11). Notice, that the ordering of arrays must be the Fortran one.

# Install:

The programm can be easily installed using CMake (version >= 3.20).

**Supported compilers and libraries versions**

The library has been compiled and tested using: ifort 2021.3.0, icc 2021.3.0, pybind 2.11.1, Python 3.7.10, (MAGMA 2.7.2), spglib 2.1.0, and boost 1.8.4

**CMake Options:**

* _GPU_: Set to ON to allow for GPU acceleration
* _MAGMA_DIR_: If _GPU_ is ON, this option provides the path to MAGMA directory (installed with make)
* _INTERFACES_: If ON, create interface for C++ and Python (this last using pybind11)  
* _PYBIND_CMAKE_DIR_: If _INTERFACES_ is ON, this provides the path to _PYBIND11_ CMake file; in case you installed it via pip (i.e. pip3 install --user pybind11
; pip3 install "pybind11[global]" --user) it will be located in "~/.local/share/cmake/pybind11".

**Testing:**

The test suit is provided via CTest, to run the test cases, go to test folder and execute "ctest".

**Documentation:**

The documentation is build using Doxygen. To compile it run "make doc".
