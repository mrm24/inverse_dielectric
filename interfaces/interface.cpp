// Copyright 2023 EXCITING developers
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

/// @file
/// This file contains a module for exposing the inverse_dielectric_t to C, and python

#include <memory>
#include <complex>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>

namespace py = pybind11;

/// External functions to call Fortran stuff
extern "C" {
    extern void allocate_inverse_dielectric_t(void** object_ptr);
    extern void deallocate_inverse_dielectric_t(void** object_ptr);
    extern void init_common(void** object_ptr, double* lattice, double* redpos, long int* elements, long nq[3]);
    extern void set_dielectric_blocks(void** object_ptr, std::complex<double>* h, std::complex<double>* wl, std::complex<double>* wu, std::complex<double>* ib);
    extern void compute_anisotropic_avg(void** object_ptr, bool hermitian);
    extern void invert_body(void** object_ptr, std::complex<double>* body);
    extern long int get_n_basis(void** object_ptr);
    extern void get_inverted_blocks(void** object_ptr, std::complex<double>* inv_head, std::complex<double>* inv_wingL, std::complex<double>* inv_wingU, std::complex<double>* inv_body);
    extern void nullify_body(void** object_ptr);
}

class inverse_dielectric_cxx {

private:
    void* inverse_dielectric_f90;

public:

    // Constructor
    inverse_dielectric_cxx() {
        allocate_inverse_dielectric_t(&inverse_dielectric_f90);
    }

    // Destructor
    ~inverse_dielectric_cxx() {
        deallocate_inverse_dielectric_t(&inverse_dielectric_f90);
    }

    // Declare copy and movement constructor deleted
    inverse_dielectric_cxx(const inverse_dielectric_cxx& A) = delete;
    inverse_dielectric_cxx(inverse_dielectric_cxx&& A)      = delete;

    void destroy(){
        deallocate_inverse_dielectric_t(&inverse_dielectric_f90);
    }

    void initialize(double* lattice, double* redpos, long* elements, long nq[3]) {
        init_common(&inverse_dielectric_f90, lattice, redpos, elements, nq);
    }

    void setDielectricBlocksFull(std::complex<double>* h, std::complex<double>* wl,
                             std::complex<double>* wu, std::complex<double>* ib) {
        set_dielectric_blocks(&inverse_dielectric_f90, h, wl, wu, ib);
    }

    void setDielectricBlocksPartial(std::complex<double>* h, std::complex<double>* wl,
                             std::complex<double>* wu) {
        set_dielectric_blocks(&inverse_dielectric_f90, h, wl, wu, nullptr);
    }

    void computeAnisotropicAvg(bool hermitian) {
        compute_anisotropic_avg(&inverse_dielectric_f90, hermitian);
    }

    void invertBody(std::complex<double>* body) {
        invert_body(&inverse_dielectric_f90, body);
    }

    long int getNBasis(){
        get_n_basis(&inverse_dielectric_f90);
    }

    void getInvertedBlocks(std::complex<double>* inv_head, std::complex<double>* inv_wingL, 
                           std::complex<double>* inv_wingU, std::complex<double>* inv_body) {
        get_inverted_blocks(&inverse_dielectric_f90, inv_head, inv_wingL, inv_wingU, inv_body);
    }

    void NullifyBody(){
        nullify_body(&inverse_dielectric_f90);
    }
};


PYBIND11_MODULE(InverseDielectric, m) {

    m.doc() = "Python bindings for InverseDielectric library";

    py::class_<inverse_dielectric_cxx>(m, "InverseDielectric")
        .def(py::init<>())
        .def("initialize", &inverse_dielectric_cxx::initialize)
        .def("setDielectricBlocks", &inverse_dielectric_cxx::setDielectricBlocksFull)
        .def("setDielectricBlocks", &inverse_dielectric_cxx::setDielectricBlocksPartial)
        .def("computeAnisotropicAvg", &inverse_dielectric_cxx::computeAnisotropicAvg)
        .def("invertBody", &inverse_dielectric_cxx::invertBody)
        .def("getNBasis", &inverse_dielectric_cxx::getNBasis)
        .def("getInvertedBlocks", &inverse_dielectric_cxx::getInvertedBlocks)
        .def("NullifyBody", &inverse_dielectric_cxx::NullifyBody)
        .def("destroy", &inverse_dielectric_cxx::destroy); 
};