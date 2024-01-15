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
/// C++ function to compute a Gauss-Legendre grid. 

#include <boost/math/special_functions/legendre.hpp>
#include <vector>

/// Computes the gauss legendre quadrature
/// @param[in]     n - the number of points
/// @param[in,out] x - the points for the quadrature
/// @param[in,out] weight - the weights for the quadrature     
extern "C" void double_compute_gauss_legendre(const int n, double* x, double* weight) {
    
    // Compute the zeros and copy to the array (only for x >= 0.0)
    auto zeros = boost::math::legendre_p_zeros<double>(n);
    zeros.reserve(n);
    for (auto &z : zeros) {
        if (z > 1.0e-8) zeros.push_back(-z);
    }
    // Efficient copy through direct copy
    std::memcpy(x,zeros.data(),n*sizeof(double));
    
    // Compute the weights
    for (int i = 0; i < n; i++) {
        auto z = x[i];
        // This function computes the derivative of the Legendre Polynomial of n-order at z
        auto p = boost::math::legendre_p_prime(n,z);
        weight[i] = 2 / ((1 - z * z) * p * p);
    }

    return;
}