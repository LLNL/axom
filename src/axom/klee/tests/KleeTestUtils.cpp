// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/tests/KleeTestUtils.hpp"

namespace axom {namespace klee { namespace test {

numerics::Matrix<double> affine(const std::array<double, 12> &values) {
    numerics::Matrix<double> m(4, 4);
    m(3, 3) = 1;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            m(i, j) = values[i * 4 + j];
        }
    }
    return m;
}

}}}
