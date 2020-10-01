// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEETESTUTILS_HPP
#define AXOM_KLEETESTUTILS_HPP

#include <array>

#include "axom/core/numerics/Matrix.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "gmock/gmock.h"

namespace axom {namespace klee { namespace test {

/**
 * Create a 3D affine transformation matrix
 *
 * \param values the values of the matrix, in row-major order
 * \return the affine transformation matrix represented by the rows
 */
numerics::Matrix<double> affine(const std::array<double, 12> &values);

/**
 * Make a vector with the given values.
 *
 * \tparam dims the number of dimensions
 * \param values the values
 * \return a vector with the given values
 */
template<int dims>
primal::Vector<double, dims> makeVector(const double (&values)[dims]) {
    return primal::Vector<double, dims>{values};
}

/**
 * Make a point with the given values.
 *
 * \tparam dims the number of dimensions
 * \param values the values
 * \return a point with the given values
 */
template<int dims>
primal::Point<double, dims> makePoint(const double (&values)[dims]) {
    return primal::Point<double, dims>{values};
}

class MockOperator : public GeometryOperator {
public:
    MOCK_METHOD(Dimensions, startDims, (), (const));
    MOCK_METHOD(Dimensions, endDims, (), (const));
    MOCK_METHOD(void, accept, (GeometryOperatorVisitor &), (const));
};

}}}

#endif //AXOM_KLEETESTUTILS_HPP
