// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file primal_bezier_curve.cpp
 * /brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/primal/geometry/Transform.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/core/numerics/matvecops.hpp"

namespace primal = axom::primal;
namespace numerics = axom::numerics;

namespace
{
/*!
 * \brief Checks if the given two vectors are equal.
 *
 * \param [in] u rhs vector to compare
 * \param [in] v lhs vector to compare
 * \param [in] N the size of the vector
 */
void expect_array_eq(const double* u, const double* v, int N)
{
  for(int i = 0; i < N; ++i)
  {
    EXPECT_DOUBLE_EQ(u[i], v[i]);
  }
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
TEST(primal_transform, twoDconstructor)
{
  constexpr int DIM = 2;
  constexpr int ARRAY = 9;
  using ValueType = double;
  using TransformType = primal::Transform<ValueType, DIM>;
  using VectorType = primal::Vector<ValueType, DIM>;

  {
    SCOPED_TRACE("Default ctor gives identity transform");
    double array_id[ARRAY] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    TransformType idtr;
    expect_array_eq(idtr.getTransform().data(), array_id, ARRAY);
  }

  {
    SCOPED_TRACE("Ctor from array");
    double array_matrix[ARRAY] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
    TransformType arraytr(array_matrix);
    expect_array_eq(arraytr.getTransform().data(), array_matrix, ARRAY);
  }

  {
    SCOPED_TRACE("Translation ctor");
    double array_translate[ARRAY] = {1, 0, 0, 0, 1, 0, 3.1, -2.4, 1};
    double vector_data[DIM] = {3.1, -2.4};
    VectorType translate(vector_data, ARRAY);
    TransformType translatetr(translate);
    expect_array_eq(translatetr.getTransform().data(), array_translate, ARRAY);
  }

  {
    SCOPED_TRACE("Rotation ctor");
    double theta = M_PI / 3.;
    double s = sin(theta);
    double c = cos(theta);
    double array_rotate[ARRAY] {c, s, 0, -s, c, 0, 0, 0, 1};
    TransformType rotatetr(0, theta);
    expect_array_eq(rotatetr.getTransform().data(), array_rotate, ARRAY);
  }
}

//------------------------------------------------------------------------------
TEST(primal_transform, threeDconstructor)
{
  constexpr int DIM = 3;
  constexpr int ARRAY = 16;
  using ValueType = double;
  using TransformType = primal::Transform<ValueType, DIM>;
  using VectorType = primal::Vector<ValueType, DIM>;

  double array_id[ARRAY] =
     {1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1};
  TransformType idtr;
  expect_array_eq(idtr.getTransform().data(), array_id, ARRAY);

  double array_matrix[ARRAY] =
     {0, 0, -1, 0,
      0, 1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, 1};
  TransformType arraytr(array_matrix);
  expect_array_eq(arraytr.getTransform().data(), array_matrix, ARRAY);

  double array_translate[ARRAY] =
     {1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      3.1, -2.4, 0.6, 1};
  double vector_data[DIM] = {3.1, -2.4, 0.6};
  VectorType translate(vector_data, ARRAY);
  TransformType translatetr(translate);
  expect_array_eq(translatetr.getTransform().data(), array_translate, ARRAY);

  double theta = M_PI / 3.;
  double s = sin(theta);
  double c = cos(theta);
  double array_rotate[ARRAY] =
     {c, s, 0, 0,
      -s, c, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1};
  TransformType rotatetr(2, theta);
  expect_array_eq(rotatetr.getTransform().data(), array_rotate, ARRAY);
}



//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();

  return result;
}
