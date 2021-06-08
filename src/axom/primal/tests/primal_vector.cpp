// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Vector.hpp"

using namespace axom;

//------------------------------------------------------------------------------
TEST(primal_vector, vector_constructors)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef primal::NumericArray<CoordType, DIM> QArray;
  typedef primal::Vector<CoordType, DIM> QVec;

  QVec vec1;
  EXPECT_EQ(vec1.dimension(), DIM);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec1[i], CoordType());
  }

  CoordType val = 5.;
  QVec vec2(val);
  EXPECT_EQ(vec2.dimension(), DIM);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec2[i], val);
  }

  CoordType val3 = 5.;
  int halfPt = DIM / 2;
  QVec vec3(val3, halfPt);
  for(int i = 0; i < DIM; ++i)
  {
    CoordType expVal = i < halfPt ? val3 : CoordType();
    EXPECT_EQ(vec3[i], expVal);
  }

  // Compare vector initialized from arbitrary
  CoordType valsArr[DIM] = {12., 23., 34., 45., 56.432};
  QVec arrVec(valsArr);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(arrVec[i], valsArr[i]);
  }

  QVec arrHalfVec(valsArr, halfPt);
  for(int i = 0; i < DIM; ++i)
  {
    CoordType expVal = i < halfPt ? valsArr[i] : CoordType();
    EXPECT_EQ(arrHalfVec[i], expVal);
  }

  // Construct from numeric array
  QArray arr(valsArr);
  QVec vFromNA(arr);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vFromNA[i], valsArr[i]);
  }

  // Initializer list constructor
  primal::Vector<int, 3> fromInitializerListRightSize = {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListRightSize[i]);
  }

  primal::Vector<int, 3> fromInitializerListTooLong = {10, 20, 30, 40};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListTooLong[i]);
  }

  primal::Vector<int, 5> fromInitializerListTooShort = {10, 20};
  for(int i = 0; i < 2; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerListTooShort[i]);
  }
  for(int i = 2; i < 5; ++i)
  {
    EXPECT_EQ(0, fromInitializerListTooShort[i]);
  }

  primal::Vector<int, 3> fromInitializerNoEqualsSign {10, 20, 30};
  for(int i = 0; i < 3; ++i)
  {
    EXPECT_EQ(10 * (i + 1), fromInitializerNoEqualsSign[i]);
  }
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_from_points_constructor)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  const double aVal = 12;
  const double bVal = 5;
  QPoint a(aVal);
  QPoint b(bVal);

  QVec vec1(a, b);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec1[i], bVal - aVal);
  }

  QVec vec2(b, a);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec2[i], aVal - bVal);
  }

  QVec vec3(b);
  QVec vec4(QPoint::zero(), b);
  for(int i = 0; i < DIM; ++i)
  {
    EXPECT_EQ(vec3[i], bVal);
    EXPECT_EQ(vec4[i], bVal);
  }
  EXPECT_EQ(vec3, vec4);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_normalize)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  EXPECT_DOUBLE_EQ(QVec().unitVector().norm(), 1.0);

  // Find the norm of an arbitrary point
  QPoint p = QPoint::make_point(543.5435, 1566.423532, -432.4);
  QVec vec2(QPoint(), p);
  EXPECT_DOUBLE_EQ(vec2.unitVector().norm(), 1.0);

  // Find the norm of the zero vector
  // Zero vector should become the unit vector when normalized
  QVec vecZero;
  QVec unitVec(1, 1);
  EXPECT_DOUBLE_EQ(vecZero.unitVector().norm(), 1.0);
  EXPECT_EQ(vecZero.unitVector(), unitVec);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_norm)
{
  static const int DIM = 2;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  QPoint p1 = QPoint::make_point(3, 0);
  QPoint p2 = QPoint::make_point(0, 4);
  QVec vec(p1, p2);
  EXPECT_DOUBLE_EQ(vec.norm(), 5.0);
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_arithmetic)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  QPoint p1 = QPoint::make_point(3, 0, 1.2);
  QPoint p2 = QPoint::make_point(0, 4, 1.2);
  QPoint pSum = QPoint::make_point(3, 4, 2.4);
  QPoint pDiff = QPoint::make_point(-3, 4, 0);

  CoordType scalar = 5.3;
  QPoint pScalar = QPoint::make_point(scalar * 3, scalar * 4, scalar * 2.4);

  QVec v1(p1);
  QVec v2(p2);

  // testing vector addition and subtraction
  EXPECT_EQ(QVec(p1, p2), v2 - v1);
  EXPECT_EQ(QVec(pDiff), v2 - v1);
  EXPECT_EQ(QVec(pSum), v1 + v2);
  EXPECT_EQ(v1 + v2, v2 + v1);

  QVec vSum(v1);
  vSum += v2;
  EXPECT_EQ(vSum, QVec(pSum));

  QVec vDiff(v2);
  vDiff -= v1;
  EXPECT_EQ(vDiff, QVec(pDiff));

  // Testing scalar multiplication
  EXPECT_EQ(QVec(pScalar), QVec(pSum) * scalar);
  EXPECT_EQ(QVec(pScalar), scalar * QVec(pSum));
  EXPECT_EQ(QVec(pScalar), (v1 + v2) * scalar);

  // Testing unary negation
  EXPECT_EQ(-v1, QVec(p1, QPoint::zero()));
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_inner_product)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  QPoint p1 = QPoint::make_point(3, 0, 1.2);
  QPoint p2 = QPoint::make_point(0, 4, 1.2);

  CoordType expDot = 1.2 * 1.2;

  QVec v1(p1);
  QVec v2(p2);
  EXPECT_EQ(expDot, v1.dot(v2));
  EXPECT_EQ(v2.dot(v1), v1.dot(v2));
  EXPECT_EQ(v1.dot(v2), QVec::dot_product(v1, v2));
  EXPECT_EQ(QVec::dot_product(v1, v2), QVec::dot_product(v2, v1));

  QVec zer;
  EXPECT_EQ(0., zer.dot(zer));
}

//------------------------------------------------------------------------------
TEST(primal_vector, vector_outer_product)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef primal::Point<CoordType, DIM> QPoint;
  typedef primal::Vector<CoordType, DIM> QVec;

  QPoint p1 = QPoint::make_point(3, 0);
  QPoint p2 = QPoint::make_point(0, 4);

  QVec v1(p1);
  QVec v2(p2);

  QVec vRes;
  vRes[2] = 12.;

  EXPECT_EQ(QVec::cross_product(v1, v2), vRes);
  EXPECT_EQ(QVec::cross_product(v2, v1), -vRes);
  EXPECT_EQ(QVec::cross_product(v1, v2), -QVec::cross_product(v2, v1));
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
