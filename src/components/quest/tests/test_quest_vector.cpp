/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "quest/Vector.hpp"

//------------------------------------------------------------------------------
TEST( quest_vector, vector_constructors)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::Vector<CoordType, DIM> QVec;

  QVec vec1;
  EXPECT_EQ(vec1.dimension(), DIM );
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(vec1[i], CoordType() );

  CoordType val = 5.;
  QVec vec2(val);
  EXPECT_EQ(vec2.dimension(), DIM );
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(vec2[i], val );


  CoordType val3 = 5.;
  int halfPt = DIM /2;
  QVec vec3(val3, halfPt);
  for(int i=0; i< DIM; ++i)
  {
      CoordType expVal=  i < halfPt ? val3 : CoordType();
      EXPECT_EQ(vec3[i], expVal );
  }

  // Compare vector initialized from arbitrary
  CoordType valsArr[DIM] = {12.,23., 34., 45., 56.432};
  QVec arrVec(valsArr);
  for(int i=0; i< DIM; ++i)
  {
      EXPECT_EQ(arrVec[i], valsArr[i] );
  }

  QVec arrHalfVec(valsArr, halfPt);
  for(int i=0; i< DIM; ++i)
  {
      CoordType expVal=  i < halfPt ? valsArr[i] : CoordType();
      EXPECT_EQ(arrHalfVec[i], expVal );
  }

}

//------------------------------------------------------------------------------
TEST( quest_vector, vector_from_points_constructor)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::Vector<CoordType, DIM> QVec;

  const double aVal = 12;
  const double bVal = 5;
  QPoint a ( aVal );
  QPoint b ( bVal );

  QVec vec1(a,b);
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(vec1[i], bVal - aVal );

  QVec vec2(b,a);
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(vec2[i], aVal - bVal );


  QVec vec3(b);
  QVec vec4(QPoint::zero(), b);
  for(int i=0; i< DIM; ++i)
  {
      EXPECT_EQ(vec3[i], bVal );
      EXPECT_EQ(vec4[i], bVal );
  }
  EXPECT_EQ( vec3, vec4);

}


//------------------------------------------------------------------------------
TEST( quest_vector, vector_normalize)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::Vector<CoordType, DIM> QVec;

  QVec vec;
  vec.normalize();
  EXPECT_DOUBLE_EQ(vec.norm(), 1.0 );

  // Find the norm of an arbitrary point
  QPoint p = QPoint::make_point(543.5435, 1566.423532, -432.4);
  QVec vec2(QPoint(), p);
  vec2.normalize();
  EXPECT_DOUBLE_EQ(vec2.norm(), 1.0 );


  // Find the norm of the zero vector
  // Zero vector should become the unit vector when normalized
  QVec vecZero;
  QVec unitVec(1,1);
  vecZero.normalize();
  EXPECT_DOUBLE_EQ(vecZero.norm(), 1.0 );
  EXPECT_EQ(vecZero, unitVec );
}

//------------------------------------------------------------------------------
TEST( quest_vector, vector_norm)
{
  static const int DIM = 2;
  typedef double CoordType;
  typedef quest::Point<CoordType, DIM> QPoint;
  typedef quest::Vector<CoordType, DIM> QVec;


  QPoint p1 = QPoint::make_point(3,0);
  QPoint p2 = QPoint::make_point(0,4);
  QVec vec(p1,p2);
  EXPECT_DOUBLE_EQ(vec.norm(), 5.0 );
}

///// Still to test
///     arithmetic operators, outer and inner product

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

