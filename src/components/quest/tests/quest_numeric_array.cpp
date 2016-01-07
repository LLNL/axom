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

#include "quest/NumericArray.hpp"

//------------------------------------------------------------------------------
TEST( quest_numeric_array, constructors)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef quest::NumericArray<CoordType, DIM> QArray;

  QArray arr1;
  EXPECT_EQ(arr1.dimension(), DIM );
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(arr1[i], CoordType() );

  CoordType val = 5.;
  QArray arr2(val);
  EXPECT_EQ(arr2.dimension(), DIM );
  for(int i=0; i< DIM; ++i)
      EXPECT_EQ(arr2[i], val );


  CoordType val3 = 5.;
  int halfPt = DIM /2;
  QArray arr3(val3, halfPt);
  for(int i=0; i< DIM; ++i)
  {
      CoordType expVal=  i < halfPt ? val3 : CoordType();
      EXPECT_EQ(arr3[i], expVal );
  }

  // Compare array initialized from arbitrary array
  CoordType valsArr[DIM] = {12.,23., 34., 45., 56.432};
  QArray arrArr(valsArr);
  for(int i=0; i< DIM; ++i)
  {
      EXPECT_EQ(arrArr[i], valsArr[i] );
  }

  QArray arrHalfVec(valsArr, halfPt);
  for(int i=0; i< DIM; ++i)
  {
      CoordType expVal=  i < halfPt ? valsArr[i] : CoordType();
      EXPECT_EQ(arrHalfVec[i], expVal );
  }

}

//------------------------------------------------------------------------------
TEST( quest_numeric_array, num_array_to_array)
{
  static const int DIM = 5;
  typedef double CoordType;
  typedef quest::NumericArray<CoordType, DIM> QArray;

  // Compare array initialized from arbitrary array
  CoordType valsArr[DIM] = {12.,23., 34., 45., 56.432};
  QArray arrArr(valsArr);

  for(int i=0; i< DIM; ++i)
  {
      EXPECT_EQ(arrArr[i], valsArr[i] );
  }

  // Copy data into new array
  CoordType arrCopy[DIM];
  arrArr.to_array(arrCopy);
  for(int i=0; i< DIM; ++i)
  {
      EXPECT_EQ(arrCopy[i], valsArr[i] );
      EXPECT_EQ(arrCopy[i], arrArr.data()[i] );
  }



}

//------------------------------------------------------------------------------
TEST( quest_numeric_array, component_wise_arithmetic)
{
  static const int DIM = 3;
  typedef double CoordType;
  typedef quest::NumericArray<CoordType, DIM> QArray;

  CoordType ca1[]    = { 3, 0, 1.2};
  CoordType ca2[]    = { 0, 4, 1.2};
  CoordType caSum[]  = { 3, 4, 2.4};
  CoordType caDiff[] = {-3, 4, 0};
  CoordType scalar = 5.3;

  CoordType caScalar[] = {scalar * 3, scalar * 4, scalar * 2.4};

  CoordType caProduct[] = {-3 * 3, 4 * 4, 0. * 2.4};

  QArray arr1(ca1);
  QArray arr2(ca2);

  // testing component-wise addition and subtraction
  EXPECT_EQ( QArray(caDiff), arr2-arr1);
  EXPECT_EQ( QArray(caSum), arr2+arr1);
  EXPECT_EQ( arr2+arr1, arr1+arr2);

  QArray aSum(arr1);
  aSum += arr2;
  EXPECT_EQ(aSum, QArray(caSum));

  QArray aDiff(arr2);
  aDiff -= arr1;
  EXPECT_EQ(aDiff, QArray(caDiff));

  // Testing scalar multiplication
  EXPECT_EQ( QArray(caScalar), QArray(caSum) * scalar  );
  EXPECT_EQ( QArray(caScalar), scalar * QArray(caSum) );
  EXPECT_EQ( QArray(caScalar), (arr1+arr2) * scalar);

  // Testing unary negation
  EXPECT_EQ( -arr1, QArray() - arr1);

  // Testing component-wise product
  EXPECT_EQ( QArray(caProduct), aSum * aDiff);
  EXPECT_EQ( QArray(caProduct), aDiff * aSum);

  // Testing component-wise product
  EXPECT_EQ( QArray(caProduct) / aSum, aDiff);

}



//------------------------------------------------------------------------------
TEST( quest_numeric_array, component_min_max)
{
  static const int DIM = 3;
  typedef int CoordType;
  typedef quest::NumericArray<CoordType, DIM> QArray;

  CoordType incSeq[]    = { 1,2,3};
  CoordType decSeq[]    = { 3,2,1};
  CoordType upDownSeq[] = { 5,8,3};
  CoordType downUpSeq[] = {6,2,5};

  QArray incArr(incSeq);
  QArray decArr(decSeq);
  QArray udArr(upDownSeq);
  QArray duArr(downUpSeq);


  // testing component-wise min and max functions
  EXPECT_EQ( incArr.max(), 3);
  EXPECT_EQ( decArr.max(), 3);
  EXPECT_EQ( incArr.argMax(), 2);
  EXPECT_EQ( decArr.argMax(), 0);

  EXPECT_EQ( udArr.max(), 8);
  EXPECT_EQ( udArr.argMax(), 1);


  EXPECT_EQ( incArr.min(), 1);
  EXPECT_EQ( decArr.min(), 1);
  EXPECT_EQ( incArr.argMin(), 0);
  EXPECT_EQ( decArr.argMin(), 2);

  EXPECT_EQ( duArr.min(), 2);
  EXPECT_EQ( duArr.argMin(), 1);


}


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

