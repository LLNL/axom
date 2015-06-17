/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.h"
#include "sidre/SidreWrapperHelpers.hpp"

using asctoolkit::sidre::DataType;
using asctoolkit::sidre::TypeID;
using asctoolkit::sidre::getTypeID;
//------------------------------------------------------------------------------
// This test verifies that the getTypeID C version and C++ (templated) version
// returns the same types. 
//------------------------------------------------------------------------------

TEST(sidre_types,get_sidre_type)
{
#if 0
  TypeID EMPTY_T = getTypeID<0>();
  EXPECT_EQ(EMPTY_T, DataType::EMPTY_T);
  EMPTY_T = getTypeID(0);
  EXPECT_EQ(EMPTY_T, DataType::EMPTY_T);

  TypeID OBJECT_T = getTypeID<1>();
  EXPECT_EQ(OBJECT_T, DataType::OBJECT_T);
  OBJECT_T = getTypeID(1);
  EXPECT_EQ(OBJECT_T, DataType::OBJECT_T);

  TypeID LIST_T = getTypeID<2>();
  EXPECT_EQ(LIST_T, DataType::LIST_T);
  LIST_T = getTypeID(2);
  EXPECT_EQ(LIST_T, DataType::LIST_T);
#endif

  TypeID INT8_T = getTypeID<ATK_INT8_T>();
  EXPECT_EQ(INT8_T, DataType::INT8_T);
  INT8_T = getTypeID(ATK_INT8_T);
  EXPECT_EQ(INT8_T, DataType::INT8_T);

  TypeID INT16_T = getTypeID<ATK_INT16_T>();
  EXPECT_EQ(INT16_T, DataType::INT16_T);
  INT16_T = getTypeID(ATK_INT16_T);
  EXPECT_EQ(INT16_T, DataType::INT16_T);

  TypeID INT32_T = getTypeID<ATK_INT32_T>();
  EXPECT_EQ(INT32_T, DataType::INT32_T);
  INT32_T = getTypeID(ATK_INT32_T);
  EXPECT_EQ(INT32_T, DataType::INT32_T);

  TypeID INT64_T = getTypeID<ATK_INT64_T>();
  EXPECT_EQ(INT64_T, DataType::INT64_T);
  INT64_T = getTypeID(ATK_INT64_T);
  EXPECT_EQ(INT64_T, DataType::INT64_T);

  TypeID UINT8_T = getTypeID<ATK_UINT8_T>();
  EXPECT_EQ(UINT8_T, DataType::UINT8_T);
  UINT8_T = getTypeID(ATK_UINT8_T);
  EXPECT_EQ(UINT8_T, DataType::UINT8_T);

  TypeID UINT16_T = getTypeID<ATK_UINT16_T>();
  EXPECT_EQ(UINT16_T, DataType::UINT16_T);
  UINT16_T = getTypeID(ATK_UINT16_T);
  EXPECT_EQ(UINT16_T, DataType::UINT16_T);

  TypeID UINT32_T = getTypeID<ATK_UINT32_T>();
  EXPECT_EQ(UINT32_T, DataType::UINT32_T);
  UINT32_T = getTypeID(ATK_UINT32_T);
  EXPECT_EQ(UINT32_T, DataType::UINT32_T);

  TypeID UINT64_T = getTypeID<ATK_UINT64_T>();
  EXPECT_EQ(UINT64_T, DataType::UINT64_T);
  UINT64_T = getTypeID(ATK_UINT64_T);
  EXPECT_EQ(UINT64_T, DataType::UINT64_T);

  TypeID FLOAT32_T = getTypeID<ATK_FLOAT32_T>();
  EXPECT_EQ(FLOAT32_T, DataType::FLOAT32_T);
  FLOAT32_T = getTypeID(ATK_FLOAT32_T);
  EXPECT_EQ(FLOAT32_T, DataType::FLOAT32_T);

  TypeID FLOAT64_T = getTypeID<ATK_FLOAT64_T>();
  EXPECT_EQ(FLOAT64_T, DataType::FLOAT64_T);
  FLOAT64_T = getTypeID(ATK_FLOAT64_T);
  EXPECT_EQ(FLOAT64_T, DataType::FLOAT64_T);

  TypeID CHAR8_STR_T = getTypeID<ATK_CHAR8_STR_T>();
  EXPECT_EQ(CHAR8_STR_T, DataType::CHAR8_STR_T);
  CHAR8_STR_T = getTypeID(ATK_CHAR8_STR_T);
  EXPECT_EQ(CHAR8_STR_T, DataType::CHAR8_STR_T);


}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
