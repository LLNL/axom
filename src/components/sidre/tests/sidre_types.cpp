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

//#include "conduit/conduit.h"

using asctoolkit::sidre::DataType;
using asctoolkit::sidre::TypeID;
using asctoolkit::sidre::getTypeID;
//------------------------------------------------------------------------------

TEST(sidre_types,get_sidre_type)
{
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

  TypeID INT8_T = getTypeID<3>();
  EXPECT_EQ(INT8_T, DataType::INT8_T);
  INT8_T = getTypeID(3);
  EXPECT_EQ(INT8_T, DataType::INT8_T);

  TypeID INT16_T = getTypeID<4>();
  EXPECT_EQ(INT16_T, DataType::INT16_T);
  INT16_T = getTypeID(4);
  EXPECT_EQ(INT16_T, DataType::INT16_T);

  TypeID INT32_T = getTypeID<5>();
  EXPECT_EQ(INT32_T, DataType::INT32_T);
  INT32_T = getTypeID(5);
  EXPECT_EQ(INT32_T, DataType::INT32_T);

  TypeID INT64_T = getTypeID<6>();
  EXPECT_EQ(INT64_T, DataType::INT64_T);
  INT64_T = getTypeID(6);
  EXPECT_EQ(INT64_T, DataType::INT64_T);

  TypeID UINT8_T = getTypeID<7>();
  EXPECT_EQ(UINT8_T, DataType::UINT8_T);
//  UINT8_T = getTypeID(7);
  EXPECT_EQ(UINT8_T, DataType::UINT8_T);

  TypeID UINT16_T = getTypeID<8>();
  EXPECT_EQ(UINT16_T, DataType::UINT16_T);
//  UINT16_T = getTypeID(8);
  EXPECT_EQ(UINT16_T, DataType::UINT16_T);

  TypeID UINT32_T = getTypeID<9>();
  EXPECT_EQ(UINT32_T, DataType::UINT32_T);
//  UINT32_T = getTypeID(9);
  EXPECT_EQ(UINT32_T, DataType::UINT32_T);

  TypeID UINT64_T = getTypeID<10>();
  EXPECT_EQ(UINT64_T, DataType::UINT64_T);
//  UINT64_T = getTypeID(10);
  EXPECT_EQ(UINT64_T, DataType::UINT64_T);

  TypeID FLOAT32_T = getTypeID<11>();
  EXPECT_EQ(FLOAT32_T, DataType::FLOAT32_T);
//  FLOAT32_T = getTypeID(11);
  EXPECT_EQ(FLOAT32_T, DataType::FLOAT32_T);

  TypeID FLOAT64_T = getTypeID<12>();
  EXPECT_EQ(FLOAT64_T, DataType::FLOAT64_T);
//  FLOAT64_T = getTypeID(12);
  EXPECT_EQ(FLOAT64_T, DataType::FLOAT64_T);

  TypeID CHAR8_STR_T = getTypeID<13>();
  EXPECT_EQ(CHAR8_STR_T, DataType::CHAR8_STR_T);
//  CHAR8_STR_T = getTypeID(13);
  EXPECT_EQ(CHAR8_STR_T, DataType::CHAR8_STR_T);


}

