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
  EXPECT_EQ(EMPTY_ID, CONDUIT_EMPTY_ID);
  EMPTY_T = getTypeID(0);
  EXPECT_EQ(EMPTY_ID, CONDUIT_EMPTY_ID);

  TypeID OBJECT_T = getTypeID<1>();
  EXPECT_EQ(OBJECT_ID, CONDUIT_OBJECT_ID);
  OBJECT_T = getTypeID(1);
  EXPECT_EQ(OBJECT_ID, CONDUIT_OBJECT_ID);

  TypeID LIST_T = getTypeID<2>();
  EXPECT_EQ(LIST_ID, CONDUIT_LIST_ID);
  LIST_T = getTypeID(2);
  EXPECT_EQ(LIST_ID, CONDUIT_LIST_ID);
#endif

  EXPECT_EQ(getTypeID<SIDRE_INT8_ID>(), CONDUIT_INT8_T);
  EXPECT_EQ(getTypeID(SIDRE_INT8_ID), CONDUIT_INT8_T);

  EXPECT_EQ(getTypeID<SIDRE_INT16_ID>(), CONDUIT_INT16_T);
  EXPECT_EQ(getTypeID(SIDRE_INT16_ID), CONDUIT_INT16_T);

  EXPECT_EQ(getTypeID<SIDRE_INT32_ID>(), CONDUIT_INT32_T);
  EXPECT_EQ(getTypeID(SIDRE_INT32_ID), CONDUIT_INT32_T);

  EXPECT_EQ(getTypeID<SIDRE_INT64_ID>(), CONDUIT_INT64_T);
  EXPECT_EQ(getTypeID(SIDRE_INT64_ID), CONDUIT_INT64_T);

  EXPECT_EQ(getTypeID<SIDRE_UINT8_ID>(), CONDUIT_UINT8_T);
  EXPECT_EQ(getTypeID(SIDRE_UINT8_ID), CONDUIT_UINT8_T);

  EXPECT_EQ(getTypeID<SIDRE_UINT16_ID>(), CONDUIT_UINT16_T);
  EXPECT_EQ(getTypeID(SIDRE_UINT16_ID), CONDUIT_UINT16_T);

  EXPECT_EQ(getTypeID<SIDRE_UINT32_ID>(), CONDUIT_UINT32_T);
  EXPECT_EQ(getTypeID(SIDRE_UINT32_ID), CONDUIT_UINT32_T);

  EXPECT_EQ(getTypeID<SIDRE_UINT64_ID>(), CONDUIT_UINT64_T);
  EXPECT_EQ(getTypeID(SIDRE_UINT64_ID), CONDUIT_UINT64_T);

  EXPECT_EQ(getTypeID<SIDRE_FLOAT32_ID>(), CONDUIT_FLOAT32_T);
  EXPECT_EQ(getTypeID(SIDRE_FLOAT32_ID), CONDUIT_FLOAT32_T);

  EXPECT_EQ(getTypeID<SIDRE_FLOAT64_ID>(), CONDUIT_FLOAT64_T);
  EXPECT_EQ(getTypeID(SIDRE_FLOAT64_ID), CONDUIT_FLOAT64_T);

  EXPECT_EQ(getTypeID<SIDRE_CHAR8_STR_ID>(), CONDUIT_CHAR8_STR_T);
  EXPECT_EQ(getTypeID(SIDRE_CHAR8_STR_ID), CONDUIT_CHAR8_STR_T);


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
