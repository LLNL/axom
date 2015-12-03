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

/* Do we need any of this type testing?  The sidre types are straight up defined to conduit types now.
  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_INT8_ID>()), 
            static_cast<int>(CONDUIT_INT8_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_INT8_ID)), 
            static_cast<int>(CONDUIT_INT8_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_INT16_ID>()), 
            static_cast<int>(CONDUIT_INT16_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_INT16_ID)), 
            static_cast<int>(CONDUIT_INT16_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_INT32_ID>()), 
            static_cast<int>(CONDUIT_INT32_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_INT32_ID)), 
            static_cast<int>(CONDUIT_INT32_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_INT64_ID>()), 
            static_cast<int>(CONDUIT_INT64_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_INT64_ID)), 
            static_cast<int>(CONDUIT_INT64_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_UINT8_ID>()), 
            static_cast<int>(CONDUIT_UINT8_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_UINT8_ID)), 
            static_cast<int>(CONDUIT_UINT8_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_UINT16_ID>()), 
            static_cast<int>(CONDUIT_UINT16_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_UINT16_ID)), 
            static_cast<int>(CONDUIT_UINT16_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_UINT32_ID>()), 
            static_cast<int>(CONDUIT_UINT32_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_UINT32_ID)), 
            static_cast<int>(CONDUIT_UINT32_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_UINT64_ID>()), 
            static_cast<int>(CONDUIT_UINT64_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_UINT64_ID)), 
            static_cast<int>(CONDUIT_UINT64_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_FLOAT32_ID>()), 
            static_cast<int>(CONDUIT_FLOAT32_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_FLOAT32_ID)), 
            static_cast<int>(CONDUIT_FLOAT32_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_FLOAT64_ID>()), 
            static_cast<int>(CONDUIT_FLOAT64_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_FLOAT64_ID)), 
            static_cast<int>(CONDUIT_FLOAT64_T));

  EXPECT_EQ(static_cast<int>(getTypeID<SIDRE_CHAR8_STR_ID>()), 
            static_cast<int>(CONDUIT_CHAR8_STR_T));
  EXPECT_EQ(static_cast<int>(getTypeID(SIDRE_CHAR8_STR_ID)), 
            static_cast<int>(CONDUIT_CHAR8_STR_T));

*/
}
