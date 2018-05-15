/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file slam_map_DynamicMap.cpp
 *
 * \brief Unit tests for slam's DynamicMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "slam/Utilities.hpp"
#include "slam/DynamicSet.hpp"
#include "slam/DynamicMap.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;


typedef axom::slam::DynamicSet<> SetType;
typedef SetType::PositionType SetPosition;
typedef SetType::ElementType SetElement;

static const SetPosition MAX_SET_SIZE = 10;
static const SetPosition ADDITIONAL_ADD_SIZE = 5;

typedef axom::slam::DynamicMap<int>    IntMap;
typedef axom::slam::DynamicMap<double> RealMap;

TEST(slam_map,construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(AXOM_NULLPTR, m.set());
}

TEST(slam_map,construct_from_set_int)
{
  SetType s(MAX_SET_SIZE);

  IntMap m(&s);
  EXPECT_TRUE(m.isValid(true));

  EXPECT_EQ(MAX_SET_SIZE, m.size());
  EXPECT_EQ(s.size(), m.size());
}

TEST(slam_map,construct_from_set_real)
{
  SetType s(MAX_SET_SIZE);

  const double DEFAULT_VAL = -1.0;
  RealMap m(&s, DEFAULT_VAL);
  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(&s, m.set());

  EXPECT_EQ(MAX_SET_SIZE, m.size());
  EXPECT_EQ(s.size(), m.size());

  for(int i=0 ; i< MAX_SET_SIZE ; ++i)
  {
    EXPECT_EQ(DEFAULT_VAL, m[i]);
  }

  // Check that we can update the data
  const double UPDATE_VAL = 5.5;
  m[2] = UPDATE_VAL;
  EXPECT_EQ(UPDATE_VAL, m[2]);
}

// Check that resizing the map to a smaller size
// TODO: This shouldn't shrink the data
TEST(slam_map,modify_size)
{
  SetType s(MAX_SET_SIZE);

  IntMap m(&s);

  // Resizing to a larger size increases the map's size
  {
    const int newSize = 2 * MAX_SET_SIZE;
    m.resize(newSize);

    EXPECT_EQ(newSize, m.size());
  }

  // Resizing to a larger size increases the map's size
  {
    const int oldSize = m.size();
    const int newSize = MAX_SET_SIZE / 2;
    m.resize(newSize);

    EXPECT_GT(oldSize, newSize);
  }
}

// Check that resizing the map
// TODO: This should resize the underlying set
TEST(slam_map,modify_size_increases_set)
{
  SetType s(MAX_SET_SIZE);
  IntMap m(&s);

  const int newSize = 2 * MAX_SET_SIZE;
  m.resize(newSize);

  EXPECT_EQ(newSize, m.size());
}



//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  int result = RUN_ALL_TESTS();

  return result;
}
