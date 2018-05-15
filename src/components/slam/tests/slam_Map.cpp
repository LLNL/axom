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
 * \file slam_Map.cpp
 *
 * \brief Unit tests for Slam's Map
 */

#include <iterator>
#include "gtest/gtest.h"

#include "slam/Utilities.hpp"
#include "slam/RangeSet.hpp"
#include "slam/Map.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;


typedef axom::slam::RangeSet SetType;
typedef axom::slam::Map<int>    IntMap;
typedef axom::slam::Map<double> RealMap;

typedef SetType::PositionType PositionType;
typedef SetType::ElementType ElementType;

static PositionType const MAX_SET_SIZE = 10;

TEST(slam_map,construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
}

template<typename T>
bool constructAndTestMap()
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("\nCreating set of size " << s.size() );

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO(
    "\nCreating "
    << axom::slam::util::TypeToString<T>::to_string() << " map on the set ");
  axom::slam::Map<T> m(&s);
  EXPECT_TRUE(m.isValid());

  SLIC_INFO( "\nSetting the elements.");
  double multFac = 1.0001;
  for(PositionType idx = 0 ; idx < m.size() ; ++idx)
  {
    m[idx] = static_cast<T>(idx * multFac);
  }

  SLIC_INFO("\nChecking the elements.");
  for(PositionType idx = 0 ; idx < m.size() ; ++idx)
  {
    EXPECT_EQ(m[idx], static_cast<T>(idx * multFac) );
  }

  EXPECT_TRUE(m.isValid(true));

  return true;
}

TEST(slam_map,construct_int_map)
{
  EXPECT_TRUE( constructAndTestMap<int>() );
}

TEST(slam_map,construct_double_map)
{
  EXPECT_TRUE( constructAndTestMap<double>());
}

TEST(slam_map,out_of_bounds)
{
  int defaultElt = 2;

  SetType s(MAX_SET_SIZE);
  IntMap m(&s, defaultElt);

  SLIC_INFO("Testing Map element access -- in bounds");
  for(PositionType idx = 0 ; idx < m.size() ; ++idx)
    EXPECT_EQ(defaultElt, m[idx]);

  // Test out of bounds
  SLIC_INFO("Testing Map element access "
            << "-- out of bounds access; Expecting the test to fail");
  #ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(  m[-1],      "")
    << " Accessed element -1 of Map -- out of bounds";
  EXPECT_DEATH_IF_SUPPORTED(  m[m.size()],"")
    << " Accessed element " << m.size() << " of Map -- out of bounds";

  #else
  SLIC_INFO("Skipped assertion failure check in release mode.");
  #endif
}


//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  int result = RUN_ALL_TESTS();

  return result;
}
