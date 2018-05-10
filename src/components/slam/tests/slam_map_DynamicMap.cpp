/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
 * \brief Unit tests for DynamicMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "slam/Utilities.hpp"
#include "slam/DynamicMap.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;


//typedef axom::slam::RangeSet    SetType;
typedef axom::slam::DynamicMap<int>    IntMap;
typedef axom::slam::DynamicMap<double> RealMap;

//typedef SetType::PositionType   PositionType;
//typedef SetType::ElementType    ElementType;

//static PositionType const MAX_SET_SIZE = 10;

TEST(gtest_slam_map,construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
}



//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
