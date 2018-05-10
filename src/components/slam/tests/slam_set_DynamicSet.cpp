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

/*
 * \file slam_set_DynamicSet.cpp
 *
 * \brief ...
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"    // for AXOM_USE_BOOST

#include "slic/slic.hpp"        // for SLIC_INFO

#include "slam/DynamicSet.hpp"
#include "slam/Utilities.hpp"


typedef axom::slam::DynamicSet<> SetType;
typedef SetType::PositionType SetPosition;
typedef SetType::ElementType SetElement;

static const SetPosition MAX_SET_SIZE = 10;
static const SetPosition ADDITIONAL_ADD_SIZE = 5;

TEST(gtest_slam_set_dynamicset,construct)
{
  SLIC_INFO("Testing constructing a DynamicSet");

  //Empty set
  SetType s0;
  EXPECT_TRUE(s0.isValid() );
  EXPECT_TRUE(s0.empty() );
  EXPECT_EQ(  0,  s0.size() );

  //Set of size 0
  SetType s1(0);
  EXPECT_TRUE(s1.isValid());
  EXPECT_EQ(0, s1.size());
  EXPECT_TRUE(s1.empty());

  EXPECT_EQ(s1,s0);

  //Construct with a size
  SetType s2(MAX_SET_SIZE);
  EXPECT_TRUE(s2.isValid() );
  EXPECT_FALSE(s2.empty() );
  EXPECT_EQ(MAX_SET_SIZE, s2.size() );

}



TEST(gtest_slam_set_dynamicset,construct_strictly_positive)
{
  SLIC_INFO("Testing that DynamicSet with negative size are invalid");

  const int NEGATIVE_NUMBER = -2;
  SetType s(NEGATIVE_NUMBER);
  EXPECT_FALSE( s.isValid());
}


TEST(gtest_slam_set_dynamicset,construct_set_builder)
{
  SLIC_INFO("Testing construction of DynamicSet using SetBuilders");

  typedef SetType::SetBuilder SetBuilder;

  SetBuilder builder = SetBuilder()
                       .size( MAX_SET_SIZE );

  SetType s(builder);
  EXPECT_TRUE( s.isValid());
  EXPECT_EQ( MAX_SET_SIZE, s.size());

  const int ZERO_OFFSET = 0;
  const int DEFAULT_STRIDE = 1;
  SetBuilder ok_builder = SetBuilder()
                          .size(MAX_SET_SIZE)
                          .offset(ZERO_OFFSET)
                          .stride(DEFAULT_STRIDE);
  SetType s2(ok_builder);
  EXPECT_TRUE( s2.isValid() );
  EXPECT_EQ(  MAX_SET_SIZE, s2.size());

  // The two sets should be equal
  EXPECT_EQ( s, s2);

#ifdef AXOM_DEBUG
  // Using inappropriate SetBuilder features generates an assert failure
  SLIC_INFO("Cannot construct a DynamicSet with an invalid SetBuilder");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  const int NON_ZERO_OFFSET = 3;
  EXPECT_DEATH_IF_SUPPORTED( SetType(SetBuilder()
                                     .size(MAX_SET_SIZE)
                                     .offset(NON_ZERO_OFFSET)), "" );
#endif
}


TEST(gtest_slam_set_dynamicset,construct_runtime)
{
  SLIC_INFO("Testing that DynamicSet can be constructed with runtime sizes");

  for(int i = 10 ; i < 15 ; ++i)
  {
    // mark as volatile to ensure compiler doesn't see this at compile time
    volatile int v = i;
    SetType s(v);

    EXPECT_TRUE( s.isValid());
    EXPECT_EQ( i, s.size());
  }
}




TEST(gtest_slam_set_dynamicset,out_of_bounds_at)
{
  SLIC_INFO(
    "Testing out of bounds access using at() "
    <<"-- code is expected to assert and die.");
  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s.at(MAX_SET_SIZE),"");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}

TEST(gtest_slam_set_dynamicset,out_of_bounds_bracket)
{
  SLIC_INFO(
    "Testing out of bounds access using operator[] "
    << "-- code is expected to assert and die.");

  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}



TEST(gtest_slam_set_dynamicset,adding_elements)
{
  SLIC_INFO("Testing adding new elements in the set");

  SetType s(MAX_SET_SIZE);

  for(SetPosition i = 0 ; i < ADDITIONAL_ADD_SIZE ; i++)
  {
    s.insert();
  }

  EXPECT_EQ(s.size(), MAX_SET_SIZE + ADDITIONAL_ADD_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("done.");
}


TEST(gtest_slam_set_dynamicset,removing_elements)
{
  SLIC_INFO("Testing removing elements in the set");

  SetType s(MAX_SET_SIZE);

  for(SetPosition i = 0 ; i < MAX_SET_SIZE ; i+=2) //remove every other entry
  {
    s.remove(i);
  }

  EXPECT_EQ(MAX_SET_SIZE, s.size());
  EXPECT_EQ(MAX_SET_SIZE / 2, s.numberOfValidEntries() );
  EXPECT_TRUE(s.isValid());

  for(SetPosition i = 1 ; i < MAX_SET_SIZE ; i+=2)
  {
    EXPECT_EQ( s[i], i );
    EXPECT_TRUE( s.isValidEntry(i) );
  }

  for(SetPosition i = 0 ; i < MAX_SET_SIZE ; i+=2)
  {
    EXPECT_EQ( s[i], (int)SetType::INVALID_ENTRY );
    EXPECT_FALSE( s.isValidEntry(i) );
  }

  SLIC_INFO("done.");
}


TEST(gtest_slam_set_dynamicset, copy_assignment)
{
  SLIC_INFO("Testing the assignment operator");

  SetType s;

  {
    SetType s_cp(MAX_SET_SIZE);

    for(SetPosition i = 0 ; i < MAX_SET_SIZE ; i+=2) //remove every other entry
    {
      s_cp.remove(i);
    }

    s = s_cp;
    EXPECT_TRUE(s.isValid());
  }

  EXPECT_EQ(s.numberOfValidEntries(), MAX_SET_SIZE / 2);
  EXPECT_EQ(s.size(), MAX_SET_SIZE );
  EXPECT_TRUE(s.isValid());

  for(SetPosition i = 1 ; i < MAX_SET_SIZE ; i+=2)
  {
    EXPECT_EQ( s[i], i );
    EXPECT_TRUE( s.isValidEntry(i) );
  }

  for(SetPosition i = 0 ; i < MAX_SET_SIZE ; i+=2)
  {
    EXPECT_EQ( s[i], (int)SetType::INVALID_ENTRY );
    EXPECT_FALSE( s.isValidEntry(i) );
  }

  SLIC_INFO("done.");
}




TEST(gtest_slam_set_dynamicset,out_of_bounds_remove)
{
  SLIC_INFO(
    "Testing out of bounds access using remove() "
    << "-- code is expected to assert and die.");

  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s.remove(MAX_SET_SIZE), "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}



TEST(gtest_slam_set_dynamicset,find_index)
{
  SLIC_INFO("Testing finding index of an element");

  SetType s(MAX_SET_SIZE);

  for(SetPosition i = 0 ; i < MAX_SET_SIZE ; i++)
  {
    EXPECT_EQ( s.findIndex(i), i );
  }

  SLIC_INFO("done.");
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  UnitTestLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
