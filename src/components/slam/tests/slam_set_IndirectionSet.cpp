/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file slam_set_indirectionset.cpp
 */


#include <iterator>
#include <sstream>                  // for std::stringstream
#include <algorithm>                // for std::next_permutation
#include <iterator>                 // for std::ostream_iterator
#include "gtest/gtest.h"

#include "axom/config.hpp"        // for AXOM_USE_BOOST
#include "axom/Types.hpp"   // for AXOM_NULLPTR

#include "slic/slic.hpp"

#include "slam/Utilities.hpp"
#include "slam/IndirectionSet.hpp"
#include "slam/RangeSet.hpp"        // for PositionSet


namespace {
  typedef axom::slam::ArrayIndirectionSet SetType;
  typedef SetType::PositionType           SetPosition;
  typedef SetType::ElementType            SetElement;

  static const SetPosition MAX_SET_SIZE = 10;

  SetPosition* allocIncrementingArray(int size)
  {
    SLIC_ASSERT(size > 0);

    SetPosition* arr = new SetPosition[size];

    for(int i = 0; i< size; ++i)
      arr[i] = i;

    return arr;
  }

  void permuteArray(SetPosition* arr, int size)
  {
    for(int i = 0; i<1000; ++i)
    {
      std::next_permutation(arr, arr + size);
    }
  }
}

TEST(gtest_slam_set_indirectionset,construct)
{
  SLIC_INFO("Testing constructors for IndirectionSet");

  SetType s_empty;
  EXPECT_TRUE(s_empty.empty());
  EXPECT_EQ(0, s_empty.size() );
  EXPECT_TRUE(s_empty.isValid());


  SetType s(MAX_SET_SIZE);
  EXPECT_FALSE(s.empty());
  EXPECT_EQ(MAX_SET_SIZE, s.size() );
  EXPECT_FALSE(s.isValid()) << "IndirectionSet not valid until we set the indirection buffer";

  s.data() = allocIncrementingArray(MAX_SET_SIZE);

  EXPECT_FALSE(s.empty());
  EXPECT_EQ(MAX_SET_SIZE, s.size() );
  EXPECT_TRUE(s.isValid());

  delete [] s.data();
}

TEST(gtest_slam_set_indirectionset,set_builder)
{
  SLIC_INFO("Testing construction of IndirectionSet using SetBuilders");

  typedef SetType::SetBuilder SetBuilder;

  bool bVerbose = true;

  {
    SetType s_empty(SetBuilder()
        .size(MAX_SET_SIZE) );

    EXPECT_FALSE(s_empty.empty());
    EXPECT_EQ(MAX_SET_SIZE, s_empty.size() );
    EXPECT_FALSE(s_empty.isValid(bVerbose));
  }

  {
    SetType s_empty2(SetBuilder()
        .data(allocIncrementingArray(MAX_SET_SIZE) ) );

    EXPECT_TRUE(s_empty2.empty());
    EXPECT_EQ(0, s_empty2.size() );
    EXPECT_TRUE(s_empty2.isValid(bVerbose));

    delete [] s_empty2.data();
  }

  {
    SetType s(SetBuilder()
        .size(MAX_SET_SIZE)
        .data(allocIncrementingArray(MAX_SET_SIZE) ) );

    EXPECT_FALSE(s.empty());
    EXPECT_EQ(MAX_SET_SIZE, s.size() );
    EXPECT_TRUE(s.isValid(bVerbose));

    delete [] s.data();
  }
}

TEST(gtest_slam_set_indirectionset,iterate)
{
  SetType s(SetType::SetBuilder()
      .size(MAX_SET_SIZE)
      .data(allocIncrementingArray(MAX_SET_SIZE) ) );

  SLIC_INFO("Iterating through set of size " << s.size());
  EXPECT_EQ(MAX_SET_SIZE, s.size());
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Using random access -- operator[] and at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement exp = static_cast<SetElement>(pos);
      EXPECT_EQ(exp,s[pos]);
      EXPECT_EQ(exp,s.at(pos));

      sstr << s[pos] << "\t";
    }
    SLIC_INFO("Data using operator[]:\t" << sstr.str());
  }

#ifdef AXOM_USE_BOOST
  SLIC_INFO("Using iterator interface");
  {
    std::stringstream sstr;
    typedef SetType::iterator SetIterator;
    for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      EXPECT_EQ( std::distance(s.begin(), it), *it );

      sstr << *it << "\t";
    }
    SLIC_INFO("Data from iterator:\t" << sstr.str());
  }

  const int NUM_PERMUTATIONS = 3;
  SLIC_INFO("Printing a few permutations of the data.");
  for(int i = 0; i< NUM_PERMUTATIONS; ++i)
  {
    permuteArray(s.data(), s.size());
    std::stringstream sstr;
    std::copy(s.begin(), s.end(),
        std::ostream_iterator<SetPosition>(sstr, "\t") );
    SLIC_INFO("Permutation " << i << ":\t" << sstr.str());
  }
#endif

  delete [] s.data();
}

TEST(gtest_slam_set_indirectionset,equality)
{
  SLIC_INFO("Testing equality/inequality for several sets");

  SetType s1(SetType::SetBuilder()
      .size(MAX_SET_SIZE)
      .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s1.isValid());

  // Test self-equality
  EXPECT_EQ(s1, s1);

  // Test equality against a set with the same array
  SetType s1b(SetType::SetBuilder()
      .size(MAX_SET_SIZE)
      .data(s1.data() ) );
  EXPECT_TRUE(s1b.isValid());
  EXPECT_EQ(s1, s1b);

  // Test equality against a set with a different array (same values)
  SetType s2(SetType::SetBuilder()
      .size(MAX_SET_SIZE)
      .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s2.isValid());
  EXPECT_EQ(s1, s2);
  EXPECT_EQ(s2, s1);

  // Check that sets are not equal if we permute the data
  permuteArray(s2.data(), s2.size());
  EXPECT_NE(s1, s2);

  // Test against a smaller set
  const int S3_SIZE = MAX_SET_SIZE / 2;
  SetType s3(SetType::SetBuilder()
      .size(S3_SIZE)
      .data(allocIncrementingArray(S3_SIZE) ) );
  EXPECT_TRUE(s3.isValid());
  EXPECT_NE(s1, s3);
  EXPECT_NE(s2, s3);


  // Check against a PositionSet with the same elements
  axom::slam::PositionSet posSet(MAX_SET_SIZE);
  EXPECT_TRUE(posSet.isValid() );
  EXPECT_EQ(posSet, s1);


  // Reclaim memory
  delete [] s1.data();
  delete [] s2.data();
  delete [] s3.data();
}


TEST(gtest_slam_set_indirectionset,out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized set-- code is expected to assert and die.");

  SetType s(SetType::SetBuilder()
      .size(MAX_SET_SIZE)
      .data(allocIncrementingArray(MAX_SET_SIZE) ) );
  EXPECT_TRUE(s.isValid());

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  delete [] s.data();
}


TEST(gtest_slam_set_indirectionset,vector_indirection)
{
  SLIC_INFO("Testing basic set operations on a VectorIndirectionSet");

  typedef axom::slam::VectorIndirectionSet  VectorIndirectionSet;
  typedef VectorIndirectionSet::SetBuilder  Builder;

  // Set up data -- an array of incrementing integers
  std::vector<int> intVec;
  intVec.reserve(MAX_SET_SIZE);
  for(int i = 0; i<MAX_SET_SIZE; ++i)
  {
    intVec.push_back(i);
  }

  // Setup the VectorIndirectionSet and test basic functionality
  VectorIndirectionSet vecSet(Builder()
      .size(MAX_SET_SIZE)
      .data(&intVec));
  EXPECT_TRUE(vecSet.isValid());
  EXPECT_FALSE(vecSet.empty());
  EXPECT_EQ(MAX_SET_SIZE, vecSet.size());
  EXPECT_TRUE(vecSet.hasIndirection());

  // Compare the values
  for(int i = 0; i<MAX_SET_SIZE; ++i)
  {
    EXPECT_EQ(  i,  vecSet[i]);
    EXPECT_EQ(  i,  vecSet.at(i));
  }

  // Test equality with an array indirection set
  {
    SetType arrSet(SetType::SetBuilder()
        .size(MAX_SET_SIZE)
        .data(allocIncrementingArray(MAX_SET_SIZE) ) );
    EXPECT_TRUE(arrSet.isValid());

    EXPECT_EQ(vecSet, arrSet);
    EXPECT_EQ(arrSet, vecSet);

    delete [] arrSet.data();
  }

  // Test equality with a position set
  {
    axom::slam::PositionSet posSet(MAX_SET_SIZE);
    EXPECT_TRUE(posSet.isValid() );
    EXPECT_EQ(vecSet, posSet);
    EXPECT_EQ(posSet, vecSet);
  }
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
