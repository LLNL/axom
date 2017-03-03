/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * \file slam_RangeSet.cpp
 *
 */

#include <iterator>
#include <sstream>       // For stringstream

#include "gtest/gtest.h"

#include "common/config.hpp"

#include "slic/slic.hpp"

#include "slam/Utilities.hpp"
#include "slam/RangeSet.hpp"


typedef axom::slam::RangeSet  SetType;
typedef SetType::PositionType       SetPosition;
typedef SetType::ElementType        SetElement;

static const SetPosition MAX_SET_SIZE = 10;


TEST(gtest_slam_range_set,construct_range_set)
{

  SetType s(MAX_SET_SIZE);

  EXPECT_TRUE(s.isValid());

  if(MAX_SET_SIZE > SetPosition())
    EXPECT_FALSE(s.empty());

  SLIC_INFO("Iterating through set of size " << s.size());
  EXPECT_EQ(s.size(), MAX_SET_SIZE);

#ifdef AXOM_USE_BOOST
  SLIC_INFO("Using begin/end");
  {
    std::stringstream sstr;
    typedef SetType::iterator SetIterator;
    for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      EXPECT_EQ( std::distance(s.begin(), it), *it )
        << "Iterator dereference should be equal to its position in the set";
      sstr << "\t" << *it << "\n";
    }

    SLIC_INFO("Element of slam set using iterators:\n" << sstr.str());
  }
#endif

  SLIC_INFO("Using random access -- operator[]");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement elt = static_cast<SetElement>(pos);
      EXPECT_EQ(elt,s[pos])
        << "Random access iterator dereference to equal its position in the set";
      sstr << "\t" << s[pos] << "\n";
    }
    SLIC_INFO("Element of slam set using operator[]:\n" << sstr.str());
  }

  SLIC_INFO("Using checked random access -- at()");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement elt = static_cast<SetElement>(pos);
      EXPECT_EQ(elt,s.at(pos))
        << "Expected checked random access iterator dereference to equal its position in the set";

      sstr << "\t" << s[pos] << "\n";
    }
    SLIC_INFO("Element of slam set using at():\n" << sstr.str());
  }

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  SLIC_INFO("Using checked random access -- at() with invalid address");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(s.at(MAX_SET_SIZE),"") << "tried to access out of range element";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif


  SLIC_INFO("done.");

}

TEST(gtest_slam_range_set,test_range_set_out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized set-- code is expected to assert and die.");


  SetType s(MAX_SET_SIZE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s[MAX_SET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
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
