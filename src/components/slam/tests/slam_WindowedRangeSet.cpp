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
 * \file slam_WindowedRangeSet.cpp
 *
 * \brief Unit tests for Slam's WindowedRangeSet
 *        A windowed range set is a contiguous range of integers from lo to high
 */

#include <iterator>
#include "gtest/gtest.h"

#include "common/config.hpp"

#include "slic/slic.hpp"

#include "slam/Utilities.hpp"
#include "slam/Set.hpp"
#include "slam/RangeSet.hpp"

typedef axom::slam::Set::PositionType                                                   SetPosition;
typedef axom::slam::Set::ElementType                                                    SetElement;

typedef axom::slam::policies::StrideOne<SetPosition>                                    StrideOnePolicy;
typedef axom::slam::policies::NoIndirection<SetPosition,SetElement>                     NoIndirectionPolicy;
typedef axom::slam::policies::VirtualParentSubset                                       SubsetPolicy;



typedef axom::slam::GenericRangeSet<StrideOnePolicy, NoIndirectionPolicy, SubsetPolicy> SetType;
static const SetPosition MAX_SET_SIZE = 20;


TEST(gtest_slam_windowed_range_set,construct_windowed_range_set)
{
  const SetElement lowerIndex = static_cast<SetElement>( .3 * MAX_SET_SIZE);
  const SetElement upperIndex = static_cast<SetElement>( .7 * MAX_SET_SIZE);

  SetType s(lowerIndex, upperIndex);

  EXPECT_TRUE(s.isValid(true));

  if(lowerIndex != upperIndex)
    EXPECT_FALSE(s.empty());


  SLIC_INFO( "Iterating through windowed range set of size "
      << s.size()
      << "\n\twith lower element " << lowerIndex << " (included in set)"
      << "\n\twith upper element " << upperIndex << " (not included in set)");
  const SetPosition expectedSize = upperIndex - lowerIndex;
  EXPECT_EQ(expectedSize, s.size() );

  SLIC_INFO("Testing random access -- operator[] and at() function");
  {
    std::stringstream sstr;
    for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
    {
      SetElement expected = pos + lowerIndex;
      EXPECT_EQ(  expected, s[pos] )
        << "Random access subscript operator on windowed range should be translated by lowerIndex";

      EXPECT_EQ(  expected, s.at(pos) )
        << "Random access at() function on windowed range should be translated by lowerIndex";

      EXPECT_EQ(  s[pos],   s.at(pos) )
        << "Random access at() function should equal value of subscript operator";

      sstr << "\t" << s[pos] << "\n";
    }
    SLIC_INFO(sstr.str());
  }

#ifdef AXOM_USE_BOOST
  SLIC_INFO("Testing iterator access");
  {
    std::stringstream sstr;

    typedef SetType::iterator SetIterator;
    for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
    {
      SetPosition position = std::distance(s.begin(), it);
      SetElement expected = position + lowerIndex;
      EXPECT_EQ( expected, *it )
        << "Iterator dereference should be equal to its translated position in the windowed range set";
      sstr << "\t" << *it << "\n";
    }
    SLIC_INFO(sstr.str());
  }
#endif

  SLIC_INFO("Using random access on invalid address -- Note: We are testing for the expected failures.");
#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( s.at(upperIndex),   "") << "tried to access out of range element (" << upperIndex << ")";
  EXPECT_DEATH_IF_SUPPORTED( s.at(MAX_SET_SIZE), "") << "tried to access out of range element (" << MAX_SET_SIZE << ")";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("done.");
}

TEST(gtest_slam_windowed_range_set,test_windowed_range_set_parents)
{
  const SetElement lowerIndex = static_cast<SetElement>( .3 * MAX_SET_SIZE);
  const SetElement upperIndex = static_cast<SetElement>( .7 * MAX_SET_SIZE);

  SLIC_INFO("Generating a parent set, one (windowed) subset and one non-windowed subset and checking validity");
  SetType parentSet(MAX_SET_SIZE);
  SetType childSet(lowerIndex, upperIndex);

  childSet.parentSet() = &parentSet;

  SetType nonChildSet(lowerIndex, upperIndex);

  EXPECT_TRUE(parentSet.isValid(true));
  EXPECT_TRUE(childSet.isValid(true));
  EXPECT_TRUE(nonChildSet.isValid(true));

  SLIC_INFO("Checking that the child is a subset, but not the parent or the non-child windowed set.");
  EXPECT_FALSE(parentSet.isSubset());
  EXPECT_TRUE(childSet.isSubset());
  EXPECT_FALSE(nonChildSet.isSubset());

  SLIC_INFO("Checking that the child set's parent is equal to the parent set (according to the equality operator==).");
  EXPECT_EQ(parentSet, *childSet.parentSet());
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
