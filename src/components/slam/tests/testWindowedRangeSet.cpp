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
 * \file testWindowedRangeSet.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <iterator>
#include "gtest/gtest.h"

#include "slam/Utilities.hpp"
#include "slam/Set.hpp"
#include "slam/RangeSet.hpp"

typedef asctoolkit::meshapi::Set::PositionType         SetPosition;
typedef asctoolkit::meshapi::Set::ElementType          SetElement;

typedef asctoolkit::meshapi::policies::StrideOne<SetPosition> StrideOnePolicy;
typedef asctoolkit::meshapi::policies::NoIndirection<SetPosition,SetElement> NoIndirectionPolicy;
typedef asctoolkit::meshapi::policies::VirtualParentSubset SubsetPolicy;



typedef asctoolkit::meshapi::GenericRangeSet<StrideOnePolicy, NoIndirectionPolicy, SubsetPolicy> SetType;
typedef SetType::iterator             SetIterator;

static const SetPosition MAX_SET_SIZE = 20;


TEST(gtest_meshapi_windowed_range_set,construct_windowed_range_set)
{
  const SetElement lowerIndex = static_cast<SetElement>( .3 * MAX_SET_SIZE);
  const SetElement upperIndex = static_cast<SetElement>( .7 * MAX_SET_SIZE);

  SetType s(lowerIndex, upperIndex);

  EXPECT_TRUE(s.isValid(true));

  if(lowerIndex != upperIndex)
    EXPECT_FALSE(s.empty());


  std::cout << "Iterating through windowed range set of size " << s.size()
            << "\n\twith lower element " << lowerIndex << " (included in set)"
            << "\n\twith upper element " << upperIndex << " (not included in set)"
            << std::endl;
  const SetPosition expectedSize = upperIndex - lowerIndex;
  EXPECT_EQ(expectedSize, s.size() );


  std::cout << "\n --Testing random access -- operator[] and at() function" << std::endl;
  for(SetPosition pos = SetPosition(); pos < s.size(); ++pos)
  {
    SetElement expected = pos + lowerIndex;
    EXPECT_EQ(  expected, s[pos] )
      << "Random access subscript operator on windowed range should be translated by lowerIndex";

    EXPECT_EQ(  expected, s.at(pos) )
      << "Random access at() function on windowed range should be translated by lowerIndex";

    EXPECT_EQ(  s[pos],   s.at(pos) )
      << "Random access at() function should equal value of subscript operator";

    std::cout << "\t" << s[pos] << "\n";
  }

  std::cout << "\n --Using begin/end" << std::endl;
  for(SetIterator it = s.begin(), itEnd = s.end(); it != itEnd; ++it)
  {
    SetPosition position = std::distance(s.begin(), it);
    SetElement expected = position + lowerIndex;
    EXPECT_EQ( expected, *it )
      << "Iterator dereference should be equal to its translated position in the windowed range set";
    std::cout << "\t" << *it << "\n";
  }

  std::cout << "\n --Using random access on invalid address -- Note: We are testing for the expected failures." << std::endl;
#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH( s.at(upperIndex),   "") << "tried to access out of range element (" << upperIndex << ")";
  ASSERT_DEATH( s.at(MAX_SET_SIZE), "") << "tried to access out of range element (" << MAX_SET_SIZE << ")";
#else
  std::cout << "Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
#endif

  std::cout << "--\ndone." << std::endl;
}

TEST(gtest_meshapi_windowed_range_set,test_windowed_range_set_parents)
{
  const SetElement lowerIndex = static_cast<SetElement>( .3 * MAX_SET_SIZE);
  const SetElement upperIndex = static_cast<SetElement>( .7 * MAX_SET_SIZE);

  std::cout << "\n-- Generating a parent set, one (windowed) subset and one non-windowed subset and checking validity" << std::endl;
  SetType parentSet(MAX_SET_SIZE);
  SetType childSet(lowerIndex, upperIndex);
  childSet.parentSet() = &parentSet;

  SetType nonChildSet(lowerIndex, upperIndex);

  EXPECT_TRUE(parentSet.isValid(true));
  EXPECT_TRUE(childSet.isValid(true));
  EXPECT_TRUE(nonChildSet.isValid(true));

  std::cout << "\n-- Checking that the child is a subset, but not the parent or the non-child windowed set." << std::endl;
  EXPECT_FALSE(parentSet.isSubset());
  EXPECT_TRUE(childSet.isSubset());
  EXPECT_FALSE(nonChildSet.isSubset());

  std::cout << "\n-- Checking that the child set's parent is equal to the parent set (according to the equality operator==)." << std::endl;
  EXPECT_EQ(parentSet, *childSet.parentSet());
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
