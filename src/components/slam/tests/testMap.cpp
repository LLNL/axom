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
 * testSet.cxx
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <iterator>
#include "gtest/gtest.h"

#include "slam/Utilities.hpp"
#include "slam/RangeSet.hpp"
#include "slam/Map.hpp"


typedef asctoolkit::slam::RangeSet    SetType;
typedef asctoolkit::slam::Map<int>    IntMap;
typedef asctoolkit::slam::Map<double> RealMap;

typedef SetType::PositionType         PositionType;
typedef SetType::ElementType          ElementType;

typedef SetType::iterator             SetIterator;
static PositionType const MAX_SET_SIZE = 10;

TEST(gtest_slam_map,construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
}

template<typename T>
bool constructAndTestMap()
{
  SetType s(MAX_SET_SIZE);

  std::cout << "\nCreating set of size " << s.size() << std::endl;

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  std::cout << "\nCreating " << asctoolkit::slam::util::TypeToString<T>::to_string() << " map on the set " << std::endl;
  asctoolkit::slam::Map<T> m(&s);
  EXPECT_TRUE(m.isValid());

  std::cout << "\nSetting the elements.";
  double multFac = 1.0001;
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    m[idx] = static_cast<T>(idx * multFac);
  }

  std::cout << "\nChecking the elements.";
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    EXPECT_EQ(m[idx], static_cast<T>(idx * multFac) );
  }

  EXPECT_TRUE(m.isValid(true));

  return true;
}

TEST(gtest_slam_map,construct_int_map)
{
  EXPECT_TRUE( constructAndTestMap<int>() );
}

TEST(gtest_slam_map,construct_double_map)
{
  EXPECT_TRUE( constructAndTestMap<double>());
}
