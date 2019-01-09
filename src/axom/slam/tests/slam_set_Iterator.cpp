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


/*
 * \file slam_Map.cpp
 *
 * \brief Unit tests for Slam's Map
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/IndirectionSet.hpp"

namespace
{
namespace slam = axom::slam;

using SetPosition = slam::PositionType;
using SetElement = slam::ElementType;

using SetBase = slam::Set<SetPosition, SetElement>;

using RangeSet = slam::RangeSet<SetPosition, SetElement>;
using VectorSet = slam::VectorIndirectionSet<SetPosition, SetElement>;


} // end anonymous namespace

TEST(slam_set_iterator,range_set)
{
  using SetType = RangeSet;
  SetType set(10, 20);

  SetType::iterator it1 = set.begin();
  SetType::iterator it2 = set.begin();
  EXPECT_EQ(it1, it2);

  ++it2;
  EXPECT_NE(it1, it2);
  EXPECT_GT(it2, it1);
  EXPECT_LT(it1, it2);

  SetType::iterator it3 = set.end();
  EXPECT_NE(it1, it3);

  it2 = it3;
  EXPECT_EQ(it2, it3);

  --it2;
  --it3;
  EXPECT_EQ(it2, it3);
}


TEST(slam_set_iterator,vector_set)
{
  using SetType = VectorSet;

  int sz = 10;
  VectorSet::IndirectionBufferType vec(sz);
  SetType set( SetType::SetBuilder()
               .size(sz)
               .data( &vec) );

  SetType set2;
  set2 = set;

  for( auto p : set.positions() )
  {
    set[p] = 2 * sz - p;
    set2[p] = 3 * sz - p;
  }

  SLIC_INFO("set 1");
  for(auto& el : set )
  {
    SLIC_INFO( el );
  }

  SLIC_INFO("set 2");
  for(auto& el : set2 )
  {
    SLIC_INFO( el );
  }

  SetType::iterator it1 = set.begin();
  SetType::iterator it2 = set.begin();
  EXPECT_EQ(it1, it2);

  ++it2;
  EXPECT_NE(it1, it2);
  EXPECT_GT(it2, it1);
  EXPECT_LT(it1, it2);

  SetType::iterator it3 = set.end();
  EXPECT_NE(it1, it3);

  it2 = it3;
  EXPECT_EQ(it2, it3);

  --it2;
  --it3;
  EXPECT_EQ(it2, it3);
}


TEST(slam_set_iterator,equality_non_const_range_set)
{
  using SetType = RangeSet;

  SetType set(10, 20);

  // check equality of iterators
  {
    SetType::iterator b1 = set.begin();
    SetType::iterator b2 = set.begin();
    SetType::iterator e1 = set.end();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);
  }

  // mix const and non-const iterators
  // to non-const set
  {
    SetType::iterator b1 = set.begin();
    SetType::iterator e1 = set.end();

    SetType::const_iterator b2 = set.begin();
    SetType::const_iterator b3 = set.begin();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);
    EXPECT_EQ(b2,b3);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);

    EXPECT_NE(b2,e1);
    EXPECT_NE(e1,b2);
  }
}


TEST(slam_set_iterator,equality_non_const_vector_set)
{
  using SetType = VectorSet;

  int sz = 10;
  VectorSet::IndirectionBufferType vec(sz);
  SetType set( SetType::SetBuilder()
               .size(sz)
               .data( &vec) );

  // check equality of iterators
  {
    SetType::iterator b1 = set.begin();
    SetType::iterator b2 = set.begin();
    SetType::iterator e1 = set.end();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);
  }

  // mix const and non-const iterators
  // to non-const set
  {
    SetType::iterator b1 = set.begin();
    SetType::iterator e1 = set.end();

    SetType::const_iterator b2 = set.begin();
    SetType::const_iterator b3 = set.begin();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);
    EXPECT_EQ(b2,b3);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);

    EXPECT_NE(b2,e1);
    EXPECT_NE(e1,b2);
  }
}

TEST(slam_set_iterator,equality_const_range_set)
{
  using SetType = const RangeSet;
  SetType set(10, 20);

  // Invalid to create a non const iterator to a const set
  // Uncommenting will result in compilation errors of the form:
  //   error: conversion const to non-const
  {
    //SetType::iterator bInvalid = set.begin();
    //SetType::iterator eInvalid = set.begin();
  }

  // Compare const iterators
  {
    auto b1 = set.begin();
    SetType::const_iterator b2 = set.begin();

    SetType::const_iterator e1 = set.end();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);
  }
}

TEST(slam_set_iterator,equality_const_vector_set)
{
  using SetType = const VectorSet;

  int sz = 10;
  VectorSet::IndirectionBufferType vec(sz);
  SetType set( SetType::SetBuilder()
               .size(sz)
               .data( &vec) );

  // Invalid to create a non const iterator to a const set
  // Uncommenting will result in compilation errors of the form:
  //   error: conversion const to non-const
  {
    //SetType::iterator bInvalid = set.begin();
    //SetType::iterator eInvalid = set.begin();
  }

  // Compare const iterators
  {
    auto b1 = set.begin();
    SetType::const_iterator b2 = set.begin();

    SetType::const_iterator e1 = set.end();

    EXPECT_EQ(b1,b2);
    EXPECT_EQ(b2,b1);

    EXPECT_NE(b1,e1);
    EXPECT_NE(e1,b1);
  }
}


#include <typeinfo>

TEST(vector, types)
{

  using Vec = std::vector<int>;

  Vec v(5);
  {
    SLIC_INFO ( typeid( v[2]).name() );
    SLIC_INFO ( typeid( &v[2]).name() );

    SLIC_INFO ( typeid( v.begin() ).name() );

    {
      auto it = v.begin();
      SLIC_INFO ( typeid( it[2] ).name() );
    }

    {
      auto it = static_cast<const Vec>(v).begin();
      SLIC_INFO ( typeid( it[2] ).name() );
    }
  }

  SLIC_INFO(" -------- ");

  const Vec cv(5);
  {
    SLIC_INFO ( typeid( cv[2]).name() );
    SLIC_INFO ( typeid( &cv[2]).name() );

    SLIC_INFO ( typeid( cv.begin() ).name() );

    {
      auto it = cv.begin();
      SLIC_INFO ( typeid( it[2] ).name() );
    }

    {
      auto it = static_cast<const Vec>(cv).begin();
      SLIC_INFO ( typeid( it[2] ).name() );
    }
  }

}




//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

  int result = RUN_ALL_TESTS();

  return result;
}
