// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_Map.cpp
 *
 * \brief Unit tests for Slam's Map
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Map.hpp"

namespace
{
namespace slam = axom::slam;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using SetType = slam::RangeSet<SetPosition, SetElement>;
using BaseSet = slam::Set<SetPosition, SetElement>;
using IntMap = slam::Map<BaseSet, int>;
using RealMap = slam::Map<BaseSet, double>;

static SetPosition const MAX_SET_SIZE = 10;

template <int S>
using CompileTimeStrideType = slam::policies::CompileTimeStride<int, S>;

using RunTimeStrideType = slam::policies::RuntimeStride<int>;
using OneStrideType = slam::policies::StrideOne<int>;

}  // end anonymous namespace

TEST(slam_map, construct_empty_map)
{
  IntMap m;

  EXPECT_TRUE(m.isValid(true));
}

template <typename T>
bool constructAndTestMap()
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("Creating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string()
                        << " map on the set ");

  slam::Map<BaseSet, T> m(&s);
  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Setting the elements.");
  double multFac = 1.0001;
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    m[idx] = static_cast<T>(idx * multFac);
  }

  SLIC_INFO("Checking the elements.");
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    EXPECT_EQ(m[idx], static_cast<T>(idx * multFac));
    EXPECT_EQ(m(idx), static_cast<T>(idx * multFac));
  }

  EXPECT_TRUE(m.isValid(true));

  m.print();
  return true;
}

TEST(slam_map, construct_int_map) { EXPECT_TRUE(constructAndTestMap<int>()); }

TEST(slam_map, construct_double_map)
{
  EXPECT_TRUE(constructAndTestMap<double>());
}

TEST(slam_map, out_of_bounds)
{
  int defaultElt = 2;

  SetType s(MAX_SET_SIZE);
  IntMap m(&s, defaultElt);

  SLIC_INFO("Testing Map element access -- in bounds");
  for(auto idx = 0; idx < m.size(); ++idx) EXPECT_EQ(defaultElt, m[idx]);

  // Test out of bounds
  SLIC_INFO("Testing Map element access "
            << "-- out of bounds access; Expecting the test to fail");
#ifdef AXOM_DEBUG
  EXPECT_DEATH_IF_SUPPORTED(m[-1], "")
    << " Accessed element -1 of Map -- out of bounds";
  EXPECT_DEATH_IF_SUPPORTED(m[m.size()], "")
    << " Accessed element " << m.size() << " of Map -- out of bounds";

#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

TEST(slam_map, map_builder)
{
  SLIC_INFO("Testing construction of Map using MapBuilders");

  using DataType = double;
  using MapType = slam::Map<slam::Set<>, DataType>;
  using MapBuilder = MapType::MapBuilder;

  MapType m(MapBuilder().set(&MapType::s_nullSet));
  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(m.size(), 0);
  EXPECT_EQ(m.stride(), 1);

  SetType s(MAX_SET_SIZE);
  std::vector<DataType> data_arr(s.size());
  for(auto i = 0u; i < data_arr.size(); ++i)
    data_arr[i] = static_cast<DataType>(i * 1.01);

  MapType m2(MapBuilder().set(&s).data(data_arr.data()));
  EXPECT_TRUE(m2.isValid());
  EXPECT_EQ(m2.size(), s.size());
  EXPECT_EQ(m2.stride(), 1);
  for(auto i = 0u; i < data_arr.size(); ++i) EXPECT_EQ(m2[i], data_arr[i]);
}

template <typename T, typename StrideType>
void constructAndTestMapWithStride(int stride)
{
  SetType s(MAX_SET_SIZE);

  SLIC_INFO("\nCreating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("\nCreating " << slam::util::TypeToString<T>::to_string()
                          << " map with stride " << stride << " on the set ");

  axom::slam::Map<BaseSet, T, StrideType> m(&s, 0, stride);
  EXPECT_TRUE(m.isValid());

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements.");
  double multFac = 100.0001;
  double multFac2 = 1.010;
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    for(auto idx2 = 0; idx2 < stride; ++idx2)
    {
      m(idx, idx2) = static_cast<T>(idx * multFac + idx2 * multFac2);
    }
  }

  EXPECT_TRUE(m.isValid());

  SLIC_INFO("\nChecking the elements.");
  for(auto idx = 0; idx < m.size(); ++idx)
  {
    for(auto idx2 = 0; idx2 < stride; ++idx2)
    {
      EXPECT_EQ(m(idx, idx2), static_cast<T>(idx * multFac + idx2 * multFac2));
    }
  }

  EXPECT_TRUE(m.isValid(true));
}

TEST(slam_map, construct_double_map_with_stride)
{
  constructAndTestMapWithStride<double, RunTimeStrideType>(1);
  constructAndTestMapWithStride<double, RunTimeStrideType>(2);
  constructAndTestMapWithStride<double, RunTimeStrideType>(3);

  constructAndTestMapWithStride<double, CompileTimeStrideType<1>>(1);
  constructAndTestMapWithStride<double, CompileTimeStrideType<2>>(2);
  constructAndTestMapWithStride<double, CompileTimeStrideType<3>>(3);

  constructAndTestMapWithStride<double, OneStrideType>(1);
}

TEST(slam_map, iterate)
{
  using IterType = RealMap::iterator;

  SLIC_INFO("Testing iterator access");

  SetType s(MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating '" << slam::util::TypeToString<double>::to_string()
                         << "' map on the set ");
  RealMap m(&s);
  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Setting the elements.");
  double multFac = 1.0001;
  {
    int idx = 0;
    for(IterType iter = m.begin(); iter != m.end(); iter++)
    {
      *iter = static_cast<double>(idx * multFac);
      idx++;
    }
  }

  SLIC_INFO("Checking the elements by iterator.");
  EXPECT_EQ(RealMap::iterator(0, &m), m.begin());

  //iter++ access
  {
    int idx = 0;
    for(IterType iter = m.begin(); iter != m.end(); iter++)
    {
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
      idx++;
    }
    EXPECT_EQ(idx, m.size());
  }

  //iter+n access
  {
    IterType beginIter = m.begin();
    for(int idx = 0; idx < m.size(); ++idx)
    {
      IterType iter = beginIter + idx;
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
    }
  }

  //iter-n access
  {
    IterType endIter = m.end();
    for(int idx = 1; idx <= m.size(); ++idx)
    {
      IterType iter = endIter - idx;
      EXPECT_EQ(*iter, static_cast<double>((m.size() - idx) * multFac));
    }
  }

  //iter+=n access
  {
    for(int idx = 0; idx < m.size(); idx++)
    {
      IterType iter = m.begin();
      iter += idx;
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
    }
  }

  //iter-=n access
  {
    for(int idx = 1; idx <= m.size(); ++idx)
    {
      IterType iter = m.end();
      iter -= idx;
      EXPECT_EQ(*iter, static_cast<double>((m.size() - idx) * multFac));
    }
  }

  //iter[n] access
  {
    IterType beginIter = m.begin();
    for(int idx = 0; idx < m.size(); idx++)
    {
      EXPECT_EQ(beginIter[idx], static_cast<double>(idx * multFac));
    }
  }

  //iter1 - iter2
  {
    IterType beginIter = m.begin();
    IterType endIter = m.end();
    for(int idx = 0; idx < m.size(); idx++)
    {
      IterType iter = beginIter + idx;
      EXPECT_EQ(idx, iter - beginIter);
      EXPECT_EQ(idx - m.size(), iter - endIter);
    }
  }

  EXPECT_TRUE(m.isValid(true));
}

template <typename StrideType>
void constructAndTestMapIteratorWithStride(int stride)
{
  using RealMap = slam::Map<slam::Set<>, double, StrideType>;
  using MapIterator = typename RealMap::MapIterator;

  SetType s(MAX_SET_SIZE);

  SLIC_INFO("\nCreating set of size " << s.size());

  EXPECT_EQ(s.size(), MAX_SET_SIZE);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating double map with stride " << stride << " on the set ");

  RealMap m(&s, 0, stride);
  EXPECT_TRUE(m.isValid());

  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements using iterator.");
  double multFac = 100.0001;
  double multFac2 = 1.010;
  {
    int idx = 0;
    for(MapIterator iter = m.begin(); iter != m.end(); ++iter)
    {
      for(auto idx2 = 0; idx2 < iter.numComp(); ++idx2)
      {
        iter(idx2) = static_cast<double>(idx * multFac + idx2 * multFac2);
      }
      ++idx;
    }
  }

  EXPECT_TRUE(m.isValid());

  SLIC_INFO("Checking the elements by iterator.");

  //iter++ access
  {
    int idx = 0;
    for(MapIterator iter = m.begin(); iter != m.end(); iter++)
    {
      EXPECT_EQ(*iter, static_cast<double>(idx * multFac));
      for(auto idx2 = 0; idx2 < iter.numComp(); ++idx2)
      {
        EXPECT_EQ(iter(idx2),
                  static_cast<double>(idx * multFac + idx2 * multFac2));
      }
      idx++;
    }
    EXPECT_EQ(idx, m.size());
  }
}

TEST(slam_map, iterate_with_stride)
{
  constructAndTestMapIteratorWithStride<OneStrideType>(1);

  constructAndTestMapIteratorWithStride<RunTimeStrideType>(1);
  constructAndTestMapIteratorWithStride<RunTimeStrideType>(2);
  constructAndTestMapIteratorWithStride<RunTimeStrideType>(3);

  constructAndTestMapIteratorWithStride<CompileTimeStrideType<1>>(1);
  constructAndTestMapIteratorWithStride<CompileTimeStrideType<2>>(2);
  constructAndTestMapIteratorWithStride<CompileTimeStrideType<3>>(3);

  SLIC_INFO("Done");
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
