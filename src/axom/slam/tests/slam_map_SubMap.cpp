// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_SubMap.cpp
 *
 * \brief Unit tests for Slam's SubMap
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/SubMap.hpp"

namespace
{
namespace slam = axom::slam;

using SetBase = slam::Set<>;
using PositionType = SetBase::PositionType;
using ElementType = SetBase::ElementType;

using RangeSetType = slam::RangeSet<PositionType, ElementType>;

template <typename T>
using Map = slam::Map<SetBase, T>;

template <typename T, typename M>
using SubMap = slam::SubMap<SetBase, T, M>;

using OrderedSetType = axom::slam::OrderedSet<
  PositionType,
  ElementType,
  slam::policies::RuntimeSize<PositionType>,
  slam::policies::ZeroOffset<PositionType>,
  slam::policies::StrideOne<PositionType>,
  slam::policies::STLVectorIndirection<PositionType, ElementType>>;

static const double multFac = 1.0001;

static PositionType const MAX_SET_SIZE = 10;

}  // namespace

TEST(slam_map, construct_empty_subsetmap)
{
  using MapType = Map<int>;
  using SubMapType = SubMap<int, MapType>;
  SubMapType m;

  EXPECT_TRUE(m.isValid(true));
}

template <typename T>
T getValue(int idx)
{
  return static_cast<T>(idx * multFac);
}

//A struct to construct the Map for testing
template <typename T>
struct MapForTest
{
  MapForTest(int size)
    : set_data(size)
    , s(OrderedSetType::SetBuilder().size(size).data(&set_data))
    , m(&s)
  {
    SLIC_INFO("\nCreating set of size " << s.size());
    for(int i = 0; i < size; i++)
    {
      set_data[i] = 100 + i;
    }

    EXPECT_EQ(s.size(), size);
    EXPECT_TRUE(s.isValid());

    SLIC_INFO("\nCreating " << slam::util::TypeToString<T>::to_string()
                            << " map on the set ");

    EXPECT_TRUE(m.isValid());

    SLIC_INFO("\nSetting the elements.");
    for(PositionType idx = 0; idx < m.size(); ++idx)
    {
      m[idx] = getValue<T>(idx);
    }
  }

  std::vector<PositionType> set_data;
  OrderedSetType s;
  Map<T> m;
};

template <typename T>
bool constructAndTestSubMap()
{
  using MapType = Map<T>;
  using SubMapType = SubMap<T, MapType>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;
  OrderedSetType& s = mft.s;

  int submapOffset = 3;
  int submapSize = 5;
  {
    SLIC_INFO("\nCreating the Subset.");
    RangeSetType ss(submapOffset, submapOffset + submapSize);
    SubMapType ssm(&m, ss);
    EXPECT_TRUE(m.isValid(true));

    SLIC_INFO("\nChecking the elements.");
    for(PositionType idx = 0; idx < submapSize; ++idx)
    {
      EXPECT_EQ(ssm[idx], getValue<T>(submapOffset + idx));
    }

    EXPECT_TRUE(ssm.isValid(true));
  }

  {
    SLIC_INFO("\nCreating Subset 2");
    std::vector<PositionType> subset_indices_data(submapSize);
    for(int i = 0; i < submapSize; i++) subset_indices_data[i] = i * 2;
    OrderedSetType subset_indices = OrderedSetType::SetBuilder()  //
                                      .size(5)                    //
                                      .offset(0)                  //
                                      .data(&subset_indices_data);

    SubMapType ssm(&m, subset_indices);

    SLIC_INFO("\nChecking the elements.");
    for(PositionType idx = 0; idx < submapSize; ++idx)
    {
      EXPECT_EQ(ssm[idx], getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.value(idx), getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.index(idx), s[subset_indices[idx]]);
    }
  }
  return true;
}

TEST(slam_map, construct_int_submap)
{
  EXPECT_TRUE(constructAndTestSubMap<int>());
}

TEST(slam_map, construct_double_submap)
{
  EXPECT_TRUE(constructAndTestSubMap<double>());
}

template <typename T>
bool constructBySubMap()
{
  //This tests modifying the values in the original map via a Submap
  //Create a Map, then create a SubMap on the Map, negate all values covered by
  //the Submap, then check the values are negative in the SubMap range.
  using MapType = Map<T>;
  using SubMapType = SubMap<T, MapType>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;

  int submapOffset = 3;
  int submapSize = 5;

  SLIC_INFO("\nCreating the Subset.");
  RangeSetType ss(submapOffset, submapOffset + submapSize);
  SubMapType ssm(&m, ss);
  EXPECT_TRUE(m.isValid(true));

  SLIC_INFO("\nNegating elements");
  for(PositionType idx = 0; idx < submapSize; ++idx)
  {
    ssm[idx] = -getValue<T>(submapOffset + idx);
  }

  SLIC_INFO("\nChecking the elements.");
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    T val = getValue<T>(idx);
    if(idx >= submapOffset && idx < submapOffset + submapSize)
      val = -getValue<T>(idx);
    EXPECT_EQ(m[idx], val);
  }

  SLIC_INFO("\nChecking elements");

  return true;
}

TEST(slam_map, construct_with_int_submap)
{
  EXPECT_TRUE(constructBySubMap<int>());
}

TEST(slam_map, construct_with_double_submap)
{
  EXPECT_TRUE(constructBySubMap<double>());
}

//Some SubMap iterator tests in BivariateMap.
//TODO add more iterator tests here

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
