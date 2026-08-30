// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/BivariateMap.hpp"

#include <type_traits>

namespace
{
namespace slam = axom::slam;

using SetBase = slam::Set<>;
using PositionType = SetBase::PositionType;
using ElementType = SetBase::ElementType;

using RangeSetType = slam::RangeSet<PositionType, ElementType>;

template <typename T>
using SuperMap = slam::Map<T, SetBase>;

using OrderedSetType =
  axom::slam::OrderedSet<PositionType,
                         ElementType,
                         slam::policies::RuntimeSize<PositionType>,
                         slam::policies::ZeroOffset<PositionType>,
                         slam::policies::StrideOne<PositionType>,
                         slam::policies::STLVectorIndirection<PositionType, ElementType>>;

constexpr double multFac = 1.0001;

PositionType const MAX_SET_SIZE = 10;

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
    SLIC_INFO("Initializing set of size " << s.size() << " and map on the set ");

    for(auto i : s.positions())
    {
      s[i] = 100 + i;
      m[i] = getValue<T>(i);
    }

    EXPECT_TRUE(s.isValid());
    EXPECT_TRUE(m.isValid());
  }

  std::vector<PositionType> set_data;
  OrderedSetType s;
  SuperMap<T> m;
};

}  // namespace

TEST(slam_map, construct_empty_subsetmap)
{
  slam::SubMap<SuperMap<int>, RangeSetType> m;
  EXPECT_TRUE(m.isValid(true));
}

template <typename T>
bool constructAndTestSubMap()
{
  using MapType = SuperMap<T>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;
  OrderedSetType& s = mft.s;

  int submapOffset = 3;
  int submapSize = 5;
  {
    using SubMapType = slam::SubMap<MapType, RangeSetType>;

    SLIC_INFO("Creating the Subset.");
    RangeSetType ss(submapOffset, submapOffset + submapSize);
    SubMapType ssm(&m, ss);
    EXPECT_TRUE(ssm.isValid(true));

    auto flatIter = ssm.begin();
    EXPECT_EQ(flatIter.operator->(), &ssm[0]);

    SLIC_INFO("Checking the elements.");
    for(auto idx = 0; idx < ssm.size(); ++idx)
    {
      auto expVal = getValue<T>(submapOffset + idx);
      EXPECT_EQ(expVal, ssm[idx]);
    }

    SLIC_INFO("Checking the elements using SubMap iterator.");
    int cnt = 0;
    for(auto it = ssm.begin(); it != ssm.end(); ++it, ++cnt)
    {
      // Check iterator's .index() function
      {
        auto subMapElt = it.index();

        auto setIdx = ss[cnt];
        auto expSetElt = s[setIdx];
        EXPECT_EQ(expSetElt, subMapElt);
      }

      // Check iterator's value access functions
      {
        auto expectedValue = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expectedValue, *it);
        EXPECT_EQ(expectedValue, it[0]);

        auto expVal = getValue<T>(submapOffset + cnt);
        EXPECT_EQ(expVal, expectedValue);
      }
    }
  }

  {
    using SubMapType = slam::SubMap<MapType, OrderedSetType>;

    SLIC_INFO("Creating Subset 2");
    std::vector<PositionType> subset_indices_data(submapSize);
    for(int i = 0; i < submapSize; i++)
    {
      subset_indices_data[i] = i * 2;
    }
    OrderedSetType subset_indices = OrderedSetType::SetBuilder().size(5).data(&subset_indices_data);

    SubMapType ssm(&m, subset_indices);

    SLIC_INFO("Checking the elements.");
    for(PositionType idx = 0; idx < submapSize; ++idx)
    {
      EXPECT_EQ(ssm[idx], getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.value(idx), getValue<T>(subset_indices[idx]));
      EXPECT_EQ(ssm.index(idx), s[subset_indices[idx]]);
    }

    SLIC_INFO("Checking the elements using SubMap range iterator.");
    {
      int cnt = 0;
      for(auto it = ssm.set_begin(); it != ssm.set_end(); ++it, ++cnt)
      {
        auto expectedValue = getValue<T>(subset_indices[cnt]);
        EXPECT_EQ(expectedValue, (*it)[0]);
        EXPECT_EQ(expectedValue, it(0));
        EXPECT_EQ(expectedValue, it.value(0));
        EXPECT_EQ(it.index(), s[subset_indices[cnt]]);
      }
      EXPECT_EQ(cnt, subset_indices.size());
    }
  }
  return true;
}

TEST(slam_map, construct_int_submap) { EXPECT_TRUE(constructAndTestSubMap<int>()); }

TEST(slam_map, construct_double_submap) { EXPECT_TRUE(constructAndTestSubMap<double>()); }

template <typename T>
bool constructBySubMap()
{
  //This tests modifying the values in the original map via a Submap
  //Create a Map, then create a SubMap on the Map, negate all values covered by
  //the Submap, then check the values are negative in the SubMap range.
  using MapType = SuperMap<T>;
  using SubMapType = slam::SubMap<MapType, RangeSetType>;

  MapForTest<T> mft(MAX_SET_SIZE);
  MapType& m = mft.m;

  int submapOffset = 3;
  int submapSize = 5;

  SLIC_INFO("Creating the Subset.");
  RangeSetType ss(submapOffset, submapOffset + submapSize);
  SubMapType ssm(&m, ss);
  EXPECT_TRUE(m.isValid(true));

  SLIC_INFO("Negating elements");
  for(PositionType idx = 0; idx < submapSize; ++idx)
  {
    ssm[idx] = -getValue<T>(submapOffset + idx);
  }

  SLIC_INFO("Checking the elements.");
  for(PositionType idx = 0; idx < m.size(); ++idx)
  {
    T val = getValue<T>(idx);
    if(idx >= submapOffset && idx < submapOffset + submapSize)
    {
      val = -getValue<T>(idx);
    }
    EXPECT_EQ(m[idx], val);
  }

  SLIC_INFO("Checking elements");

  return true;
}

TEST(slam_map, submap_of_submap)
{
  // A SubMap is itself a map, so we can take a SubMap from it
  using ProductSetType = typename slam::ProductSet<RangeSetType, RangeSetType>::ConcreteSet;
  using BMapType = slam::BivariateMap<double, ProductSetType>;
  using RowSubMap = typename BMapType::SubMapType;
  using NestedSubMap = slam::SubMap<RowSubMap, typename RowSubMap::IndexSetType>;

  static_assert(slam::SubMappable<BMapType>);
  static_assert(slam::SubMappable<RowSubMap>, "a SubMap can serve as a super-map");
  static_assert(slam::SubMappable<NestedSubMap>, "and composition does not bottom out");

  constexpr PositionType NROWS = 4, NCOLS = 5;
  RangeSetType rows(NROWS), cols(NCOLS);
  ProductSetType prod(&rows, &cols);
  BMapType bmap(prod, 0.0);

  // each value encodes its own coordinate, so a misindex cannot alias
  for(PositionType i = 0; i < NROWS; ++i)
  {
    for(PositionType j = 0; j < NCOLS; ++j)
    {
      bmap(i, j) = 100.0 * i + j;
    }
  }

  RowSubMap row2 = bmap(2);
  ASSERT_EQ(NCOLS, row2.size());

  // columns 1..3 of row 2
  constexpr PositionType FIRST = 1, COUNT = 3;
  auto inner = typename RowSubMap::IndexSetType::SetBuilder().size(COUNT).offset(FIRST);
  NestedSubMap mid(&row2, inner);

  ASSERT_EQ(COUNT, mid.size());
  for(PositionType k = 0; k < mid.size(); ++k)
  {
    EXPECT_EQ(200.0 + (k + FIRST), mid[k]);
  }

  // element iteration
  double sum = 0.0;
  PositionType visited = 0;
  for(auto it = mid.begin(); it != mid.end(); ++it, ++visited)
  {
    sum += *it;
  }
  EXPECT_EQ(COUNT, visited);
  EXPECT_EQ(201.0 + 202.0 + 203.0, sum);

  // Range iteration. This is a regression test for an index-space bug in
  // SubMap::RangeIterator::advance(): stepping the parent iterator by a
  // difference computed from the parent's flatIndex() only works when the
  // parent is a Map or BivariateMap.
  double rangeSum = 0.0;
  PositionType rangeVisited = 0;
  for(auto it = mid.set_begin(); it != mid.set_end(); ++it, ++rangeVisited)
  {
    rangeSum += (*it)[0];
  }
  EXPECT_EQ(COUNT, rangeVisited);
  EXPECT_EQ(201.0 + 202.0 + 203.0, rangeSum);

  // index() projects exactly one level: a nested subset index becomes a
  // position in the parent SubMap's index set, not a bivariate coordinate.
  for(PositionType k = 0; k < mid.size(); ++k)
  {
    EXPECT_EQ(row2.set()->at(k + FIRST), mid.index(k));
  }
  // ... and projecting once more recovers the coordinate the value encodes
  for(PositionType k = 0; k < mid.size(); ++k)
  {
    const auto coordinate = row2.index(k + FIRST);
    EXPECT_EQ(2, coordinate.first);
    EXPECT_EQ(k + FIRST, coordinate.second);
  }
}

TEST(slam_map, construct_with_int_submap) { EXPECT_TRUE(constructBySubMap<int>()); }

TEST(slam_map, construct_with_double_submap) { EXPECT_TRUE(constructBySubMap<double>()); }

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
#ifdef AXOM_DEBUG
  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
#endif

  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
