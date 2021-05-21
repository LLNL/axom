// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_map_BivariateMap.cpp
 *
 * \brief Unit tests for Slam's Bivariate Map
 */

#include <iterator>
#include "gtest/gtest.h"

#include "axom/slic.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/BivariateMap.hpp"
#include "axom/slam/RelationSet.hpp"
#include "axom/slam/ProductSet.hpp"
#include "axom/slam/StaticRelation.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;
using SetType = slam::RangeSet<SetPosition, SetElement>;

using StrideOneType = policies::StrideOne<SetPosition>;

template <unsigned int S>
using CompileTimeStrideType = policies::CompileTimeStride<SetPosition, S>;

using RuntimeStrideType = policies::RuntimeStride<SetPosition>;

using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;
using VariableCardinality =
  policies::VariableCardinality<SetPosition, STLIndirection>;

using RelationType =
  slam::StaticRelation<SetPosition, SetElement, VariableCardinality, STLIndirection, SetType, SetType>;

using BivariateSetType = axom::slam::BivariateSet<SetPosition, SetElement>;
using ProductSetType = axom::slam::ProductSet<SetPosition, SetElement>;
using RelationSetType = axom::slam::RelationSet<RelationType>;

template <typename T, typename S>
using BivariateMapType = axom::slam::BivariateMap<SetType, T, S>;

static const SetPosition MAX_SET_SIZE1 = 10;
static const SetPosition MAX_SET_SIZE2 = 15;

static double const multFac3 = 0000.1;
static double const multFac1 = 1000.0;
static double const multFac2 = 0010.0;

}  // end anonymous namespace

TEST(slam_bivariate_map, construct_empty_map)
{
  BivariateMapType<int, StrideOneType> m;

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(m.totalSize(), 0);
  EXPECT_EQ(m.firstSetSize(), 0);
  EXPECT_EQ(m.secondSetSize(), 0);
}

template <typename T>
inline T getVal(SetPosition idx1, SetPosition idx2, SetPosition idx3 = 0)
{
  return static_cast<T>(idx1 * multFac1 + idx2 * multFac2 + idx3 * multFac3);
}

template <typename T, typename S>
void constructAndTestCartesianMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on ProductSet with stride " << stride);

  SLIC_INFO("\nCreating set");
  using MapType = BivariateMapType<T, S>;
  using SubMapType = typename MapType::SubMapType;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("\nCreating " << axom::slam::util::TypeToString<T>::to_string()
                          << " map on the set ");

  MapType m(&s, (T)0, stride);

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements in the map.");

  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<T>(idx1, idx2, i);
      }

  SLIC_INFO("\nChecking the elements with findValue().");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
      }

  SLIC_INFO("\nChecking the elements with SubMap.");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    SubMapType sm = m(idx1);
    for(auto idx2 = 0; idx2 < sm.size(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T v = sm.value(idx2, i);
        EXPECT_EQ(v, getVal<T>(idx1, idx2, i));
        EXPECT_EQ(sm.index(idx2), idx2);
      }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, construct_int_map)
{
  constructAndTestCartesianMap<int, RuntimeStrideType>(1);
  constructAndTestCartesianMap<int, RuntimeStrideType>(2);
  constructAndTestCartesianMap<int, RuntimeStrideType>(3);

  constructAndTestCartesianMap<int, StrideOneType>(1);

  constructAndTestCartesianMap<int, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<int, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<int, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_map)
{
  constructAndTestCartesianMap<double, StrideOneType>(1);

  constructAndTestCartesianMap<double, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<double, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<double, CompileTimeStrideType<3>>(3);

  constructAndTestCartesianMap<double, RuntimeStrideType>(1);
  constructAndTestCartesianMap<double, RuntimeStrideType>(2);
  constructAndTestCartesianMap<double, RuntimeStrideType>(3);
}

template <typename T, typename S>
void constructAndTestRelationSetMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on RelationSet with stride " << stride);

  SLIC_INFO("\nCreating set");
  using MapType = BivariateMapType<T, S>;
  using SubMapType = typename MapType::SubMapType;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);

  RelationType rel(&s1, &s2);

  std::vector<SetPosition> begin_vec(MAX_SET_SIZE1 + 1, 0);
  std::vector<SetPosition> indice_vec;

  auto curIdx = SetPosition();

  for(auto i = 0; i < MAX_SET_SIZE1; ++i)
  {
    begin_vec[i] = curIdx;
    if(MAX_SET_SIZE1 / 4 <= i && i <= MAX_SET_SIZE1 / 4 * 3)
    {
      for(auto j = MAX_SET_SIZE2 / 4; j < MAX_SET_SIZE2 / 4 * 3; ++j)
      {
        indice_vec.push_back(j);
        ++curIdx;
      }
    }
  }
  begin_vec[MAX_SET_SIZE1] = curIdx;

  rel.bindBeginOffsets(MAX_SET_SIZE1, &begin_vec);
  rel.bindIndices(indice_vec.size(), &indice_vec);
  RelationSetType s(&rel);

  const SetPosition indice_size = indice_vec.size();
  EXPECT_EQ(indice_size, s.totalSize());
  EXPECT_TRUE(s.isValid(true));

  SLIC_INFO("\nCreating " << axom::slam::util::TypeToString<T>::to_string()
                          << " map on the set ");

  MapType m(&s, (T)0, stride);

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(indice_size, m.totalSize());
  EXPECT_EQ(rel.fromSetSize(), m.firstSetSize());
  EXPECT_EQ(rel.toSetSize(), m.secondSetSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements in the map.");

  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    for(auto si = 0; si < relsubset.size(); ++si)
    {
      auto idx2 = relsubset[si];
      for(auto i = 0; i < stride; i++)
      {
        T* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<T>(idx1, idx2, i);
      }
    }
  }

  SLIC_INFO("\nChecking the elements with findValue().");
  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    auto rel_idx = 0;
    for(auto idx2 = 0; idx2 < rel.toSetSize(); ++idx2)
    {
      bool isInRel = relsubset.size() > rel_idx && relsubset[rel_idx] == idx2;
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        if(isInRel)
        {
          EXPECT_NE(ptr, nullptr);
          EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
        }
        else
        {
          EXPECT_EQ(ptr, nullptr);
        }
      }
      if(isInRel) rel_idx++;
    }
  }

  SLIC_INFO("\nChecking the elements with SubMap.");
  for(auto idx1 = 0; idx1 < rel.fromSetSize(); idx1++)
  {
    auto relsubset = rel[idx1];
    SubMapType sm = m(idx1);
    for(auto idx2 = 0; idx2 < sm.size(); ++idx2)
    {
      ASSERT_EQ(relsubset[idx2], sm.index(idx2));
      for(auto i = 0; i < stride; i++)
      {
        T v = sm.value(idx2, i);
        EXPECT_EQ(v, getVal<T>(idx1, sm.index(idx2), i));
      }
    }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, construct_int_relset_map)
{
  constructAndTestRelationSetMap<int, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<int, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<int, RuntimeStrideType>(3);

  constructAndTestRelationSetMap<int, StrideOneType>(1);

  constructAndTestRelationSetMap<int, CompileTimeStrideType<1>>(1);
  constructAndTestRelationSetMap<int, CompileTimeStrideType<2>>(2);
  constructAndTestRelationSetMap<int, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_relset_map)
{
  constructAndTestRelationSetMap<double, StrideOneType>(1);

  constructAndTestRelationSetMap<double, CompileTimeStrideType<1>>(1);
  constructAndTestRelationSetMap<double, CompileTimeStrideType<2>>(2);
  constructAndTestRelationSetMap<double, CompileTimeStrideType<3>>(3);

  constructAndTestRelationSetMap<double, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<double, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<double, RuntimeStrideType>(3);
}

template <typename StridePolicy>
void constructAndTestBivariateMapIterator(int stride)
{
  SLIC_INFO("\nCreating set");
  using DataType = double;
  using MapType = BivariateMapType<DataType, StridePolicy>;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("\nCreating " << axom::slam::util::TypeToString<DataType>::to_string()
                          << " map on the set ");
  MapType m(&s, 0.0, stride);
  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("\nSetting the elements in the map.");
  //currently can't set value using iterator
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        DataType* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<DataType>(idx1, idx2, i);
      }

  SLIC_INFO("\nChecking the elements with SubMap iterator.");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
  {
    int idx2 = 0;
    auto begin_iter = m.begin(idx1);
    for(auto iter = m.begin(idx1); iter != m.end(idx1); ++iter, ++idx2)
    {
      EXPECT_EQ(begin_iter[idx2], getVal<DataType>(idx1, idx2));
      EXPECT_EQ(*iter, getVal<DataType>(idx1, idx2));
      for(auto i = 0; i < iter.numComp(); i++)
      {
        EXPECT_EQ(iter(i), getVal<DataType>(idx1, idx2, i));
      }
    }
  }

  SLIC_INFO("\nChecking the elements with BivariateMap iterator.");
  {
    auto iter = m.begin();
    auto begin_iter = m.begin();
    auto end_iter = m.end();
    auto inval_iter = end_iter + 1;
    auto flat_idx = 0;
    for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    {
      for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      {
        EXPECT_EQ(*iter, getVal<DataType>(idx1, idx2, 0));
        EXPECT_EQ(iter.firstIndex(), idx1);
        EXPECT_EQ(iter.secondIndex(), idx2);
        for(auto i = 0; i < stride; i++)
        {
          DataType val = getVal<DataType>(idx1, idx2, i);
          EXPECT_EQ(iter.value(i), val);
          EXPECT_EQ(iter(i), val);
          EXPECT_EQ((begin_iter + flat_idx).value(i), val);
          EXPECT_EQ((end_iter - (m.totalSize() - flat_idx)).value(i), val);
          EXPECT_EQ((inval_iter - (m.totalSize() - flat_idx + 1)).value(i), val);
        }
        iter++;
        flat_idx++;
      }
    }
  }

  EXPECT_TRUE(m.isValid());
}

TEST(slam_bivariate_map, iterate)
{
  constructAndTestBivariateMapIterator<RuntimeStrideType>(1);
  constructAndTestBivariateMapIterator<RuntimeStrideType>(2);
  constructAndTestBivariateMapIterator<RuntimeStrideType>(3);

  constructAndTestBivariateMapIterator<CompileTimeStrideType<1>>(1);
  constructAndTestBivariateMapIterator<CompileTimeStrideType<2>>(2);
  constructAndTestBivariateMapIterator<CompileTimeStrideType<3>>(3);

  constructAndTestBivariateMapIterator<StrideOneType>(1);
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
