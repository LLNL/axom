// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
#include "axom/slam.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;
namespace traits = axom::slam::traits;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;
using SetType = slam::RangeSet<SetPosition, SetElement>;

using StrideOneType = policies::StrideOne<SetPosition>;

template <unsigned int S>
using CompileTimeStrideType = policies::CompileTimeStride<SetPosition, S>;

using RuntimeStrideType = policies::RuntimeStride<SetPosition>;

template <typename T>
using STLIndirection = policies::STLVectorIndirection<SetPosition, T>;
using VariableCardinality =
  policies::VariableCardinality<SetPosition, STLIndirection<SetElement>>;

using RelationType = slam::StaticRelation<SetPosition,
                                          SetElement,
                                          VariableCardinality,
                                          STLIndirection<SetElement>,
                                          SetType,
                                          SetType>;

using BivariateSetType = slam::BivariateSet<>;
using ProductSetType = slam::ProductSet<>;
using RelationSetType = slam::RelationSet<RelationType>;

template <typename T, typename B, typename I, typename S>
using BivariateMapType = slam::BivariateMap<T, B, I, S>;

static const SetPosition MAX_SET_SIZE1 = 10;
static const SetPosition MAX_SET_SIZE2 = 15;

static double const multFac3 = 0000.1;
static double const multFac1 = 1000.0;
static double const multFac2 = 0010.0;

}  // end anonymous namespace

TEST(slam_bivariate_map, construct_empty_map)
{
  slam::BivariateMap<int> m;

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

template <typename T, typename B, typename I, typename S>
void constructAndTestCartesianMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on ProductSet with stride " << stride);

  SLIC_INFO("Creating set");
  using BMapType = BivariateMapType<T, B, I, S>;
  using SubMapType = typename BMapType::SubMapType;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string()
                        << " map on the set ");

  BMapType m(&s, static_cast<T>(0), stride);

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");

  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<T>(idx1, idx2, i);
      }

  SLIC_INFO("Checking the elements with findValue().");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
      }

  SLIC_INFO("Checking the elements with SubMap.");
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
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<int>;

  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestCartesianMap<int, BSet, IndPol, RuntimeStrideType>(3);

  constructAndTestCartesianMap<int, BSet, IndPol, StrideOneType>(1);

  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<int, BSet, IndPol, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  constructAndTestCartesianMap<double, BSet, IndPol, StrideOneType>(1);

  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestCartesianMap<double, BSet, IndPol, CompileTimeStrideType<3>>(3);

  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestCartesianMap<double, BSet, IndPol, RuntimeStrideType>(3);
}

template <typename T, typename B, typename I, typename S>
void constructAndTestRelationSetMap(int stride)
{
  SLIC_INFO("Testing BivariateMap on RelationSet with stride " << stride);

  SLIC_INFO("Creating set");
  using MapType = BivariateMapType<T, B, I, S>;
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

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string()
                        << " map on the set ");

  MapType m(&s, (T)0, stride);

  EXPECT_TRUE(m.isValid(true));
  EXPECT_EQ(indice_size, m.totalSize());
  EXPECT_EQ(rel.fromSetSize(), m.firstSetSize());
  EXPECT_EQ(rel.toSetSize(), m.secondSetSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");

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

  SLIC_INFO("Checking the elements with findValue().");
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

  SLIC_INFO("Checking the elements with SubMap.");
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
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<int>;

  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<int, BSet, IndPol, RuntimeStrideType>(3);

  constructAndTestRelationSetMap<int, BSet, IndPol, StrideOneType>(1);

  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<1>>(1);
  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<2>>(2);
  constructAndTestRelationSetMap<int, BSet, IndPol, CompileTimeStrideType<3>>(3);
}

TEST(slam_bivariate_map, construct_double_relset_map)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  constructAndTestRelationSetMap<double, BSet, IndPol, StrideOneType>(1);

  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<1>>(
    1);
  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<2>>(
    2);
  constructAndTestRelationSetMap<double, BSet, IndPol, CompileTimeStrideType<3>>(
    3);

  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(1);
  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(2);
  constructAndTestRelationSetMap<double, BSet, IndPol, RuntimeStrideType>(3);
}

template <typename BSet, typename StridePolicy>
void constructAndTestBivariateMapIterator(int stride)
{
  SLIC_INFO("Creating set");
  using DataType = double;
  using IndPol = STLIndirection<DataType>;
  using MapType = BivariateMapType<DataType, BSet, IndPol, StridePolicy>;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<DataType>::to_string()
                        << " map on the set ");
  MapType m(&s, 0.0, stride);
  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Setting the elements in the map.");
  //currently can't set value using iterator
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        DataType* valPtr = m.findValue(idx1, idx2, i);
        EXPECT_NE(valPtr, nullptr);
        *valPtr = getVal<DataType>(idx1, idx2, i);
      }

  SLIC_INFO("Checking the elements with SubMap iterator.");
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

  SLIC_INFO("Checking the elements with BivariateMap iterator.");
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
  using BSet = BivariateSetType;

  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(1);
  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(2);
  constructAndTestBivariateMapIterator<BSet, RuntimeStrideType>(3);

  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<1>>(1);
  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<2>>(2);
  constructAndTestBivariateMapIterator<BSet, CompileTimeStrideType<3>>(3);

  constructAndTestBivariateMapIterator<BSet, StrideOneType>(1);
}

template <typename T, typename B, typename I, typename S>
void testScopedCopyBehavior(int stride)
{
  using BMapType = BivariateMapType<T, B, I, S>;

  SetType s1(MAX_SET_SIZE1);
  SetType s2(MAX_SET_SIZE2);
  ProductSetType s(&s1, &s2);

  EXPECT_EQ(s.size(), MAX_SET_SIZE1 * MAX_SET_SIZE2);
  EXPECT_TRUE(s.isValid());

  SLIC_INFO("Creating " << slam::util::TypeToString<T>::to_string()
                        << " map on the set ");

  BMapType m;
  {
    BMapType m_inner(&s, static_cast<T>(0), stride);

    EXPECT_TRUE(m_inner.isValid());
    EXPECT_EQ(s.size(), m_inner.totalSize());
    EXPECT_EQ(m_inner.stride(), stride);

    SLIC_INFO("Setting the elements in the map.");

    for(auto idx1 = 0; idx1 < m_inner.firstSetSize(); ++idx1)
      for(auto idx2 = 0; idx2 < m_inner.secondSetSize(); ++idx2)
        for(auto i = 0; i < stride; i++)
        {
          T* valPtr = m_inner.findValue(idx1, idx2, i);
          EXPECT_NE(valPtr, nullptr);
          *valPtr = getVal<T>(idx1, idx2, i);
        }

    m = m_inner;
  }

  EXPECT_TRUE(m.isValid());
  EXPECT_EQ(s.size(), m.totalSize());
  EXPECT_EQ(m.stride(), stride);

  SLIC_INFO("Checking the elements with findValue().");
  for(auto idx1 = 0; idx1 < m.firstSetSize(); ++idx1)
    for(auto idx2 = 0; idx2 < m.secondSetSize(); ++idx2)
      for(auto i = 0; i < stride; i++)
      {
        T* ptr = m.findValue(idx1, idx2, i);
        EXPECT_NE(ptr, nullptr);
        EXPECT_EQ(*ptr, getVal<T>(idx1, idx2, i));
      }
}

TEST(slam_bivariate_map, testScopedMapBehavior)
{
  using BSet = BivariateSetType;
  using IndPol = STLIndirection<double>;

  testScopedCopyBehavior<double, BSet, IndPol, StrideOneType>(1);

  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<1>>(1);
  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<2>>(2);
  testScopedCopyBehavior<double, BSet, IndPol, CompileTimeStrideType<3>>(3);

  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(1);
  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(2);
  testScopedCopyBehavior<double, BSet, IndPol, RuntimeStrideType>(3);
}

TEST(slam_bivariate_map, traits)
{
  EXPECT_TRUE(traits::indices_use_indirection<RelationSetType>::value);
  EXPECT_TRUE(traits::indices_use_indirection<BivariateSetType>::value);
  EXPECT_FALSE(traits::indices_use_indirection<ProductSetType>::value);
}

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
