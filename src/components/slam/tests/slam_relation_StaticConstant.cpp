/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file slam_StaticConstantRelation.cpp
 *
 * \brief Unit tests for Slam's StaticConstantRelation class
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for AXOM_USE_BOOST

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/StaticConstantRelation.hpp"
#include "slam/SizePolicies.hpp"
#include "slam/ModularInt.hpp"


// TODO: This currently only tests RuntimeStride relations.
//       We should also test the CompileTimeStride relations.

namespace {
  using axom::slam::RangeSet;
  using axom::slam::StaticConstantRelation;
  using axom::slam::ModularInt;

  typedef RangeSet::ElementType                 ElementType;
  typedef RangeSet::PositionType                PositionType;
  typedef StaticConstantRelation<>::RelationVec IndexVec;

  const PositionType FROMSET_SIZE = 5;
  const PositionType TOSET_SIZE = 6;
  const PositionType ELEM_STRIDE = 5;

  // Use a slam::ModularInt type for more interesting test data
  typedef axom::slam::policies::CompileTimeSizeHolder<int, ELEM_STRIDE> SizePolicy;
  typedef ModularInt< SizePolicy >                                            FixedModularInt;

  template<typename StrType, typename VecType>
  void printVector(StrType const& msg, VecType const& vec)
  {
    std::stringstream sstr;

    sstr << "Array of size " << vec.size() << ": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<ElementType>(sstr, " "));

    SLIC_INFO(msg << ": " << sstr.str());
  }

  /**
   * \brief Sets the value at relation element (i,j) to (i + j) % ELEM_SIZE using slam::ModularInt
   */
  template<typename VecType>
  void generateIncrementingRelations(PositionType stride, VecType* offsets)
  {
    VecType& offsetsVec = *offsets;

    PositionType curIdx = PositionType();

    for(PositionType i = 0; i < FROMSET_SIZE; ++i)
    {
      FixedModularInt modInt(i);
      for(PositionType j = 0; j < stride; ++j)
      {
        offsetsVec.push_back( modInt + j);
        ++curIdx;
      }
    }
  }

}

TEST(gtest_slam_relation_static_constant,construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be true.");

  StaticConstantRelation<> emptyRel;
  EXPECT_TRUE(emptyRel.isValid());
}

TEST(gtest_slam_relation_static_constant,construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be TRUE since stride is 0 by default.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> emptyRel(&fromSet, &toSet);
  EXPECT_TRUE(emptyRel.isValid(true));
}


TEST(gtest_slam_relation_static_constant,construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  IndexVec offsets;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);
  printVector("offsets vector for relation", offsets);

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

  const bool VERBOSE_OUTPUT = true;
  EXPECT_TRUE(incrementingRel.isValid(VERBOSE_OUTPUT));
}


/// Tests for data access

TEST(gtest_slam_relation_static_constant,iterate_relation)
{
  SLIC_INFO("Testing simple iteration of relation data.");

  // Construct relation
  IndexVec offsets;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

  // Test several data access patterns

  SLIC_INFO(".. access via double subscript.");
  {
    for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
    {
      const PositionType fromSize = incrementingRel.size(fromPos);
      EXPECT_EQ(ELEM_STRIDE, fromSize);

      for(PositionType toPos = 0; toPos < fromSize; ++toPos)
      {
        PositionType actualVal = incrementingRel[fromPos][toPos]; // double subscript

        PositionType expectedVal =  FixedModularInt(fromPos + toPos);
        EXPECT_EQ( expectedVal, actualVal);
      }
    }
  }

  SLIC_INFO(".. access via delayed double subscript.");
  {
    typedef StaticConstantRelation<>::RelationSet RelSet;
    for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
    {
      RelSet rSet = incrementingRel[fromPos];   // first subscript

      const PositionType rSize = rSet.size();
      EXPECT_EQ(ELEM_STRIDE, rSize);

      for(PositionType toPos = 0; toPos < rSize; ++toPos)
      {
        PositionType actualVal = rSet[toPos];   // second subscript

        PositionType expectedVal =  FixedModularInt(fromPos + toPos);
        EXPECT_EQ( expectedVal, actualVal);
      }
    }

  }

  SLIC_INFO(".. access via iterators.");
#ifdef AXOM_USE_BOOST
  {
    typedef RangeSet::iterator                                  SetIter;
    typedef StaticConstantRelation<>::RelationVecConstIterator  RelSetConstIter;

    SLIC_INFO("\t using iterator begin()/end() functions");
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
      // Test cardinality of inner relation
      PositionType actualSize = incrementingRel.size( *sIt);
      PositionType expectedSize = ELEM_STRIDE;
      EXPECT_EQ( expectedSize, actualSize );

      PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      // Get iterators to inner set of relation
      RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
      RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
      for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
      {
        PositionType actualVal = *innerIt;

        PositionType toSetEltNum = std::distance(toSetBegin, innerIt);
        PositionType expectedVal =  FixedModularInt(fromSetEltNum + toSetEltNum);
        ASSERT_EQ( expectedVal, actualVal);
      }
    }

    SLIC_INFO("\t  using iterator range() function");
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
      PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      typedef StaticConstantRelation<>::RelationVecConstIterator      RelSetConstIter;
      typedef StaticConstantRelation<>::RelationVecConstIteratorPair  RelSetConstIterPair;

      RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
      for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
      {
        PositionType toSetEltNum = std::distance(toSetItPair.first, it);
        PositionType expectedVal =  FixedModularInt(fromSetEltNum + toSetEltNum);
        ASSERT_EQ( expectedVal, *it) << "incrementing relation's value was incorrect";
      }
    }
  }
#else
  SLIC_WARNING("Must compile with boost enabled to test iterators")
#endif
}


/// Tests for out-of-bounds access


TEST(gtest_slam_relation_static_constant,out_of_bounds_empty)
{
  SLIC_INFO("Testing access on empty relation -- code is expected to assert and die.");

  StaticConstantRelation<> emptyRel;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_ASSSERT is disabled in release mode, so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH( emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif

}

TEST(gtest_slam_relation_static_constant,out_of_bounds_initialized)
{
  SLIC_INFO("Testing out of bounds access on initialized relation.  Code is expected to assert and die.");

  IndexVec offsets;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
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

  // create & initialize test logger. finalized when exiting main scope
  UnitTestLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
