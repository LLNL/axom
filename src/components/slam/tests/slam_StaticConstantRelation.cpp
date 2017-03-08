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


#include "common/config.hpp"  // for AXOM_USE_BOOST

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/StaticConstantRelation.hpp"

using axom::slam::RangeSet;
using axom::slam::StaticConstantRelation;

typedef RangeSet::ElementType   ElementType;
typedef RangeSet::PositionType  PositionType;

const PositionType FROMSET_SIZE = 5;
const PositionType TOSET_SIZE = 6;


TEST(gtest_slam_static_constant_relation,empty_relation)
{
  SLIC_INFO("Testing empty relation.  isValid() should be true.");

  StaticConstantRelation<> emptyRel;

  EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not valid";

  SLIC_INFO("done.");
}

TEST(gtest_slam_static_constant_relation,empty_relation_out_of_bounds)
{
  SLIC_INFO("Testing access on empty relation -- code is expected to assert and die.");

  StaticConstantRelation<> emptyRel;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif

}

TEST(gtest_slam_static_constant_relation,test_uninitialized_relation)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be TRUE since stride is 0 by default.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> emptyRel(&fromSet, &toSet);

  EXPECT_TRUE(emptyRel.isValid(true)) << "Constant relation with stride 0 should be valid";

  SLIC_INFO("done.");
}

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<ElementType>(sstr, " "));

  SLIC_INFO(msg << ":" << sstr.str());
}

template<typename VecType>
void generateIncrementingRelations(PositionType stride, VecType* offsets)
{
  VecType& offsetsVec = *offsets;

  PositionType curIdx = PositionType();

  for(PositionType i = 0; i < FROMSET_SIZE; ++i)
  {
    for(PositionType j = 0; j < stride; ++j)
    {
      offsetsVec.push_back( (i + j) % TOSET_SIZE );
      ++curIdx;
    }
  }
}

TEST(gtest_slam_static_constant_relation,simple_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);

  typedef StaticConstantRelation<>::RelationVec IndexVec;
  IndexVec offsets;

  printVector("offsets vector", offsets);

  PositionType const ELEM_STRIDE = 5;

  generateIncrementingRelations(ELEM_STRIDE, &offsets);

  printVector("offsets vector", offsets);


  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);


  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

#ifdef AXOM_USE_BOOST
  typedef RangeSet::iterator                                  SetIter;
  typedef StaticConstantRelation<>::RelationVecConstIterator  RelSetConstIter;

  SLIC_INFO("Looking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    std::stringstream sstr;
    sstr << "Inspecting element " << *sIt << " of first set.";

    PositionType actualSize = incrementingRel.size( *sIt);
    PositionType expectedSize = ELEM_STRIDE;

    sstr << "\n\t\tExpected: " << expectedSize;
    sstr << "\n\t\tActual: " <<  actualSize << "\n";

    EXPECT_EQ( expectedSize, actualSize ) << "relation for this element was incorrect size.";

    PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

    RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
    RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
    for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
    {
      PositionType toSetEltNum = std::distance(toSetBegin, innerIt);

      sstr << "\n\t\t " << toSetEltNum << ": " << *innerIt;

      PositionType expectedVal =  (fromSetEltNum + toSetEltNum) % TOSET_SIZE;
      PositionType actualVal = *innerIt;
      ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }
    SLIC_INFO( sstr.str() );
  }
#endif
  SLIC_INFO("done.");
}

TEST(gtest_slam_static_constant_relation,initialized_rel_out_of_bounds)
{
  SLIC_INFO("Testing out of bounds access on initialized relation.  Code is expected to assert and die.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);

  typedef StaticConstantRelation<>::RelationVec IndexVec;
  IndexVec offsets;

  PositionType const ELEM_STRIDE = 5;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif

}


TEST(gtest_slam_static_constant_relation,test_iterator_range)
{
  SLIC_INFO("Testing range function on incrementing relation.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);

  typedef StaticConstantRelation<>::RelationVec IndexVec;
  IndexVec offsets;

  PositionType const ELEM_STRIDE = 5;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

#ifdef AXOM_USE_BOOST
  typedef RangeSet::iterator                                      SetIter;
  typedef StaticConstantRelation<>::RelationVecConstIterator      RelSetConstIter;
  typedef StaticConstantRelation<>::RelationVecConstIteratorPair  RelSetConstIterPair;

  SLIC_INFO("Looking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    std::stringstream sstr;
    sstr << "\n\tInspecting element " << *sIt << " of first set.";
    PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

    RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
    for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
    {
      PositionType toSetEltNum = std::distance(toSetItPair.first, it);

      sstr << "\n\t\t " << toSetEltNum << ": " << *it;

      PositionType expectedVal =  (fromSetEltNum + toSetEltNum) % TOSET_SIZE;
      ASSERT_EQ( expectedVal, *it) << "incrementing relation's value was incorrect";
    }
    SLIC_INFO( sstr.str() );
  }
#endif
  SLIC_INFO("done." );
}


TEST(gtest_slam_static_constant_relation,double_subscript_test)
{
  SLIC_INFO("Testing access via double subscript.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);

  typedef StaticConstantRelation<>::RelationVec IndexVec;
  IndexVec offsets;

  PositionType const ELEM_STRIDE = 5;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

#ifdef AXOM_USE_BOOST
  typedef RangeSet::iterator SetIter;

  SLIC_INFO("Looking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    SLIC_INFO("\tInspecting element " << *sIt << " of first set.");

    for(PositionType idx = 0; idx< incrementingRel.size(*sIt); ++idx)
    {
      PositionType expectedVal =  (*sIt + idx) % TOSET_SIZE;
      PositionType actualVal = incrementingRel[*sIt][idx];
      EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }
  }
#endif

  SLIC_INFO("done.");
}


TEST(gtest_slam_static_constant_relation,delayed_double_subscript_test)
{
  SLIC_INFO("Testing access via delayed double subscript.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticConstantRelation<> incrementingRel(&fromSet, &toSet);

  typedef StaticConstantRelation<>::RelationVec IndexVec;
  IndexVec offsets;

  PositionType const ELEM_STRIDE = 5;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);
  incrementingRel.bindRelationData(offsets, ELEM_STRIDE);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

  typedef StaticConstantRelation<>::RelationSet RelSet;
  SLIC_INFO("Looking at relation's stored values...");
  for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
  {
    SLIC_INFO("\tInspecting element " << fromSet[fromPos] << " of first set (in position " << fromPos << ").");

    RelSet rSet = incrementingRel[fromPos];
    for(PositionType toPos = 0; toPos < rSet.size(); ++toPos)
    {
      PositionType expectedVal =  (fromPos + toPos) % TOSET_SIZE;
      PositionType actualVal = rSet[toPos];
      EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }
  }

  SLIC_INFO("done.");
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
