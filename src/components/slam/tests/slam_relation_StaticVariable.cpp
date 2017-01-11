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
 * \file testStaticVariableRelation.cxx
 *
 * \brief Unit tests for Slam's StaticConstantRelation class
 */


#include <iostream>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for AXOM_USE_BOOST

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/StaticVariableRelation.hpp"


namespace {
  using axom::slam::RangeSet;
  using axom::slam::StaticVariableRelation;
  typedef StaticVariableRelation::RelationVec IndexVec;

  typedef RangeSet::ElementType               ElementType;
  typedef RangeSet::PositionType              PositionType;

  const PositionType FROMSET_SIZE = 5;
  const PositionType TOSET_SIZE = 8;


  template<typename StrType, typename VecType>
  void printVector(StrType const& msg, VecType const& vec)
  {
    std::stringstream sstr;

    sstr << "\n** " << msg << "\n\t";
    sstr << "Array of size " << vec.size() << ": ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<PositionType>(sstr, " "));

    SLIC_INFO( sstr.str() );
  }


  template<typename VecType>
  void generateIncrementingRelations(VecType* begins, VecType* offsets)
  {
    VecType& beginsVec = *begins;
    VecType& offsetsVec = *offsets;

    PositionType curIdx = PositionType();

    for(PositionType i = 0; i < FROMSET_SIZE; ++i)
    {
      beginsVec[i] = curIdx;
      for(PositionType j = 0; j <= i; ++j)
      {
        offsetsVec.push_back( j % TOSET_SIZE );
        ++curIdx;
      }
    }
    beginsVec[FROMSET_SIZE] = curIdx;
  }


}

TEST(gtest_slam_static_variable_relation,construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be true.");

  StaticVariableRelation emptyRel;

  EXPECT_TRUE(emptyRel.isValid(true));
}

TEST(gtest_slam_static_variable_relation,construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation emptyRel(&fromSet, &toSet);

  EXPECT_FALSE(emptyRel.isValid(true));
}


TEST(gtest_slam_static_variable_relation,construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;

  SLIC_INFO("Uninitialized relation data");
  printVector("begins vector",  begins);
  printVector("offsets vector", offsets);

  generateIncrementingRelations(&begins, &offsets);

  SLIC_INFO("Initialized relation data");
  printVector("begins vector",  begins);
  printVector("offsets vector", offsets);


  incrementingRel.bindRelationData(begins, offsets);


  EXPECT_TRUE(incrementingRel.isValid(true));
}


TEST(gtest_slam_static_variable_relation,iterate_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  // Initialize the relation
  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;
  generateIncrementingRelations(&begins, &offsets);
  incrementingRel.bindRelationData(begins, offsets);

  EXPECT_TRUE(incrementingRel.isValid());


  // Test several data access patterns

  SLIC_INFO(".. access via double subscript.");
  {
    for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
    {
      const PositionType fromSize = incrementingRel.size(fromPos);
      EXPECT_EQ(fromPos + 1, fromSize);

      for(PositionType toPos = 0; toPos < fromSize; ++toPos)
      {
        PositionType actualVal = incrementingRel[fromPos][toPos]; // double subscript

        PositionType expectedVal =  toPos % TOSET_SIZE;
        EXPECT_EQ( expectedVal, actualVal);
      }
    }
  }

  SLIC_INFO(".. access via delayed double subscript.");
  {
    typedef StaticVariableRelation::RelationSet RelSet;
    for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
    {
      RelSet rSet = incrementingRel[fromPos];   // first subscript

      const PositionType rSize = rSet.size();
      EXPECT_EQ(fromPos + 1, rSize);

      for(PositionType toPos = 0; toPos < rSize; ++toPos)
      {
        PositionType actualVal = rSet[toPos];   // second subscript

        PositionType expectedVal =  toPos % TOSET_SIZE;
        EXPECT_EQ( expectedVal, actualVal);
      }
    }
  }

  SLIC_INFO(".. access via iterators.");
#ifdef AXOM_USE_BOOST
  {
    typedef RangeSet::iterator                                SetIter;
    typedef StaticVariableRelation::RelationVecConstIterator  RelSetConstIter;

    SLIC_INFO("\t using iterator begin()/end() functions");
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
      PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      // Test cardinality of inner relation
      PositionType actualSize = incrementingRel.size( *sIt);
      PositionType expectedSize = fromSetEltNum + 1;
      EXPECT_EQ( expectedSize, actualSize );

      RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
      RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
      for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
      {
        PositionType toSetEltNum = std::distance(toSetBegin, innerIt);

        PositionType expectedVal =  (toSetEltNum ) % TOSET_SIZE;
        PositionType actualVal = *innerIt;
        ASSERT_EQ( expectedVal, actualVal);
      }
    }

    SLIC_INFO("\t  using iterator range() function");
    typedef StaticVariableRelation::RelationVecConstIteratorPair RelSetConstIterPair;
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
    {
      // PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
      for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
      {
        PositionType toSetEltNum = std::distance(toSetItPair.first, it);
        PositionType expectedVal =  toSetEltNum % TOSET_SIZE;
        ASSERT_EQ( expectedVal, *it) << "incrementing relation's value was incorrect";
      }
    }
  }
#else
  SLIC_WARNING("Must compile with boost enabled to test iterators")
#endif
}




TEST(gtest_slam_static_variable_relation,empty_relation_out_of_bounds)
{
  StaticVariableRelation emptyRel;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}


/// Tests for out-of-bounds access

TEST(gtest_slam_static_variable_relation,initialized_rel_out_of_bounds)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);
  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  typedef StaticVariableRelation::RelationVec IndexVec;
  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;

  generateIncrementingRelations(&begins, &offsets);
  incrementingRel.bindRelationData(begins, offsets);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
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
