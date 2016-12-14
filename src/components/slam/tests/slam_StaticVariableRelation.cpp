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
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>
#include <sstream>      // for ATK_USE_BOOST

#include "gtest/gtest.h"

#include "common/config.hpp"

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/StaticVariableRelation.hpp"

using asctoolkit::slam::RangeSet;
using asctoolkit::slam::StaticVariableRelation;

typedef RangeSet::ElementType   ElementType;
typedef RangeSet::PositionType  PositionType;

const PositionType FROMSET_SIZE = 5;
const PositionType TOSET_SIZE = 8;


TEST(gtest_slam_static_variable_relation,empty_relation)
{
  SLIC_INFO("Testing empty relation.  isValid() should be true.");

  StaticVariableRelation emptyRel;

  EXPECT_TRUE(emptyRel.isValid(true)) << "Empty relation was not valid";

  SLIC_INFO("done.");
}



TEST(gtest_slam_static_variable_relation,empty_relation_out_of_bounds)
{

  StaticVariableRelation emptyRel;

#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH( emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

}



TEST(gtest_slam_static_variable_relation,test_uninitialized_relation)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation emptyRel(&fromSet, &toSet);

  EXPECT_FALSE(emptyRel.isValid(true)) << "Empty relation was not initialized";

  SLIC_INFO("done.");
}


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


TEST(gtest_slam_static_variable_relation,simple_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  typedef StaticVariableRelation::RelationVec IndexVec;
  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;

  printVector("begins vector",  begins);
  printVector("offsets vector", offsets);

  generateIncrementingRelations(&begins, &offsets);

  printVector("begins vector",  begins);
  printVector("offsets vector", offsets);


  incrementingRel.bindRelationData(begins, offsets);


  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

#ifdef ATK_USE_BOOST
  typedef RangeSet::iterator                                SetIter;
  typedef StaticVariableRelation::RelationVecConstIterator  RelSetConstIter;

  SLIC_INFO("Looking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    std::stringstream sstr;
    sstr << "\tInspecting element " << *sIt << " of first set.";

    PositionType actualSize = incrementingRel.size( *sIt);
    PositionType expectedSize = std::distance(fromSet.begin(), sIt) + 1;

    sstr << "\n\t\tExpected: " << expectedSize;
    sstr << "\n\t\tActual: " <<  actualSize << "\n";

    EXPECT_EQ( expectedSize, actualSize ) << "relation for this element was incorrect size.";

    RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
    RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
    for(RelSetConstIter innerIt = toSetBegin; innerIt != toSetEnd; ++innerIt)
    {
      PositionType eltNum = std::distance(toSetBegin, innerIt);

      sstr << "\t\t " << eltNum << ": " << *innerIt << "\n";

      PositionType expectedVal =  (eltNum ) % TOSET_SIZE;
      PositionType actualVal = *innerIt;
      ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }

    SLIC_INFO(sstr.str());
  }
#endif

  SLIC_INFO("done.");
}

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

#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH( incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

}

TEST(gtest_slam_static_variable_relation,test_iterator_range)
{
  SLIC_INFO("Testing range function on incrementing relation.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  typedef StaticVariableRelation::RelationVec IndexVec;
  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;
  generateIncrementingRelations(&begins, &offsets);
  incrementingRel.bindRelationData(begins, offsets);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";


#ifdef ATK_USE_BOOST
  typedef RangeSet::iterator                                    SetIter;
  typedef StaticVariableRelation::RelationVecConstIterator      RelSetConstIter;
  typedef StaticVariableRelation::RelationVecConstIteratorPair  RelSetConstIterPair;

  SLIC_INFO("\tLooking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    std::stringstream sstr;

    sstr << "\n\tInspecting element " << *sIt << " of first set.";

    RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
    for(RelSetConstIter it = toSetItPair.first; it < toSetItPair.second; ++it)
    {
      PositionType eltNum = std::distance(toSetItPair.first, it);

      sstr << "\t\t " << eltNum << ": " << *it << "\n";

      PositionType expectedVal =  (eltNum ) % TOSET_SIZE;
      PositionType actualVal = *it;
      ASSERT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }

    SLIC_INFO(sstr.str());
  }
#endif

  SLIC_INFO("done.");
}


TEST(gtest_slam_static_variable_relation,double_subscript_test)
{
  SLIC_INFO("Testing access via double subscript.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  typedef StaticVariableRelation::RelationVec IndexVec;
  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;
  generateIncrementingRelations(&begins, &offsets);
  incrementingRel.bindRelationData(begins, offsets);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

#ifdef ATK_USE_BOOST
  typedef RangeSet::iterator SetIter;

  SLIC_INFO("\tLooking at relation's stored values...");
  for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end(); sIt != sItEnd; ++sIt)
  {
    SLIC_INFO("\tInspecting element " << *sIt << " of first set.");

    for(PositionType idx = 0; idx< incrementingRel.size(*sIt  ); ++idx)
    {
      PositionType actualVal = incrementingRel[*sIt][idx];
      PositionType expectedVal = idx;
      EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }
  }
#endif
  SLIC_INFO("done.");
}


TEST(gtest_slam_static_variable_relation,delayed_double_subscript_test)
{
  SLIC_INFO("Testing access via delayed double subscript.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  StaticVariableRelation incrementingRel(&fromSet, &toSet);

  typedef StaticVariableRelation::RelationVec IndexVec;
  IndexVec begins(FROMSET_SIZE + 1);
  IndexVec offsets;
  generateIncrementingRelations(&begins, &offsets);
  incrementingRel.bindRelationData(begins, offsets);

  EXPECT_TRUE(incrementingRel.isValid(true)) << "Incrementing relation was not valid";

  typedef StaticVariableRelation::RelationSet RelSet;

  SLIC_INFO("\tLooking at relation's stored values...");
  for(PositionType fromPos = 0; fromPos < fromSet.size(); ++fromPos)
  {
    SLIC_INFO("\tInspecting element " << fromSet[fromPos] << " of first set (in position " << fromPos << ").");

    RelSet rSet = incrementingRel[fromPos];
    for(PositionType toPos = 0; toPos < rSet.size(); ++toPos)
    {
      PositionType expectedVal =  toPos;
      PositionType actualVal = rSet[toPos];
      EXPECT_EQ( expectedVal, actualVal) << "incrementing relation's value was incorrect";
    }
  }

  SLIC_INFO("done." );
}
