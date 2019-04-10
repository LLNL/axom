// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file slam_relation_DynamicVariable.cpp
 *
 * \brief Unit tests for Slam's DynamicVariableRelation class
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"

#include "axom/config.hpp"        // for AXOM_USE_CXX11

#include "axom/slic/interface/slic.hpp"

#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Relation.hpp"
#include "axom/slam/DynamicVariableRelation.hpp"


namespace
{
using axom::slam::RangeSet;
using axom::slam::DynamicVariableRelation;

typedef RangeSet::PositionType PositionType;
typedef RangeSet::ElementType ElementType;
// typedef axom::slam::Set::SetPosition PositionType;

const PositionType FROMSET_SIZE = 5;
const PositionType TOSET_SIZE = 8;

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::cout << "\n** " << msg << "\n\t";
  std::cout << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(),
            std::ostream_iterator<PositionType>(std::cout, " "));
}

void generateIncrementingRelations(DynamicVariableRelation* rel)
{
  PositionType curIdx = PositionType();

  for(PositionType i = 0 ; i < FROMSET_SIZE ; ++i)
  {
    for(PositionType j = 0 ; j <= i ; ++j)
    {
      rel->insert(i, j % TOSET_SIZE );
      ++curIdx;
    }
  }
}
}

TEST(slam_relation_dynamic_variable,construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be true.");

  DynamicVariableRelation emptyRel;

  EXPECT_TRUE(emptyRel.isValid(true));
}

TEST(slam_relation_dynamic_variable,construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  DynamicVariableRelation emptyRel(&fromSet, &toSet);

  EXPECT_TRUE(emptyRel.isValid(true));
}



TEST(slam_relation_dynamic_variable,construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  DynamicVariableRelation incrementingRel(&fromSet, &toSet);
  generateIncrementingRelations(&incrementingRel);

  const PositionType sz = static_cast<PositionType>(fromSet.size());
  for(PositionType idx = 0 ; idx< sz ; ++idx)
  {
    std::stringstream sstr;
    sstr << "Related to index " << idx;
    printVector(sstr.str(), incrementingRel[idx]);
  }

  EXPECT_TRUE(incrementingRel.isValid(true));
}


/// Tests for data access

TEST(slam_relation_dynamic_variable,iterate_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  // Construct the relation
  RangeSet fromSet(FROMSET_SIZE);
  RangeSet toSet(TOSET_SIZE);

  DynamicVariableRelation incrementingRel(&fromSet, &toSet);

  generateIncrementingRelations(&incrementingRel);
  EXPECT_TRUE(incrementingRel.isValid());



  // Test several data access patterns

  SLIC_INFO(".. access via double subscript.");
  {
    for(PositionType fromPos = 0 ; fromPos < fromSet.size() ; ++fromPos)
    {
      const PositionType fromSize = incrementingRel.size(fromPos);
      EXPECT_EQ(fromPos + 1, fromSize);

      for(PositionType toPos = 0 ; toPos < fromSize ; ++toPos)
      {
        PositionType actualVal = incrementingRel[fromPos][toPos]; // double
                                                                  // subscript

        PositionType expectedVal =  toPos % TOSET_SIZE;
        EXPECT_EQ( expectedVal, actualVal);
      }
    }
  }


  SLIC_INFO(".. access via delayed double subscript.");
  {
    typedef DynamicVariableRelation::RelationVec RelSet;
    for(PositionType fromPos = 0 ; fromPos < fromSet.size() ; ++fromPos)
    {
      RelSet rSet = incrementingRel[fromPos];   // first subscript

      const PositionType rSize = rSet.size();
      EXPECT_EQ(fromPos + 1, rSize);

      for(PositionType toPos = 0 ; toPos < rSize ; ++toPos)
      {
        PositionType actualVal = rSet[toPos];   // second subscript

        PositionType expectedVal =  toPos % TOSET_SIZE;
        EXPECT_EQ( expectedVal, actualVal);
      }
    }
  }

  SLIC_INFO(".. access via iterators.");
#ifdef AXOM_USE_CXX11
  {
    typedef RangeSet::iterator SetIter;
    typedef DynamicVariableRelation::RelationVecConstIterator RelSetConstIter;

    SLIC_INFO("\t using iterator begin()/end() functions");
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end() ;
        sIt != sItEnd ; ++sIt)
    {
      PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      // Test cardinality of inner relation
      PositionType actualSize = incrementingRel.size( *sIt);
      PositionType expectedSize = fromSetEltNum + 1;
      EXPECT_EQ( expectedSize, actualSize );

      RelSetConstIter toSetBegin = incrementingRel.begin(*sIt);
      RelSetConstIter toSetEnd = incrementingRel.end(*sIt);
      for(RelSetConstIter innerIt = toSetBegin ;
          innerIt != toSetEnd ; ++innerIt)
      {
        PositionType eltNum = std::distance(toSetBegin, innerIt);

        PositionType expectedVal =  (eltNum ) % TOSET_SIZE;
        PositionType actualVal = *innerIt;
        ASSERT_EQ( expectedVal, actualVal);
      }
    }

    SLIC_INFO("\t  using iterator range() function");
    typedef DynamicVariableRelation::
      RelationVecConstIteratorPair RelSetConstIterPair;
    for(SetIter sIt = fromSet.begin(), sItEnd = fromSet.end() ;
        sIt != sItEnd ; ++sIt)
    {
      // PositionType fromSetEltNum = std::distance(fromSet.begin(), sIt);

      RelSetConstIterPair toSetItPair = incrementingRel.range(*sIt);
      for(RelSetConstIter it = toSetItPair.first ;
          it < toSetItPair.second ; ++it)
      {
        PositionType toSetEltNum = std::distance(toSetItPair.first, it);
        PositionType expectedVal =  toSetEltNum % TOSET_SIZE;
        ASSERT_EQ( expectedVal, *it)
          << "incrementing relation's value was incorrect";
      }
    }
  }
#else
  SLIC_INFO("Skipping iterator tests when not using C++11");
#endif
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  UnitTestLogger logger;
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info);

  int result = RUN_ALL_TESTS();

  return result;
}
