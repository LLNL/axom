/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file slam_relation_DynamicConstant.cpp
 *
 * \brief Unit tests for Slam's DynamicConstantRelation class
 */


#include "gtest/gtest.h"

#include "axom/config.hpp"        // for AXOM_USE_CXX11

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/SizePolicies.hpp"
#include "slam/IndirectionPolicies.hpp"
#include "slam/CardinalityPolicies.hpp"
#include "slam/ModularInt.hpp"


#include "slam/StridePolicies.hpp"
#include "slam/DynamicConstantRelation.hpp"



namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;


typedef slam::RangeSet::ElementType ElementType;
typedef slam::RangeSet::PositionType PositionType;

const PositionType FROMSET_SIZE = 5;
const PositionType TOSET_SIZE = 6;
const PositionType ELEM_STRIDE = 6;

typedef policies::CompileTimeStride<PositionType, ELEM_STRIDE>      CTStride;
typedef policies::RuntimeStride<PositionType>                       RTStride;

typedef policies::
  ConstantCardinality<PositionType,CTStride>       ConstantCardinalityCT;
typedef policies::
  ConstantCardinality<PositionType, RTStride>       ConstantCardinalityRT;

typedef policies::
  STLVectorIndirection<PositionType, PositionType>  STLIndirection;
typedef policies::
  ArrayIndirection<PositionType, PositionType>      ArrayIndirection;

typedef std::vector<PositionType>                                   IndexVec;

typedef slam::DynamicConstantRelation<ConstantCardinalityCT>        RelationType;

} //end anonymous namespace


/** \brief Utility function to set the value at rel[i][j] to i */
template<typename DynamicSetType>
RelationType generateRelation(DynamicSetType& fromSet,
                              DynamicSetType& toSet)
{
  RelationType rel( &fromSet, &toSet);

  const int fromSize = fromSet.size();
  const int toSize = toSet.size();

  for(int i=0 ; i< fromSize ; ++i)
  {
    for(int j=0 ; j< toSize ; ++j)
    {
      rel.modify(i,j,i);
    }
  }

  return rel;
}


TEST(slam_relation_dynamic_constant,construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  RelationType emptyRel;
  EXPECT_FALSE(emptyRel.isValid());
}


TEST(slam_relation_dynamic_constant, assignment)
{
  SLIC_INFO("Testing assignment of values to relation");

  slam::DynamicSet<> fromSet(FROMSET_SIZE);
  slam::DynamicSet<> toSet(TOSET_SIZE);

  RelationType rel = generateRelation(fromSet,toSet);

  EXPECT_TRUE( rel.isValid(true) );

  const int fromSize = rel.size();
  EXPECT_EQ( FROMSET_SIZE, fromSize);

  // Check that accessed relation values are as expected
  for(int i=0 ; i< fromSize ; ++i)
  {
    EXPECT_EQ( rel[i].size(), rel.size(i));

    // Test double indirection
    {
      const int relSize = rel.size(i);
      EXPECT_EQ( ELEM_STRIDE, relSize );
      for(int j=0 ; j < relSize ; j++)
      {
        EXPECT_EQ( i, rel[i][j]);
      }
    }

    // Test RelationSet
    {
      RelationType::RelationSet set = rel[i];
      const int relSize = set.size();
      for(int j=0 ; j < relSize ; ++j)
      {
        EXPECT_EQ( i, set[j]);
      }

      EXPECT_EQ(set, rel.at(i));
    }
  }

}

#ifdef AXOM_USE_CXX11
TEST(slam_relation_dynamic_constant, iterators)
{
  SLIC_INFO("Testing iterator interface");

  // Add tests for relation iterators
  slam::DynamicSet<> fromSet(FROMSET_SIZE);
  slam::DynamicSet<> toSet(TOSET_SIZE);

  {
    RelationType rel = generateRelation(fromSet,toSet);
    const int fromSize = rel.size();

    for(int i=0 ; i< fromSize ; ++i)
    {
      ElementType val = i;
      for(RelationType::RelationIterator it = rel.begin(i), itEnd = rel.end(i) ;
          it != itEnd ;
          ++it)
      {
        EXPECT_EQ(val, *it);
      }

      for(RelationType::RelationIteratorPair itPair = rel.range(i) ;
          itPair.first != itPair.second ;
          ++itPair.first)
      {
        EXPECT_EQ(val, *itPair.first);
      }
    }
  }
}

TEST(slam_relation_dynamic_constant, const_iterators)
{
  SLIC_INFO("Testing iterator interface");

  // Add tests for relation iterators
  slam::DynamicSet<> fromSet(FROMSET_SIZE);
  slam::DynamicSet<> toSet(TOSET_SIZE);

  {
    const RelationType rel = generateRelation(fromSet,toSet);
    const int fromSize = rel.size();

    for(int i=0 ; i< fromSize ; ++i)
    {
      const ElementType val = i;
      for(RelationType::RelationConstIterator it = rel.begin(i),
          itEnd = rel.end(i) ;
          it != itEnd ;
          ++it)
      {
        EXPECT_EQ(val, *it);
      }

      for(RelationType::RelationConstIteratorPair itPair = rel.range(i) ;
          itPair.first != itPair.second ;
          ++itPair.first)
      {
        EXPECT_EQ(val, *itPair.first);
      }
    }
  }
}
#endif // AXOM_USE_CXX11

TEST(slam_relation_dynamic_constant, remove)
{
  SLIC_INFO("Testing ability to remove relations");

  slam::DynamicSet<> fromSet(FROMSET_SIZE);
  slam::DynamicSet<> toSet(TOSET_SIZE);

  RelationType rel = generateRelation(fromSet,toSet);
  EXPECT_EQ(rel.size(), rel.numberOfValidEntries());

  const int fromSize = rel.size();
  int removeCount = 0;
  for(int i=1 ; i< fromSize ; i+=2)
  {
    rel.remove(i);
    ++removeCount;
  }

  // Check counts
  EXPECT_TRUE(rel.isValid());
  EXPECT_EQ(rel.size() - removeCount, rel.numberOfValidEntries());

  // Test individual relations
  for(int i=0 ; i< fromSize ; i+=2)
  {
    EXPECT_TRUE(rel.isValidEntry(i));
  }

  for(int i=1 ; i< fromSize ; i+=2)
  {
    EXPECT_FALSE(rel.isValidEntry(i));
  }



}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
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
