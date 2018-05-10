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


#include <iostream>
#include <iterator>

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
const PositionType ELEM_STRIDE = 5;

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

  RelationType rel;

  // Modify the relation values
  {
    RelationType rel_t( &fromSet, &toSet);

    for(int i=0 ; i<FROMSET_SIZE ; i++)
    {
      for(int j=0 ; j<ELEM_STRIDE ; j++)
      {
        rel_t.modify(i,j,i);
      }
    }

    rel = rel_t;
  }

  EXPECT_TRUE( rel.isValid(true) );
  EXPECT_EQ( FROMSET_SIZE, rel.size() );

  // Check that accessed relation values are as expected
  for(int i=0 ; i<rel.size() ; i++)
  {
    EXPECT_EQ( ELEM_STRIDE, (int)rel[i].size() );
    for(int j=0 ; j < (int)rel[i].size() ; j++)
    {
      EXPECT_EQ( i, rel[i][j]);
    }
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
