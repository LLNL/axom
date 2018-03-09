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
 * \file slam_AccessingRelationDataInMap.cpp
 *
 * \brief This file tests Sets, Relations and Maps working together.
 *
 * \note There are no maps.
 */


#include <iostream>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "slam/ModularInt.hpp"
#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/IndirectionPolicies.hpp"
#include "slam/CardinalityPolicies.hpp"

#include "slam/StaticRelation.hpp"

#include "slam/Map.hpp"


namespace
{

namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using slam::RangeSet;
using slam::Relation;

typedef RangeSet::ElementType ElementType;
typedef RangeSet::PositionType PositionType;
typedef PositionType SetPosition;
typedef std::vector<SetPosition>  IndexVec;

const PositionType FROMSET_SIZE = 10;
const PositionType TOSET_SIZE = 8;

typedef policies::STLVectorIndirection<PositionType, PositionType>
  STLIndirection;
typedef policies::ArrayIndirection<PositionType, PositionType>
  ArrayIndirection;

typedef policies::VariableCardinality<PositionType, STLIndirection>
  VariableCardinality;

typedef slam::StaticRelation<VariableCardinality, STLIndirection,
                             slam::RangeSet,slam::RangeSet>
  StaticVariableRelationType;


// Use a slam::ModularInt type for more interesting test data
typedef policies::CompileTimeSize<PositionType, TOSET_SIZE >  CTSize;
typedef slam::ModularInt< CTSize >                            FixedModularInt;


PositionType elementCardinality(PositionType fromPos)
{
  return fromPos;
}

PositionType relationData(PositionType fromPos, PositionType toPos)
{
  return FixedModularInt(fromPos + toPos);
}


template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "\n** " << msg << "\n\t";
  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(),
            std::ostream_iterator<PositionType>(sstr, " "));

  SLIC_INFO( sstr.str() );
}


template<typename VecType>
void generateIncrementingRelations(VecType* begins, VecType* offsets)
{
  VecType& beginsVec = *begins;
  VecType& offsetsVec = *offsets;

  PositionType curIdx = PositionType();

  for(PositionType i = 0 ; i < FROMSET_SIZE ; ++i)
  {
    beginsVec.push_back( curIdx );
    for(PositionType j = 0 ; j < elementCardinality(i) ; ++j)
    {
      offsetsVec.push_back( relationData(i,j) );
      ++curIdx;
    }
  }
  beginsVec.push_back ( curIdx );
}

}

TEST(slam_set_relation_map,access_pattern)
{
  SLIC_INFO("Testing accessing relation data.");

  IndexVec relOffsets, relIndices;
  generateIncrementingRelations(&relOffsets, &relIndices);

  RangeSet fromSet(FROMSET_SIZE), toSet(TOSET_SIZE);
  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);


  // Note: Nothing requires the relations elements to be unique
  //  -- the relation can still be valid with duplicates
  EXPECT_TRUE(incrementingRel.isValid(true));

  SLIC_INFO("-- Looking at relation's stored values...");
  for(SetPosition fromPos = SetPosition() ; fromPos < fromSet.size() ;
      ++fromPos)
  {
    SLIC_INFO(
      "--Inspecting element "
      << fromSet[fromPos]
      << " in position " << fromPos << " of first set.");

    for(SetPosition idx = 0 ; idx< incrementingRel.size( fromPos ) ; ++idx)
    {
      SetPosition posInToSet_actual = incrementingRel[fromPos][idx];
      SetPosition posInToSet_expected = relationData(fromPos,idx);
      EXPECT_EQ( posInToSet_expected, posInToSet_actual);

      SLIC_INFO(
        "-- \t pos: "
        << idx
        << " ToSet position: " << incrementingRel[fromPos][idx]
        << " ToSet element " << toSet[ incrementingRel[fromPos][idx] ] );
      ;
    }
  }
  SLIC_INFO("done.");
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // create & initialize test logger. finalized when exiting main scope
  UnitTestLogger logger;

  // axom::slic::setLoggingMsgLevel( axom::slic::message::Debug);

  result = RUN_ALL_TESTS();

  return result;
}
