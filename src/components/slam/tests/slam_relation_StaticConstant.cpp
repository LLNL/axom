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
 * \file slam_relation_StaticConstant.cpp
 *
 * \brief Unit tests for Slam's StaticRelation class
 * configured with constant per-element cardinality
 *
 * Exercises several different variants of the class.
 * Namely:
 * * with a runtime or compile time striding of offsets per element
 * * where the underlying indirection array uses STL vectors or arrays
 *
 */


#include <iostream>
#include <iterator>

#include "gtest/gtest.h"

#include "axom/config.hpp"  // for AXOM_USE_CXX11

#include "slic/slic.hpp"

#include "slam/RangeSet.hpp"
#include "slam/Relation.hpp"
#include "slam/SizePolicies.hpp"
#include "slam/IndirectionPolicies.hpp"
#include "slam/CardinalityPolicies.hpp"
#include "slam/ModularInt.hpp"


#include "slam/StridePolicies.hpp"
#include "slam/StaticRelation.hpp"

namespace
{

namespace slam = axom::slam;
namespace policies = axom::slam::policies;


typedef slam::RangeSet::ElementType ElementType;
typedef slam::RangeSet::PositionType PositionType;

const PositionType FROMSET_SIZE = 5;
const PositionType TOSET_SIZE = 6;
const PositionType ELEM_STRIDE = 5;

typedef policies::
  CompileTimeStride<PositionType, ELEM_STRIDE>      CTStride;
typedef policies::
  RuntimeStride<PositionType>                       RTStride;

typedef policies::
  ConstantCardinality<PositionType, CTStride>       ConstantCardinalityCT;
typedef policies::
  ConstantCardinality<PositionType, RTStride>       ConstantCardinalityRT;

typedef policies::
  STLVectorIndirection<PositionType, PositionType>  STLIndirection;
typedef policies::
  ArrayIndirection<PositionType, PositionType>      ArrayIndirection;

typedef std::vector<PositionType>                   IndexVec;


typedef slam::StaticRelation<ConstantCardinalityCT, STLIndirection,
                             slam::RangeSet,
                             slam::RangeSet>
  StaticConstantRelationType;


// Use a slam::ModularInt type for more interesting test data
typedef policies::CompileTimeSize<int, ELEM_STRIDE> CTSize;
typedef slam::ModularInt< CTSize >                  FixedModularInt;

template<typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(),
            std::ostream_iterator<ElementType>(sstr, " "));

  SLIC_INFO(msg << ": " << sstr.str());
}

PositionType elementCardinality(PositionType AXOM_NOT_USED(fromPos) )
{
  return ELEM_STRIDE;
}

/**
 * \brief Sets the value at relation element (i,j) to (i + j) % ELEM_SIZE using
 *  slam::ModularInt
 */
PositionType relationData(PositionType fromPos, PositionType toPos)
{
  return FixedModularInt(fromPos + toPos);
}


template<typename VecType>
void generateIncrementingRelations(PositionType stride, VecType* offsets)
{
  VecType& offsetsVec = *offsets;

  PositionType curIdx = PositionType();

  for(PositionType i = 0 ; i < FROMSET_SIZE ; ++i)
  {
    EXPECT_EQ(elementCardinality(i), stride);

    for(PositionType j = 0 ; j < elementCardinality(i) ; ++j)
    {
      offsetsVec.push_back( relationData(i,j) );
      ++curIdx;
    }
  }
}


/**
 * \brief Traverses the relation's entities using the double subscript access
 *  API
 *
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
 */
template<typename RelationType>
void traverseRelation_doubleSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using double subscript: " );
  for(PositionType fromPos = 0 ; fromPos < rel.fromSet()->size() ; ++fromPos)
  {
    const PositionType fromSize = rel.size(fromPos);
    EXPECT_EQ( elementCardinality(fromPos), fromSize );

    for(int toPos = 0 ; toPos < fromSize ; ++toPos)
    {
      PositionType actualVal = rel[fromPos][toPos];
      EXPECT_EQ(relationData(fromPos,toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using separated subscript operators
 *
 * The first subscript operator gets the set of entities in the ToSet that are
 * mapped to the given element of the relation's FromSet
 *
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
 */
template<typename RelationType>
void traverseRelation_delayedSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using delayed second subscript: " );
  for(PositionType fromPos = 0 ; fromPos < rel.fromSet()->size() ; ++fromPos)
  {
    const PositionType fromSize = rel.size(fromPos);
    EXPECT_EQ( elementCardinality(fromPos), fromSize );

    typename RelationType::RelationSubset set = rel[fromPos];
    for(int toPos = 0 ; toPos < set.size() ; ++toPos)
    {
      PositionType actualVal = rel[fromPos][toPos];
      EXPECT_EQ(relationData(fromPos,toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using the iterator API (begin()/end() )
 *
 * \note The iterator API depends on C++11
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
 */
template<typename RelationType>
void iterateRelation_begin_end(RelationType& rel)
{
#ifdef AXOM_USE_CXX11
  typedef typename RelationType::FromSetType FromSet;
  typedef typename FromSet::iterator FromSetIter;

  typedef typename RelationType::RelationIterator RelIter;

  SLIC_INFO("Traversing relation data using iterator begin()/end() functions");
  for(FromSetIter sIt = rel.fromSet()->begin(), sItEnd = rel.fromSet()->end() ;
      sIt != sItEnd ; ++sIt)
  {
    PositionType actualSize = rel.size( *sIt);

    PositionType fromSetEltNum = std::distance(rel.fromSet()->begin(), sIt);
    EXPECT_EQ( elementCardinality(fromSetEltNum), actualSize );

    RelIter toSetBegin = rel.begin(*sIt);
    RelIter toSetEnd   = rel.end(*sIt);
    for(RelIter relIt = toSetBegin ; relIt != toSetEnd ; ++relIt)
    {
      PositionType toSetEltNum = std::distance(toSetBegin, relIt);
      ASSERT_EQ( relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
#else
  AXOM_DEBUG_VAR(rel);
  SLIC_INFO("Skipping iterator tests when not using C++11");
#endif
}

/**
 * \brief Traverses relation using the iterator range API
 *
 * \note The iterator API depends on C++11
 * \note Expects cardinality and relation data for each element to match the
 *  results of the elementCardinality() and relationData() functions above,
 *  respectively.
 */
template<typename RelationType>
void iterateRelation_range(RelationType& rel)
{
#ifdef AXOM_USE_CXX11
  typedef typename RelationType::FromSetType FromSet;
  typedef typename FromSet::iterator FromSetIter;
  typedef typename FromSet::iterator_pair FromSetIterPair;

  typedef typename RelationType::RelationIterator RelIter;
  typedef typename RelationType::RelationIteratorPair RelIterPair;

  SLIC_INFO("Traversing relation data using iterator range() functions");
  FromSetIterPair itPair = rel.fromSet()->range();
  for(FromSetIter sIt = itPair.first ; sIt != itPair.second ; ++sIt)
  {
    PositionType fromSetEltNum = std::distance(itPair.first, sIt);

    RelIterPair toSetItPair = rel.range(*sIt);
    for(RelIter relIt = toSetItPair.first ;
        relIt != toSetItPair.second ; ++relIt)
    {
      PositionType toSetEltNum = std::distance(toSetItPair.first, relIt);
      ASSERT_EQ( relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
#else
  AXOM_DEBUG_VAR(rel);
  SLIC_INFO("Skipping iterator tests when not using C++11");
#endif
}

} // end anonymous namespace


TEST(slam_relation_static_constant,construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  StaticConstantRelationType emptyRel;
  EXPECT_FALSE(emptyRel.isValid());
}

TEST(slam_relation_static_constant,construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  StaticConstantRelationType emptyRel(&fromSet, &toSet);
  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_relation_static_constant,construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  printVector("Relation indices vector ", relIndices);

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  StaticConstantRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);    // init the
                                                                    // begins
                                                                    // set
  incrementingRel.bindIndices(relIndices.size(), &relIndices);  // init the
                                                                // relation
                                                                // indices set

  EXPECT_TRUE(incrementingRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("CTStride_vector_simple_construct");
  traverseRelation_doubleSubscript(incrementingRel);
  traverseRelation_delayedSubscript(incrementingRel);
  iterateRelation_begin_end(incrementingRel);
  iterateRelation_range(incrementingRel);

}


TEST(slam_relation_static_constant,construct_builder)
{
  SLIC_INFO(
    "Checking if we can instantiate a concrete StaticRelation (constant).");

  typedef StaticConstantRelationType::RelationBuilder RelationBuilder;

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  IndexVec offsets;
  generateIncrementingRelations(ELEM_STRIDE, &offsets);

  StaticConstantRelationType relation =
    RelationBuilder()
    .fromSet( &fromSet)
    .toSet( &toSet)
    .begins( RelationBuilder::BeginsSetBuilder()
             .stride(ELEM_STRIDE))
    .indices( RelationBuilder::IndicesSetBuilder()
              .size(offsets.size())
              .data(&offsets));

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("CTStride_vector_builder_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);
}


/// Tests for out-of-bounds access

TEST(slam_relation_static_constant,out_of_bounds_empty)
{
  SLIC_INFO("Testing access on empty relation "
            <<"-- code is expected to assert and die.");

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);
  StaticConstantRelationType emptyRel(&fromSet, &toSet);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif

}

TEST(slam_relation_static_constant,out_of_bounds_initialized)
{
  SLIC_INFO("Testing out of bounds access on initialized relation.  "
            << "Code is expected to assert and die.");

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);

  StaticConstantRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), ELEM_STRIDE); // init the
                                                                 // begins set
  incrementingRel.bindIndices(relIndices.size(), &relIndices);   // init the
                                                                 // relation
                                                                 // indices set

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED( incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure checks in release mode.");
#endif

}

/// Tests that vary the relation's policies (STL vs. array indirection and
// compile time vs. runtime offset striding)

TEST(slam_relation_static_constant,runtime_stride_STLIndirection)
{
  SLIC_INFO("Tests for Static Relation with runtime stride and STL Indirection");

  typedef slam::StaticRelation<
      ConstantCardinalityRT, STLIndirection,
      slam::RangeSet, slam::RangeSet >
    StaticConstantRelation_RT_STL;

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);


  // --  Construct empty and uninitialized relation
  StaticConstantRelation_RT_STL emptyRel;
  EXPECT_FALSE( emptyRel.isValid(true));

  StaticConstantRelation_RT_STL uninitRel(&fromSet, &toSet);
  EXPECT_FALSE( uninitRel.isValid(true));


  // -- Simple relation construction
  StaticConstantRelation_RT_STL relation(&fromSet, &toSet);
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  relation.bindIndices(relIndices.size(), &relIndices);

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_vector_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  typedef StaticConstantRelation_RT_STL::RelationBuilder RelationBuilder;
  StaticConstantRelation_RT_STL builderRel =
    RelationBuilder()
    .fromSet( &fromSet)
    .toSet( &toSet)
    .begins( RelationBuilder::BeginsSetBuilder()
             .stride(ELEM_STRIDE))
    .indices(RelationBuilder::IndicesSetBuilder()
             .size(relIndices.size())
             .data(&relIndices) )
  ;
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_vector_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);
}

TEST(slam_relation_static_constant,runtime_stride_ArrayIndirection)
{
  SLIC_INFO("Tests for Static Relation "
            <<" with runtime stride and array Indirection");

  typedef slam::StaticRelation<
      ConstantCardinalityRT, ArrayIndirection,
      slam::RangeSet, slam::RangeSet>
    StaticConstantRelation_RT_Array;

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  PositionType* data = &relIndices[0];    // Get a pointer to the data


  // --  Construct empty and uninitialized relation
  StaticConstantRelation_RT_Array emptyRel;
  EXPECT_FALSE( emptyRel.isValid());

  StaticConstantRelation_RT_Array uninitRel(&fromSet, &toSet);
  EXPECT_FALSE( uninitRel.isValid());


  // -- Simple relation construction
  StaticConstantRelation_RT_Array relation(&fromSet, &toSet);
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  relation.bindIndices(relIndices.size(), data);

  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  typedef StaticConstantRelation_RT_Array::RelationBuilder RelationBuilder;
  StaticConstantRelation_RT_Array builderRel =
    RelationBuilder()
    .fromSet( &fromSet)
    .toSet( &toSet)
    .begins(RelationBuilder::BeginsSetBuilder()
            .stride(ELEM_STRIDE))
    .indices(RelationBuilder::IndicesSetBuilder()
             .size(relIndices.size())
             .data(data));
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);
}


TEST(slam_relation_static_constant,compileTime_stride_ArrayIndirection)
{
  SLIC_INFO("Tests for Static Relation with "
            << " runtime stride and array Indirection");

  typedef slam::StaticRelation<
      ConstantCardinalityCT, ArrayIndirection,
      slam::RangeSet, slam::RangeSet>
    StaticConstantRelation_CT_Array;

  slam::RangeSet fromSet(FROMSET_SIZE);
  slam::RangeSet toSet(TOSET_SIZE);

  IndexVec relIndices;
  generateIncrementingRelations(ELEM_STRIDE, &relIndices);
  PositionType* data = &relIndices[0];    // Get a pointer to the data


  // --  Construct empty and uninitialized relation
  StaticConstantRelation_CT_Array emptyRel;
  EXPECT_FALSE( emptyRel.isValid());

  StaticConstantRelation_CT_Array uninitRel(&fromSet, &toSet);
  EXPECT_FALSE( uninitRel.isValid());


  // -- Simple relation construction
  StaticConstantRelation_CT_Array relation(&fromSet, &toSet);
  relation.bindIndices(relIndices.size(), data);

  // Note: Since the striding is a compile time constant,
  //       we don't need to set the striding,
  EXPECT_TRUE(relation.isValid(true));

  //        .. but we can, if we'd like to
  relation.bindBeginOffsets(fromSet.size(), ELEM_STRIDE);
  EXPECT_TRUE(relation.isValid(true));



  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_simple_construct");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);

  // -- Construction using set and relation builder objects
  typedef StaticConstantRelation_CT_Array::RelationBuilder RelationBuilder;
  StaticConstantRelation_CT_Array builderRel =
    RelationBuilder()
    .fromSet( &fromSet)
    .toSet( &toSet)
    .begins(
      RelationBuilder::BeginsSetBuilder()
      .stride(ELEM_STRIDE))
    .indices(
      RelationBuilder::IndicesSetBuilder()
      .size(relIndices.size())
      .data(data) )
  ;
  EXPECT_TRUE(builderRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel);
  traverseRelation_delayedSubscript(builderRel);
  iterateRelation_begin_end(builderRel);
  iterateRelation_range(builderRel);


  // Similarly here, it is ok to omit the striding in the setup
  // due to the compile time striding policy
  StaticConstantRelation_CT_Array builderRel_implicitStride =
    RelationBuilder()
    .fromSet( &fromSet)
    .toSet( &toSet)
//  .begins( RelationBuilder::BeginsSetBuilder()
//          .stride(ELEM_STRIDE))
    .indices(RelationBuilder::IndicesSetBuilder()
             .size(relIndices.size())
             .data(data) );
  EXPECT_TRUE(builderRel_implicitStride.isValid(true));


  // Test traversal of the relation data
  SCOPED_TRACE("RTStride_array_builder construct");
  traverseRelation_doubleSubscript(builderRel_implicitStride);
  traverseRelation_delayedSubscript(builderRel_implicitStride);
  iterateRelation_begin_end(builderRel_implicitStride);
  iterateRelation_range(builderRel_implicitStride);
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
