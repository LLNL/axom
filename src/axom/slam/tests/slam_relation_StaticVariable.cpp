// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file slam_relation_StaticVariable.cpp
 *
 * \brief Unit tests for Slam's StaticRelation class
 *  configured with variable per-element cardinality
 */

#include <iostream>
#include <iterator>
#include <sstream>

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/ModularInt.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/Relation.hpp"

#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/RelationSet.hpp"
#include "axom/slam/StaticRelation.hpp"

namespace
{
namespace slam = axom::slam;
namespace policies = axom::slam::policies;

using SetPosition = slam::DefaultPositionType;
using SetElement = slam::DefaultElementType;

using RangeSetType = slam::RangeSet<SetPosition, SetElement>;
using RelationType = slam::Relation<SetPosition, SetElement>;

using IndexVec = std::vector<SetPosition>;

constexpr SetPosition FROMSET_SIZE = 7;
constexpr SetPosition TOSET_SIZE = 8;

using STLIndirection = policies::STLVectorIndirection<SetPosition, SetElement>;

using VariableCardinality = policies::VariableCardinality<SetPosition, STLIndirection>;

using StaticVariableRelationType =
  slam::StaticRelation<SetPosition, SetElement, VariableCardinality, STLIndirection, RangeSetType, RangeSetType>;

using MappedVariableCardinality = policies::MappedVariableCardinality<SetPosition, STLIndirection>;

using StaticMappedVariableRelationType = slam::StaticRelation<SetPosition,
                                                               SetElement,
                                                               MappedVariableCardinality,
                                                               STLIndirection,
                                                               RangeSetType,
                                                               RangeSetType>;

// Use a slam::ModularInt type for more interesting test data
using CTSize = policies::CompileTimeSize<int, TOSET_SIZE>;
using FixedModularInt = slam::ModularInt<CTSize>;

SetPosition elementCardinality(SetPosition fromPos) { return fromPos; }

SetPosition relationData(SetPosition fromPos, SetPosition toPos)
{
  return FixedModularInt(fromPos + toPos);
}

template <typename StrType, typename VecType>
void printVector(StrType const& msg, VecType const& vec)
{
  std::stringstream sstr;

  sstr << "\n** " << msg << "\n\t";
  sstr << "Array of size " << vec.size() << ": ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<SetPosition>(sstr, " "));

  SLIC_INFO(sstr.str());
}

template <typename VecType>
void generateIncrementingRelations(VecType* begins, VecType* indices)
{
  VecType& beginsVec = *begins;
  VecType& indicesVec = *indices;

  SetPosition curIdx = SetPosition();

  for(SetPosition i = 0; i < FROMSET_SIZE; ++i)
  {
    beginsVec.push_back(curIdx);
    for(SetPosition j = 0; j < elementCardinality(i); ++j)
    {
      indicesVec.push_back(relationData(i, j));
      ++curIdx;
    }
  }
  beginsVec.push_back(curIdx);
}

/**
 * \brief Traverses the relation's entities using the double subscript access API
 *
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void traverseRelation_doubleSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using double subscript: ");
  for(SetPosition fromPos = 0; fromPos < rel.fromSet()->size(); ++fromPos)
  {
    const SetPosition fromSize = rel.size(fromPos);
    EXPECT_EQ(elementCardinality(fromPos), fromSize);

    for(int toPos = 0; toPos < fromSize; ++toPos)
    {
      SetPosition actualVal = rel[fromPos][toPos];

      EXPECT_EQ(relationData(fromPos, toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using separated subscript operators
 *
 * The first subscript operator gets the set of entities in the ToSet
 * that are mapped to the given element of the relation's FromSet
 *
 * \note Assumes that the relation data has been set using
 * the generateIncrementingRelations() function
 */
template <typename RelationType>
void traverseRelation_delayedSubscript(RelationType& rel)
{
  SLIC_INFO("Traversing relation data using delayed second subscript: ");
  for(SetPosition fromPos = 0; fromPos < rel.fromSet()->size(); ++fromPos)
  {
    const SetPosition fromSize = rel.size(fromPos);
    EXPECT_EQ(elementCardinality(fromPos), fromSize);

    typename RelationType::RelationSubset set = rel[fromPos];
    for(int toPos = 0; toPos < set.size(); ++toPos)
    {
      SetPosition actualVal = rel[fromPos][toPos];
      EXPECT_EQ(relationData(fromPos, toPos), actualVal);
    }
  }
}

/**
 * \brief Traverses relation using the iterator API (begin()/end() )
 *
 * \note The iterator API depends on C++11
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void iterateRelation_begin_end(RelationType& rel)
{
  using FromSet = typename RelationType::FromSetType;
  using FromSetIter = typename FromSet::iterator;
  using RelIter = typename RelationType::RelationIterator;

  SLIC_INFO("Traversing relation data using iterator begin()/end() functions");
  for(FromSetIter sIt = rel.fromSet()->begin(), sItEnd = rel.fromSet()->end(); sIt != sItEnd; ++sIt)
  {
    SetPosition actualSize = rel.size(*sIt);

    SetPosition fromSetEltNum = std::distance(rel.fromSet()->begin(), sIt);
    EXPECT_EQ(elementCardinality(fromSetEltNum), actualSize);

    RelIter toSetBegin = rel.begin(*sIt);
    RelIter toSetEnd = rel.end(*sIt);
    for(RelIter relIt = toSetBegin; relIt != toSetEnd; ++relIt)
    {
      SetPosition toSetEltNum = std::distance(toSetBegin, relIt);
      ASSERT_EQ(relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
}

/**
 * \brief Traverses relation using the iterator range API
 *
 * \note The iterator API depends on C++11
 * \note Assumes that the relation data has been set
 * using the generateIncrementingRelations() function
 */
template <typename RelationType>
void iterateRelation_range(RelationType& rel)
{
  using FromSet = typename RelationType::FromSetType;
  using FromSetIter = typename FromSet::iterator;
  using FromSetIterPair = typename FromSet::iterator_pair;

  using RelIter = typename RelationType::RelationIterator;
  using RelIterPair = typename RelationType::RelationIteratorPair;

  SLIC_INFO("Traversing relation data using iterator range() functions");
  FromSetIterPair itPair = rel.fromSet()->range();
  for(FromSetIter sIt = itPair.first; sIt != itPair.second; ++sIt)
  {
    SetPosition fromSetEltNum = std::distance(itPair.first, sIt);

    RelIterPair toSetItPair = rel.range(*sIt);
    for(RelIter relIt = toSetItPair.first; relIt != toSetItPair.second; ++relIt)
    {
      SetPosition toSetEltNum = std::distance(toSetItPair.first, relIt);
      ASSERT_EQ(relationData(fromSetEltNum, toSetEltNum), *relIt);
    }
  }
}

}  // end anonymous namespace

TEST(slam_static_variable_relation, construct_empty)
{
  SLIC_INFO("Testing empty relation.  isValid() should be false.");

  StaticVariableRelationType emptyRel;

  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_static_variable_relation, construct_uninitialized)
{
  SLIC_INFO("Testing uninitialized relation.  isValid() should be false.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  StaticVariableRelationType emptyRel(&fromSet, &toSet);

  EXPECT_FALSE(emptyRel.isValid(true));
}

TEST(slam_static_variable_relation, construct_relation)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec relOffsets;
  IndexVec relIndices;

  SLIC_INFO("Uninitialized relation data");
  printVector("begins vector", relOffsets);
  printVector("indices vector", relIndices);
  generateIncrementingRelations(&relOffsets, &relIndices);

  SLIC_INFO("Initialized relation data");
  printVector("begins vector", relOffsets);
  printVector("indices vector", relIndices);

  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

  EXPECT_TRUE(incrementingRel.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("Vector_simple_construct");
  traverseRelation_doubleSubscript(incrementingRel);
  traverseRelation_delayedSubscript(incrementingRel);
  iterateRelation_begin_end(incrementingRel);
  iterateRelation_range(incrementingRel);
}

TEST(slam_static_variable_relation, construct_builder)
{
  SLIC_INFO("Testing construction using builder interface.");

  RangeSetType fromSet(FROMSET_SIZE);
  RangeSetType toSet(TOSET_SIZE);

  IndexVec offsets;
  IndexVec relIndices;
  generateIncrementingRelations(&offsets, &relIndices);

  using RelationBuilder = StaticVariableRelationType::RelationBuilder;
  StaticVariableRelationType relation = RelationBuilder()
                                          .fromSet(&fromSet)
                                          .toSet(&toSet)
                                          .begins(RelationBuilder::BeginsSetBuilder()  //
                                                    .size(offsets.size())              //
                                                    .data(&offsets))
                                          .indices(RelationBuilder::IndicesSetBuilder()  //
                                                     .size(relIndices.size())            //
                                                     .data(&relIndices));
  EXPECT_TRUE(relation.isValid(true));

  // Test traversal of the relation data
  SCOPED_TRACE("Vector_builder");
  traverseRelation_doubleSubscript(relation);
  traverseRelation_delayedSubscript(relation);
  iterateRelation_begin_end(relation);
  iterateRelation_range(relation);
}

TEST(slam_static_variable_relation, empty_relation_out_of_bounds)
{
  StaticVariableRelationType emptyRel;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(emptyRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

/// Tests for out-of-bounds access

TEST(slam_static_variable_relation, initialized_rel_out_of_bounds)
{
  SLIC_INFO("Testing simple incrementing relation.  isValid() should be true.");

  IndexVec relOffsets, relIndices;
  generateIncrementingRelations(&relOffsets, &relIndices);

  RangeSetType fromSet(FROMSET_SIZE), toSet(TOSET_SIZE);
  StaticVariableRelationType incrementingRel(&fromSet, &toSet);
  incrementingRel.bindBeginOffsets(fromSet.size(), &relOffsets);
  incrementingRel.bindIndices(relIndices.size(), &relIndices);

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode,
  // so this test will only fail in debug mode
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(incrementingRel[FROMSET_SIZE], "");
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif
}

//----------------------------------------------------------------------

template <typename Relation>
void expectConsistentRelationIndexing(Relation& relation,
                                     const IndexVec& begins,
                                     const IndexVec& indices)
{
  ASSERT_TRUE(relation.isValid(true));
  slam::RelationSet<Relation> relationSet(&relation);
  ASSERT_TRUE(relationSet.isValid(true));

  for(SetPosition row = 0; row < relation.fromSetSize(); ++row)
  {
    ASSERT_EQ(relation.size(row), begins[row + 1] - begins[row]);
    for(SetPosition flat = begins[row]; flat < begins[row + 1]; ++flat)
    {
      const auto coordinate = relationSet.at(flat);
      EXPECT_EQ(relation.firstIndex(flat), row);
      EXPECT_EQ(coordinate.first, row);
      EXPECT_EQ(coordinate.second, indices[flat]);
      EXPECT_EQ(relation[row][flat - begins[row]], indices[flat]);
    }
  }

  EXPECT_EQ(relation.firstIndex(SetPosition {-1}), SetPosition {-1});
  EXPECT_EQ(relation.firstIndex(static_cast<SetPosition>(indices.size())), SetPosition {-1});
}

TEST(slam_static_variable_relation, consistent_indexing_plain_and_mapped)
{
  RangeSetType fromSet(6), toSet(10);
  IndexVec begins {0, 0, 3, 3, 4, 7, 7};  // rows 0, 2 and 5 are empty
  IndexVec indices {2, 4, 6, 1, 3, 5, 7};

  StaticVariableRelationType plain(&fromSet, &toSet);
  plain.bindBeginOffsets(fromSet.size(), &begins);
  plain.bindIndices(static_cast<SetPosition>(indices.size()), &indices);
  expectConsistentRelationIndexing(plain, begins, indices);

  IndexVec firstIndices(indices.size(), SetPosition {-1});
  StaticMappedVariableRelationType mapped(&fromSet, &toSet);
  mapped.bindBeginOffsets(fromSet.size(), &begins);
  mapped.bindIndices(static_cast<SetPosition>(indices.size()), &indices);
  mapped.bindFirstIndices(static_cast<SetPosition>(firstIndices.size()), &firstIndices);
  expectConsistentRelationIndexing(mapped, begins, indices);
}

TEST(slam_static_variable_relation, empty_and_single_row_relations)
{
  RangeSetType emptyFromSet(0), toSet(4);
  IndexVec zeroRowsBegins {0};
  IndexVec emptyIndices;

  StaticVariableRelationType zeroRows(&emptyFromSet, &toSet);
  zeroRows.bindBeginOffsets(emptyFromSet.size(), &zeroRowsBegins);
  zeroRows.bindIndices(0, &emptyIndices);
  expectConsistentRelationIndexing(zeroRows, zeroRowsBegins, emptyIndices);

  IndexVec emptyFirstIndices;
  StaticMappedVariableRelationType mappedZeroRows(&emptyFromSet, &toSet);
  mappedZeroRows.bindBeginOffsets(emptyFromSet.size(), &zeroRowsBegins);
  mappedZeroRows.bindIndices(0, &emptyIndices);
  mappedZeroRows.bindFirstIndices(0, &emptyFirstIndices);
  expectConsistentRelationIndexing(mappedZeroRows, zeroRowsBegins, emptyIndices);

  RangeSetType threeRowsSet(3);
  IndexVec allEmptyBegins {0, 0, 0, 0};
  StaticVariableRelationType allEmpty(&threeRowsSet, &toSet);
  allEmpty.bindBeginOffsets(threeRowsSet.size(), &allEmptyBegins);
  allEmpty.bindIndices(0, &emptyIndices);
  expectConsistentRelationIndexing(allEmpty, allEmptyBegins, emptyIndices);

  StaticMappedVariableRelationType mappedAllEmpty(&threeRowsSet, &toSet);
  mappedAllEmpty.bindBeginOffsets(threeRowsSet.size(), &allEmptyBegins);
  mappedAllEmpty.bindIndices(0, &emptyIndices);
  mappedAllEmpty.bindFirstIndices(0, &emptyFirstIndices);
  expectConsistentRelationIndexing(mappedAllEmpty, allEmptyBegins, emptyIndices);

  RangeSetType oneRowSet(1);
  IndexVec oneRowBegins {0, 2};
  IndexVec oneRowIndices {1, 3};
  StaticVariableRelationType oneRow(&oneRowSet, &toSet);
  oneRow.bindBeginOffsets(oneRowSet.size(), &oneRowBegins);
  oneRow.bindIndices(static_cast<SetPosition>(oneRowIndices.size()), &oneRowIndices);
  expectConsistentRelationIndexing(oneRow, oneRowBegins, oneRowIndices);
}

TEST(slam_static_variable_relation, rejects_malformed_begin_offsets)
{
  RangeSetType fromSet(3), toSet(8);
  IndexVec indices {0, 1, 2};

  const std::vector<IndexVec> malformedBegins {
    {1, 1, 2, 3},   // leading gap
    {0, -1, 1, 3},  // negative offset
    {0, 2, 1, 3},   // decreasing offsets
    {0, 1, 2, 2},   // terminal does not equal the index count
  };

  for(auto begins : malformedBegins)
  {
    StaticVariableRelationType relation(&fromSet, &toSet);
    relation.bindBeginOffsets(fromSet.size(), &begins);
    relation.bindIndices(static_cast<SetPosition>(indices.size()), &indices);
    EXPECT_FALSE(relation.isValid()) << "accepted malformed begins array";
  }

  IndexVec nonzeroTerminal {0, 0, 0, 1};
  IndexVec noIndices;
  StaticVariableRelationType emptyData(&fromSet, &toSet);
  emptyData.bindBeginOffsets(fromSet.size(), &nonzeroTerminal);
  emptyData.bindIndices(0, &noIndices);
  EXPECT_FALSE(emptyData.isValid());

  IndexVec leadingGap {2, 2, 2, 3};
  StaticVariableRelationType gapped(&fromSet, &toSet);
  gapped.bindBeginOffsets(fromSet.size(), &leadingGap);
  gapped.bindIndices(static_cast<SetPosition>(indices.size()), &indices);
  EXPECT_EQ(gapped.firstIndex(0), SetPosition {-1});
  EXPECT_EQ(gapped.firstIndex(1), SetPosition {-1});
}

TEST(slam_static_variable_relation, mapped_inverse_must_agree_with_rows)
{
  RangeSetType fromSet(3), toSet(8);
  IndexVec begins {0, 1, 1, 3};
  IndexVec indices {2, 4, 6};
  IndexVec wrongFirstIndices {0, 0, 0};

  StaticMappedVariableRelationType relation(&fromSet, &toSet);
  relation.bindBeginOffsets(fromSet.size(), &begins);
  relation.bindIndices(static_cast<SetPosition>(indices.size()), &indices);
  relation.bindFirstIndices(static_cast<SetPosition>(wrongFirstIndices.size()),
                            &wrongFirstIndices,
                            false);
  EXPECT_FALSE(relation.isValid());

  IndexVec shortFirstIndices(2, 0);
  relation.bindFirstIndices(static_cast<SetPosition>(shortFirstIndices.size()),
                            &shortFirstIndices,
                            false);
  EXPECT_FALSE(relation.isValid());
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
