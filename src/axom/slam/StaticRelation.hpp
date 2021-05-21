// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file StaticRelation.hpp
 *
 * \brief API for a topological relation between two sets where the
 *        relation does not change after it is initialized
 *
 */

#ifndef SLAM_STATIC_RELATION_HPP_
#define SLAM_STATIC_RELATION_HPP_

#include "axom/config.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/slam/OrderedSet.hpp"
#include "axom/slam/Relation.hpp"

namespace axom
{
namespace slam
{
template <typename PosType,   // = slam::DefaultPositionType,
          typename ElemType,  // = slam::DefaultElementType,
          typename RelationCardinalityPolicy,
          typename RelationIndicesIndirectionPolicy,
          typename TheFromSet = Set<PosType, ElemType>,
          typename TheToSet = Set<PosType, ElemType>>
class StaticRelation : public /*Relation,*/ RelationCardinalityPolicy
{
public:
  using SetPosition = PosType;
  using SetElement = ElemType;

  using FromSetType = TheFromSet;
  using ToSetType = TheToSet;

  using CardinalityPolicy = RelationCardinalityPolicy;
  using BeginsSizePolicy = typename CardinalityPolicy::RelationalOperatorSizeType;

  using IndicesIndirectionPolicy = RelationIndicesIndirectionPolicy;

  using RelationSubset = OrderedSet<SetPosition,
                                    SetElement,
                                    BeginsSizePolicy,
                                    policies::RuntimeOffset<SetPosition>,
                                    policies::StrideOne<SetPosition>,
                                    IndicesIndirectionPolicy>;

  using IndicesSet = OrderedSet<SetPosition,
                                SetElement,
                                policies::RuntimeSize<SetPosition>,
                                policies::ZeroOffset<SetPosition>,
                                policies::StrideOne<SetPosition>,
                                IndicesIndirectionPolicy>;

  using IndirectionBufferType =
    typename IndicesIndirectionPolicy::IndirectionBufferType;

  // types for iterator
  using RelationIterator = typename RelationSubset::iterator;
  using RelationIteratorPair = typename RelationSubset::iterator_pair;

  using RelationConstIterator = typename RelationSubset::const_iterator;
  using RelationConstIteratorPair = typename RelationSubset::const_iterator_pair;

public:
  struct RelationBuilder;

  StaticRelation()
    : m_fromSet(EmptySetTraits<FromSetType>::emptySet())
    , m_toSet(EmptySetTraits<ToSetType>::emptySet())
  { }

  StaticRelation(FromSetType* fromSet, ToSetType* toSet)
    : CardinalityPolicy(
        EmptySetTraits<FromSetType>::isEmpty(fromSet) ? 0 : fromSet->size())
    , m_fromSet(fromSet)
    , m_toSet(toSet)
  { }

  StaticRelation(const RelationBuilder& builder)
    : CardinalityPolicy(builder.m_cardPolicy)
    , m_fromSet(builder.m_fromSet)
    , m_toSet(builder.m_toSet)
    , m_relationIndices(builder.m_indBuilder)
  { }

  struct RelationBuilder
  {
    friend class StaticRelation;

    using BeginsSetBuilder =
      typename StaticRelation::CardinalityPolicy::BeginsSet::SetBuilder;
    using IndicesSetBuilder = typename StaticRelation::IndicesSet::SetBuilder;

    RelationBuilder()
      : m_fromSet(EmptySetTraits<FromSetType>::emptySet())
      , m_toSet(EmptySetTraits<ToSetType>::emptySet())
    { }

    RelationBuilder& fromSet(FromSetType* pFromSet)
    {
      m_fromSet = pFromSet;
      if(m_cardPolicy.totalSize() == 0 &&
         !EmptySetTraits<FromSetType>::isEmpty(m_fromSet))
      {
        m_cardPolicy = CardinalityPolicy(m_fromSet->size());
      }
      return *this;
    }

    RelationBuilder& toSet(ToSetType* pToSet)
    {
      m_toSet = pToSet;
      return *this;
    }

    RelationBuilder& begins(BeginsSetBuilder& beginsBuilder)
    {
      SLIC_ASSERT_MSG(
        !EmptySetTraits<FromSetType>::isEmpty(m_fromSet),
        "Must set the 'fromSet' pointer before setting the begins set");

      m_cardPolicy = CardinalityPolicy(m_fromSet->size(), beginsBuilder);
      return *this;
    }
    RelationBuilder& indices(const IndicesSetBuilder& indicesBuilder)
    {
      m_indBuilder = indicesBuilder;
      return *this;
    }

  private:
    FromSetType* m_fromSet;
    ToSetType* m_toSet;
    CardinalityPolicy m_cardPolicy;
    IndicesSetBuilder m_indBuilder;
  };

public:
  const RelationSubset operator[](SetPosition fromSetInd) const
  {
    SLIC_ASSERT(m_relationIndices.isValid(true));

    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(CardinalityPolicy::size(fromSetInd))
      .offset(CardinalityPolicy::offset(fromSetInd))
      .data(m_relationIndices.data());
  }

  RelationSubset operator[](SetPosition fromSetInd)
  {
    SLIC_ASSERT(m_relationIndices.isValid(true));

    using SetBuilder = typename RelationSubset::SetBuilder;
    return SetBuilder()
      .size(CardinalityPolicy::size(fromSetInd))
      .offset(CardinalityPolicy::offset(fromSetInd))
      .data(m_relationIndices.data());
  }

  bool isValid(bool verboseOutput = false) const;

  RelationIterator begin(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].begin();
  }

  RelationConstIterator begin(SetPosition fromSetInd) const
  {
    return (*this)[fromSetInd].begin();
  }

  RelationIterator end(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].end();
  }

  RelationConstIterator end(SetPosition fromSetInd) const
  {
    return (*this)[fromSetInd].end();
  }

  RelationIteratorPair range(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].range();
  }

  RelationConstIteratorPair range(SetPosition fromSetInd) const
  {
    return (*this)[fromSetInd].range();
  }

  bool hasFromSet() const
  {
    return !EmptySetTraits<FromSetType>::isEmpty(m_fromSet);
  }
  FromSetType* fromSet() { return m_fromSet; }
  const FromSetType* fromSet() const { return m_fromSet; }

  bool hasToSet() const { return !EmptySetTraits<ToSetType>::isEmpty(m_toSet); }
  ToSetType* toSet() { return m_toSet; }
  const ToSetType* toSet() const { return m_toSet; }

  SetPosition fromSetSize() { return m_fromSet->size(); }

  SetPosition toSetSize() { return m_toSet->size(); }

  void bindIndices(SetPosition size, IndirectionBufferType* data)
  {
    m_relationIndices = typename IndicesSet::SetBuilder().size(size).data(data);
  }

  const IndirectionBufferType* relationData() const
  {
    return m_relationIndices.data();
  }

  IndirectionBufferType* relationData() { return m_relationIndices.data(); }

private:
  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  IndicesSet m_relationIndices;
};

/**
 * \brief Checks whether the relation is valid
 *
 * A relation is valid when:
 * * Its fromSet and toSet are not null
 * * The CardinalityPolicy is valid.
 *   This implies that for each element, pos, of the fromSet,
 *   it is valid to call rel.size(pos), rel.offset()
 *   It is also valid to call rel.totalSize()
 *
 *
 * @return True if the relation is valid, false otherwise
 */
template <typename PosType,
          typename ElemType,
          typename RelationCardinalityPolicy,
          typename RelationIndicesIndirectionPolicy,
          typename FromSetType,
          typename ToSetType>
bool StaticRelation<PosType,
                    ElemType,
                    RelationCardinalityPolicy,
                    RelationIndicesIndirectionPolicy,
                    FromSetType,
                    ToSetType>::isValid(bool verboseOutput) const
{
  std::stringstream errSstr;

  bool setsAreValid = true;
  bool cardinalityIsValid = true;
  bool relationdataIsValid = true;

  // Step 1: Check if the sets are valid
  bool isFromSetNull = (m_fromSet == nullptr);
  bool isToSetNull = (m_toSet == nullptr);

  if(isFromSetNull || isToSetNull)
  {
    if(verboseOutput)
    {
      errSstr << "\n\t Static relations require both the fromSet"
              << " and toSet to be non-null"
              << "\n\t -- fromSet was " << (isFromSetNull ? "" : " not ")
              << "null"
              << "\n\t -- toSet was " << (isToSetNull ? "" : " not ") << "null";
    }

    setsAreValid = false;
  }

  // Step 2: Check if the cardinality is valid
  if(setsAreValid)
  {
    if(!CardinalityPolicy::isValid(m_fromSet, verboseOutput))
    {
      if(verboseOutput)
      {
        errSstr << "\n\t Invalid cardinality state.";
        //  TODO -- improve this message ;
      }
      cardinalityIsValid = false;
    }
  }

  // Step 3: Check if the relation data is valid
  if(cardinalityIsValid)
  {
    if(m_relationIndices.size() != this->totalSize())
    {
      if(verboseOutput)
      {
        errSstr << "\n\t* relation indices has the wrong size."
                << "\n\t-- from set size is: " << m_fromSet->size()
                << "\n\t-- expected relation size: " << this->totalSize()
                << "\n\t-- actual size: " << m_relationIndices.size();
      }
      relationdataIsValid = false;
    }

    if(!m_relationIndices.empty())
    {
      // Check that all begins offsets are in the right range
      // Specifically, they must be in the index space of m_relationIndices
      for(SetPosition pos = 0; pos < m_fromSet->size(); ++pos)
      {
        SetPosition off = this->offset(pos);
        if(!m_relationIndices.isValidIndex(off) && off != m_relationIndices.size())
        {
          if(verboseOutput)
          {
            errSstr << "\n\t* Begin offset for index " << pos
                    << " was out of range."
                    << "\n\t-- value: " << this->offset(pos)
                    << " needs to be within range [0,"
                    << m_relationIndices.size() << "]";
          }
          relationdataIsValid = false;
        }
      }
    }

    // Check that all relation indices are in range for m_toSet
    for(SetPosition pos = 0; pos < m_relationIndices.size(); ++pos)
    {
      if(!m_toSet->isValidIndex(m_relationIndices[pos]))
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* Relation index was out of range."
                  << "\n\t-- value: " << m_relationIndices[pos]
                  << " needs to be in range [0," << m_toSet->size() << ")";
        }
        relationdataIsValid = false;
      }
    }
  }

  // We are done.  Output the messages if applicable and return
  bool bValid = setsAreValid && cardinalityIsValid && relationdataIsValid;

  if(verboseOutput && !bValid)
  {
    SLIC_DEBUG(errSstr.str());
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_STATIC_RELATION_HPP_
