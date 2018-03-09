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
 * \file CardinalityPolicies.hpp
 *
 * \brief Cardinality policies for Slam
 *
 * Cardinality policies are meant to represent the cardinality of a relation
 * with respect to an element of a an OrderedSet, i.e., the number of elements
 * of a FromSet to which each element of a ToSet maps.
 *
 * This file implements two concrete cardinality policies:
 * * ConstantCardinality, in which every member of the FromSet maps to a fixed
 *   number of entries in the ToSet
 * * VariableCardinality, in which members of the FromSet map to an arbitrary
 *   number of entries in the ToSet
 *
 * A valid cardinality policy must support the following interface:
 *  * RelationalOperatorSizeType
 *    -- A public type that indicates the SizePolicy for the each entry in
 *        the Cardinality relation \see SizePolicies.hpp
 *  * size(ElementType idx) : const ElementType
 *    -- returns the cardinality of the relation for element with index
 *        idx of the FromSet
 *  * offset(ElementType idx) : const ElementType
 *     -- returns the offset to the first element of the ToSet for element
 *        with index idx of the FromSet
 *  * totalSize(): int
 *     -- returns the total number of elements in this relation.
 *        That is, the sum of size(idx) for each element (with index idx)
 *        of the from set.
 *  * isValid(): bool
 *     -- indicates whether the CardinalityPolicy instance is valid
 *
 */

#ifndef SLAM_POLICIES_CARDINALITY_H_
#define SLAM_POLICIES_CARDINALITY_H_

#include "axom/config.hpp"
#include "axom/Macros.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/StridePolicies.hpp"
#include "slam/OffsetPolicies.hpp"
#include "slam/IndirectionPolicies.hpp"
#include "slam/PolicyTraits.hpp"

#include "slam/OrderedSet.hpp"  // Note: Not a circular dependency since
                                // CardinalityPolicies are for relations

namespace axom
{
namespace slam
{

namespace policies
{

template<
  typename ElementType = int,
  typename StridePolicy = RuntimeStride<ElementType> >
struct ConstantCardinality
{
  typedef RuntimeSize<ElementType>                BeginsSizePolicy;
  typedef ZeroOffset<ElementType>                 BeginsOffsetPolicy;
  typedef StridePolicy BeginsStridePolicy;
  typedef NoIndirection<ElementType,ElementType>  BeginsIndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  typedef OrderedSet<
      BeginsSizePolicy,
      BeginsOffsetPolicy,
      BeginsStridePolicy >                      BeginsSet;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  typedef typename StrideToSize<
      BeginsStridePolicy,
      ElementType,
      BeginsStridePolicy::DEFAULT_VALUE>::SizeType RelationalOperatorSizeType;


  ConstantCardinality() : m_begins() {}
  ConstantCardinality(BeginsSet begins) : m_begins(begins) {}
  ConstantCardinality(ElementType fromSetSize)
  {
    m_begins = BeginsSet(fromSetSize);
  }

  ConstantCardinality(ElementType fromSetSize,
                      typename BeginsSet::SetBuilder& builder )
  {
    // needs a size and a stride (when runtime)
    builder.size(fromSetSize);
    m_begins = builder;
  }

  const ElementType size(ElementType AXOM_NOT_USED(fromPos) ) const
  {
    return m_begins.stride();
  }

  const ElementType offset(ElementType fromPos) const
  {
    return m_begins[fromPos];
  }

  void bindBeginOffsets(ElementType fromSetSize, ElementType stride)
  {
    m_begins = typename BeginsSet::SetBuilder()
               .size(fromSetSize)
               .stride(stride);
  }

  ElementType totalSize() const
  {
    return m_begins.stride() * m_begins.size();
  }

  template<typename FromSetType>
  bool isValid(const FromSetType* fromSet,
               bool AXOM_NOT_USED(vertboseOutput) = false) const
  {
    return m_begins.size() == fromSet->size();
  }


  BeginsSet m_begins;
};

template<
  typename ElementType = int,
  typename IndirectionPolicy = STLVectorIndirection<ElementType, ElementType>
  >
struct VariableCardinality
{
  typedef RuntimeSize<ElementType>  BeginsSizePolicy;
  typedef ZeroOffset<ElementType>   BeginsOffsetPolicy;
  typedef StrideOne<ElementType>    BeginsStridePolicy;
  typedef IndirectionPolicy BeginsIndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  typedef OrderedSet<
      BeginsSizePolicy,
      BeginsOffsetPolicy,
      BeginsStridePolicy,
      IndirectionPolicy>                BeginsSet;


  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  typedef BeginsSizePolicy RelationalOperatorSizeType;

  typedef typename IndirectionPolicy::IndirectionBufferType
    IndirectionBufferType;

  VariableCardinality() : m_begins() {}
  VariableCardinality(BeginsSet begins) : m_begins(begins) {}
  VariableCardinality(ElementType fromSetSize,
                      typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionBufferType* data)
  {
    m_begins = typename BeginsSet::SetBuilder()
               .size(fromSetSize + 1)
               .data(data);
  }

  const ElementType size(ElementType fromPos) const
  {
    return offset(fromPos + 1) - offset(fromPos);
  }

  const ElementType offset(ElementType fromPos) const
  {
    return m_begins[fromPos];
  }

  ElementType totalSize() const
  {
    return m_begins.empty()
           ? ElementType()
           : offset(m_begins.size() - 1);
  }

  template<typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    return m_begins.size() == (fromSet->size() + 1)
           && static_cast<IndirectionPolicy>(m_begins).isValid(
      m_begins.size(), m_begins.offset(), m_begins.stride(), verboseOutput);
  }


  BeginsSet m_begins;
};

} // end namespace policies

} // end namespace slam
} // end namespace axom


#endif // SLAM_POLICIES_CARDINALITY_H_
