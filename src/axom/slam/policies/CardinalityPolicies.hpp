// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
#include "axom/core/Macros.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/slam/OrderedSet.hpp"  // Note: Not a circular dependency since
// CardinalityPolicies are for relations

namespace axom
{
namespace slam
{
namespace policies
{
template <typename ElementType = int, typename StridePolicy = RuntimeStride<ElementType>>
struct ConstantCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StridePolicy;
  using BeginsIndirectionPolicy = NoIndirection<ElementType, ElementType>;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  using BeginsSet =
    OrderedSet<ElementType, ElementType, BeginsSizePolicy, BeginsOffsetPolicy, BeginsStridePolicy>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType =
    typename StrideToSize<BeginsStridePolicy,
                          ElementType,
                          BeginsStridePolicy::DEFAULT_VALUE>::SizeType;

  ConstantCardinality() : m_begins() { }
  ConstantCardinality(BeginsSet begins) : m_begins(begins) { }
  ConstantCardinality(ElementType fromSetSize)
  {
    m_begins = BeginsSet(fromSetSize);
  }

  ConstantCardinality(ElementType fromSetSize,
                      typename BeginsSet::SetBuilder& builder)
  {
    // needs a size and a stride (when runtime)
    builder.size(fromSetSize);
    m_begins = builder;
  }

  const ElementType size(ElementType AXOM_NOT_USED(fromPos)) const
  {
    return m_begins.stride();
  }

  const ElementType offset(ElementType fromPos) const
  {
    return m_begins[fromPos];
  }

  void bindBeginOffsets(ElementType fromSetSize, ElementType stride)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize).stride(stride);
  }

  ElementType totalSize() const { return m_begins.stride() * m_begins.size(); }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet,
               bool AXOM_NOT_USED(vertboseOutput) = false) const
  {
    return m_begins.size() == fromSet->size();
  }

  BeginsSet m_begins;
};

template <typename ElementType = int,
          typename IndirectionPolicy = STLVectorIndirection<ElementType, ElementType>>
struct VariableCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StrideOne<ElementType>;
  using BeginsIndirectionPolicy = IndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset
  using BeginsSet = OrderedSet<ElementType,
                               ElementType,
                               BeginsSizePolicy,
                               BeginsOffsetPolicy,
                               BeginsStridePolicy,
                               IndirectionPolicy>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType = BeginsSizePolicy;

  using IndirectionBufferType = typename IndirectionPolicy::IndirectionBufferType;

  VariableCardinality() : m_begins() { }
  VariableCardinality(BeginsSet begins) : m_begins(begins) { }
  VariableCardinality(ElementType fromSetSize,
                      typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionBufferType* data)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize + 1).data(data);
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
    return m_begins.empty() ? ElementType() : offset(m_begins.size() - 1);
  }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    return m_begins.size() == (fromSet->size() + 1) &&
      static_cast<IndirectionPolicy>(m_begins).isValid(m_begins.size(),
                                                       m_begins.offset(),
                                                       m_begins.stride(),
                                                       verboseOutput);
  }

  BeginsSet m_begins;
};

}  // end namespace policies

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_CARDINALITY_H_
