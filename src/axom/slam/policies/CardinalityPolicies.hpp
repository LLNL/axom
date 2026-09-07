// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file CardinalityPolicies.hpp
 *
 * \brief Cardinality and begin-offset policies for Slam relations.
 *
 * A cardinality policy records how many to-set positions are associated with
 * each from-set position and where those entries begin in the relation's storage.
 * ConstantCardinality computes these offsets from a fixed stride.
 * VariableCardinality reads begin offsets from a buffer. MappedVariableCardinality
 * also binds a buffer for constant-time lookup of the from-set position.
 *
 * The policies provide:
 * - RelationalOperatorSizeType, the size policy for each related subset.
 * - size(from), the number of entries associated with a from-set position.
 * - offset(from), the begin offset in the relation's flat storage.
 * - firstIndex(flat), the from-set position associated with a flat position.
 * - totalSize(), the sum of the per-element cardinalities.
 * - isValid(fromSet, verbose), a check of the policy's stored state.
 *
 * Begin offsets use the relation's flat position type. The connectivity buffer
 * separately stores to-set positions, which may have a different type.
 *
 */

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/PolicyTraits.hpp"

// Note: Not a circular dependency since CardinalityPolicies are for relations
#include "axom/slam/OrderedSet.hpp"

#include <limits>
#include <utility>

namespace axom::slam::policies
{
namespace detail
{
/// \brief Begin offsets must start at zero and be nonnegative and nondecreasing.
template <typename BeginsSet>
bool hasValidBeginOffsetValues(const BeginsSet& begins)
{
  using PositionType = typename BeginsSet::PositionType;

  if(begins.empty() || begins[PositionType {}] != PositionType {})
  {
    return false;
  }

  PositionType previous = begins[PositionType {}];
  for(PositionType pos = PositionType {1}; pos < begins.size(); ++pos)
  {
    const PositionType current = begins[pos];
    if(current < PositionType {} || current < previous)
    {
      return false;
    }
    previous = current;
  }

  return true;
}

/// \brief Check begin-offset storage and values against the from-set size.
template <typename BeginsSet, typename FromSet>
bool hasValidBeginOffsets(const BeginsSet& begins, const FromSet* fromSet, bool verboseOutput)
{
  using PositionType = typename BeginsSet::PositionType;

  if(fromSet == nullptr)
  {
    return false;
  }

  const auto fromSize = fromSet->size();
  if(fromSize < decltype(fromSize) {} || !std::in_range<PositionType>(fromSize))
  {
    return false;
  }

  const PositionType flatFromSize = static_cast<PositionType>(fromSize);
  return flatFromSize != std::numeric_limits<PositionType>::max() &&
    begins.size() == flatFromSize + PositionType {1} && begins.isValid(verboseOutput) &&
    hasValidBeginOffsetValues(begins);
}
}  // namespace detail

/*!
 * \class ConstantCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to a fixed number of elements in the second set
 *
 * \tparam ElementType The flat position type used for counts and offsets.
 * \tparam StridePolicy policy for number of elements being mapped
 */
template <typename ElementType = int, typename StridePolicy = RuntimeStride<ElementType>>
struct ConstantCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StridePolicy;
  using BeginsIndirectionPolicy = NoIndirection<ElementType, ElementType>;

  // runtime size (fromSet.size()), striding from template parameter, no offset.
  // The concrete interface avoids virtual dispatch for begin-offset access.
  using BeginsSet = OrderedSet<ElementType,
                               ElementType,
                               BeginsSizePolicy,
                               BeginsOffsetPolicy,
                               BeginsStridePolicy,
                               BeginsIndirectionPolicy,
                               NoSubset,
                               ConcreteInterface>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType =
    typename StrideToSize<BeginsStridePolicy, ElementType, BeginsStridePolicy::DEFAULT_VALUE>::SizeType;

  using IndirectionPtrType = typename BeginsIndirectionPolicy::IndirectionPtrType;

  ConstantCardinality() : m_begins() { }
  ConstantCardinality(BeginsSet begins) : m_begins(begins) { }
  ConstantCardinality(ElementType fromSetSize) { m_begins = BeginsSet(fromSetSize); }

  ConstantCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    // needs a size and a stride (when runtime)
    builder.size(fromSetSize);
    m_begins = builder;
  }

  AXOM_HOST_DEVICE ElementType size(ElementType AXOM_UNUSED_PARAM(fromPos)) const
  {
    return m_begins.stride();
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  AXOM_HOST_DEVICE ElementType firstIndex(ElementType offset) const
  {
    return offset / m_begins.stride();
  }

  IndirectionPtrType offsetData() { return m_begins.ptr(); }

  const IndirectionPtrType offsetData() const { return m_begins.ptr(); }

  void bindBeginOffsets(ElementType fromSetSize, ElementType stride)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize).stride(stride);
  }

  ElementType totalSize() const { return m_begins.stride() * m_begins.size(); }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    if(fromSet == nullptr || m_begins.size() != fromSet->size() ||
       !m_begins.isValid(verboseOutput) || m_begins.stride() <= ElementType {})
    {
      return false;
    }

    const ElementType size = m_begins.size();
    const ElementType stride = m_begins.stride();
    return size == ElementType {} || stride <= std::numeric_limits<ElementType>::max() / size;
  }

  BeginsSet m_begins;
};

/*!
 * \class VariableCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to an arbitrary number of elements in the second set.
 *
 * The begin-offset buffer has one entry per from-set element plus a final offset
 * equal to the total relation size. Offsets start at zero and are non-decreasing.
 * Equal adjacent offsets describe an element with no related to-set positions.
 *
 * \tparam ElementType The flat position type used for counts and offsets.
 * \tparam IndirectionPolicy How begin offsets are accessed. ArrayIndirection
 *  borrows an axom::Array object, ArrayViewIndirection stores a borrowed view,
 *  and STLVectorIndirection borrows a host-side std::vector.
 */
template <typename ElementType = int, typename IndirectionPolicy = ArrayIndirection<ElementType, ElementType>>
struct VariableCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StrideOne<ElementType>;
  using BeginsIndirectionPolicy = IndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset.
  // The concrete interface avoids virtual dispatch for begin-offset access.
  using BeginsSet = OrderedSet<ElementType,
                               ElementType,
                               BeginsSizePolicy,
                               BeginsOffsetPolicy,
                               BeginsStridePolicy,
                               IndirectionPolicy,
                               NoSubset,
                               ConcreteInterface>;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType = BeginsSizePolicy;

  using IndirectionBufferType = typename IndirectionPolicy::IndirectionBufferType;
  using IndirectionPtrType = typename IndirectionPolicy::IndirectionPtrType;

  VariableCardinality() : m_begins() { }
  VariableCardinality(BeginsSet begins) : m_begins(begins) { }
  VariableCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionPtrType data)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize + 1).data(data);
  }

  AXOM_HOST_DEVICE ElementType size(ElementType fromPos) const
  {
    return offset(fromPos + 1) - offset(fromPos);
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  /*!
   * \brief Returns the from-set position owning \a relationOffset, or -1 if none does.
   *
   * \pre Begin offsets start at zero and are nondecreasing, with one offset per
   *      from-set element followed by the total number of relation indices.
   * \note Negative positions and positions at or beyond totalSize() return -1.
   *
   * \note Binary search takes O(log(fromSetSize)) time. It finds the first
   *  from-set position i for which offset(i+1) > relationOffset.
   *  MappedVariableCardinality provides O(1) lookup using an additional buffer
   *  of totalSize() from-set positions.
   */
  AXOM_HOST_DEVICE ElementType firstIndex(ElementType relationOffset) const
  {
    const ElementType numRows = m_begins.size() - 1;
    if(numRows <= ElementType {} || relationOffset < ElementType {} ||
       relationOffset < offset(ElementType {}) || offset(numRows) <= relationOffset)
    {
      return ElementType(-1);
    }

    ElementType lo = ElementType {};
    ElementType hi = numRows - 1;
    while(lo < hi)
    {
      const ElementType mid = lo + (hi - lo) / 2;
      if(offset(mid + 1) > relationOffset)
      {
        hi = mid;
      }
      else
      {
        lo = mid + 1;
      }
    }
    return lo;
  }

  IndirectionPtrType offsetData() { return m_begins.data(); }

  const IndirectionPtrType offsetData() const { return m_begins.data(); }

  ElementType totalSize() const
  {
    return m_begins.empty() ? ElementType() : offset(m_begins.size() - 1);
  }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    return detail::hasValidBeginOffsets(m_begins, fromSet, verboseOutput);
  }

  BeginsSet m_begins;
};

/*!
 * \class MappedVariableCardinality
 * \brief Represents a mapping between two sets, where each element in the
 *  first set maps to an arbitrary number of elements in the second set.
 *
 * Uses the same begin-offset contract as VariableCardinality and an additional
 * buffer of totalSize() from-set positions. firstIndex() reads that buffer in O(1).
 * The caller supplies both buffers and keeps them valid while the policy is used.
 *
 * \tparam ElementType The flat position type used for counts, offsets and stored from-set positions.
 * \tparam IndirectionPolicy How begin offsets and from-set positions are accessed.
 */
template <typename ElementType = int, typename IndirectionPolicy = ArrayIndirection<ElementType, ElementType>>
struct MappedVariableCardinality
{
  using BeginsSizePolicy = RuntimeSize<ElementType>;
  using BeginsOffsetPolicy = ZeroOffset<ElementType>;
  using BeginsStridePolicy = StrideOne<ElementType>;
  using BeginsIndirectionPolicy = IndirectionPolicy;

  // runtime size (fromSet.size()), striding from template parameter, no offset.
  // The concrete interface avoids virtual dispatch for index access.
  using IndexSet = OrderedSet<ElementType,
                              ElementType,
                              BeginsSizePolicy,
                              BeginsOffsetPolicy,
                              BeginsStridePolicy,
                              IndirectionPolicy,
                              NoSubset,
                              ConcreteInterface>;
  using BeginsSet = IndexSet;

  // The cardinality of each relational operator is determined by the
  // StridePolicy of the relation
  using RelationalOperatorSizeType = BeginsSizePolicy;

  using IndirectionBufferType = typename IndirectionPolicy::IndirectionBufferType;
  using IndirectionPtrType = typename IndirectionPolicy::IndirectionPtrType;

  MappedVariableCardinality() : m_begins() { }
  MappedVariableCardinality(BeginsSet begins) : m_begins(begins) { }
  MappedVariableCardinality(ElementType fromSetSize, typename BeginsSet::SetBuilder& builder)
  {
    builder.size(fromSetSize + 1);
    m_begins = builder;
  }

  void bindBeginOffsets(ElementType fromSetSize, IndirectionPtrType data)
  {
    m_begins = typename BeginsSet::SetBuilder().size(fromSetSize + 1).data(data);
  }

  /// \brief Bind the flat-to-from-set lookup buffer, optionally filling it from the begin offsets.
  /// \pre relationSize == totalSize(). Filling requires host-accessible buffers.
  void bindFirstIndices(ElementType relationSize, IndirectionPtrType data, bool fillIndices = true)
  {
    m_firstIndexes = typename IndexSet::SetBuilder().size(relationSize).data(data);
    if(fillIndices)
    {
      // Construct the flat-to-first mapping.
      for(ElementType fromIdx = ElementType {}; fromIdx < m_begins.size() - 1; ++fromIdx)
      {
        const ElementType beginIdx = offset(fromIdx);
        for(ElementType slotIdx = ElementType {}; slotIdx < size(fromIdx); ++slotIdx)
        {
          m_firstIndexes[slotIdx + beginIdx] = fromIdx;
        }
      }
    }
  }

  AXOM_HOST_DEVICE ElementType size(ElementType fromPos) const
  {
    return offset(fromPos + 1) - offset(fromPos);
  }

  AXOM_HOST_DEVICE ElementType offset(ElementType fromPos) const { return m_begins[fromPos]; }

  AXOM_HOST_DEVICE ElementType firstIndex(ElementType relationOffset) const
  {
    if(relationOffset < ElementType {} || relationOffset >= totalSize() ||
       relationOffset >= m_firstIndexes.size())
    {
      return ElementType(-1);
    }
    return m_firstIndexes[relationOffset];
  }

  IndirectionPtrType offsetData() { return m_begins.data(); }

  const IndirectionPtrType offsetData() const { return m_begins.data(); }

  ElementType totalSize() const
  {
    return m_begins.empty() ? ElementType() : offset(m_begins.size() - 1);
  }

  template <typename FromSetType>
  bool isValid(const FromSetType* fromSet, bool verboseOutput = false) const
  {
    if(!detail::hasValidBeginOffsets(m_begins, fromSet, verboseOutput) ||
       m_firstIndexes.size() != totalSize() || !m_firstIndexes.isValid(verboseOutput))
    {
      return false;
    }

    for(ElementType fromIdx = ElementType {}; fromIdx < fromSet->size(); ++fromIdx)
    {
      for(ElementType flatIdx = offset(fromIdx); flatIdx < offset(fromIdx + 1); ++flatIdx)
      {
        if(m_firstIndexes[flatIdx] != fromIdx)
        {
          return false;
        }
      }
    }

    return true;
  }

  IndexSet m_firstIndexes;
  BeginsSet m_begins;
};

}  // end namespace axom::slam::policies
