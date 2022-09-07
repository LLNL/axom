// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file OrderedSet.hpp
 *
 * \brief Basic API for an ordered set of entities in a simulation
 * \note We are actually storing (ordered) multisets, since elements can be
 *  repeated an arbitrary number of times (e.g. for indirection sets)
 *
 */

#ifndef SLAM_ORDERED_SET_H_
#define SLAM_ORDERED_SET_H_

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic.hpp"

#include "axom/slam/Set.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/SubsettingPolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/core/IteratorBase.hpp"

#include <type_traits>
#include <cstddef>
#include <vector>
#include <iterator>

namespace axom
{
namespace slam
{
/**
 * \class OrderedSet
 *
 * \brief Models a set whose elements can be defined as strided offsets
 * of the position, possibly with a level of indirection.
 *
 * In an OrderedSet, the element at position pos can be defined as:
 *     static_cast<ElementType>( indirection[ pos * stride + offset ] )
 */
template <typename PosType = slam::DefaultPositionType,
          typename ElemType = slam::DefaultElementType,
          typename SizePolicy = policies::RuntimeSize<PosType>,
          typename OffsetPolicy = policies::ZeroOffset<PosType>,
          typename StridePolicy = policies::StrideOne<PosType>,
          typename IndirectionPolicy = policies::NoIndirection<PosType, ElemType>,
          typename SubsettingPolicy = policies::NoSubset>
struct OrderedSet : SizePolicy,
                    OffsetPolicy,
                    StridePolicy,
                    IndirectionPolicy,
                    SubsettingPolicy
{
public:
  using PositionType = PosType;
  using ElementType = ElemType;

  using SizePolicyType = SizePolicy;
  using OffsetPolicyType = OffsetPolicy;
  using StridePolicyType = StridePolicy;
  using IndirectionPolicyType = IndirectionPolicy;
  using SubsettingPolicyType = SubsettingPolicy;

  using ModularIntType = ModularInt<SizePolicy>;

  using PositionSet = OrderedSet<PositionType>;

  struct SetBuilder;

  // types for OrderedSet iterator
  template <typename T, bool C>
  class OrderedSetIterator;

  using const_iterator = OrderedSetIterator<ElementType, true>;
  using const_iterator_pair = std::pair<const_iterator, const_iterator>;

  using iterator = OrderedSetIterator<ElementType, false>;
  using iterator_pair = std::pair<iterator, iterator>;

public:
  OrderedSet(PositionType size = SizePolicyType::DEFAULT_VALUE,
             PositionType offset = OffsetPolicyType::DEFAULT_VALUE,
             PositionType stride = StridePolicyType::DEFAULT_VALUE
             // Note: constructor does not yet take an
             // indirection type pointer...
             // const Set* parentSet = &s_nullSet
             )
    : SizePolicyType(size)
    , OffsetPolicyType(offset)
    , StridePolicyType(stride)
  //, SubsettingPolicyType(parentSet)
  { }

  OrderedSet(const SetBuilder& builder)
    : SizePolicyType(builder.m_size)
    , OffsetPolicyType(builder.m_offset)
    , StridePolicyType(builder.m_stride)
    , IndirectionPolicyType(builder.m_data)
    , SubsettingPolicyType(builder.m_parent)
  { }

  template <typename OtherIndirectionPolicy>
  OrderedSet(const OrderedSet<PosType,
                              ElemType,
                              SizePolicyType,
                              OffsetPolicyType,
                              StridePolicyType,
                              OtherIndirectionPolicy,
                              SubsettingPolicyType>& other)
    : SizePolicyType(other)
    , OffsetPolicyType(other)
    , StridePolicyType(other)
    , IndirectionPolicyType(other)
    , SubsettingPolicyType(other)
  { }

  OrderedSet(const OrderedSet& other) = default;
  OrderedSet& operator=(const OrderedSet& other) = default;

  OrderedSet(OrderedSet&& other) = default;
  OrderedSet& operator=(OrderedSet&& other) = default;

  ~OrderedSet() = default;

public:
  template <typename OtherSet>
  inline bool operator==(const OtherSet& other) const
  {
    using OtherPosType = typename OtherSet::PositionType;
    using OtherElemType = typename OtherSet::ElementType;

    using CmnElemType = typename std::common_type<ElemType, OtherElemType>::type;

    if(this->size() != static_cast<PosType>(other.size()))
    {
      return false;
    }

    for(auto pos = PosType {}; pos < this->size(); pos++)
    {
      auto&& e1 = static_cast<CmnElemType&&>(this->at(pos));
      auto&& e2 =
        static_cast<CmnElemType&&>(other.at(static_cast<OtherPosType>(pos)));
      if(e1 != e2)
      {
        return false;
      }
    }
    return true;
  }

  template <typename OtherSet>
  inline bool operator!=(const OtherSet& other) const
  {
    return !(this->operator==(other));
  }

public:
  /**
   * \class SetBuilder
   * \brief Helper class for constructing an ordered set.
   *
   *  Uses named parameter idiom to enable function chaining and for better code
   *  self-documentation
   * */
  struct SetBuilder
  {
    friend struct OrderedSet;

    using DataType = typename IndirectionPolicyType::IndirectionPtrType;
    using ParentSetType = typename SubsettingPolicyType::ParentSetType;

    SetBuilder& size(PositionType sz)
    {
      SLIC_ASSERT_MSG(
        !m_hasRange,
        "Cannot call both SetBuilder::size() and SetBuilder::range()");

      m_size = SizePolicyType(sz);
      return *this;
    }

    SetBuilder& offset(PositionType off)
    {
      SLIC_ASSERT(!m_hasRange || off == m_offset.offset());

      m_offset = OffsetPolicyType(off);

      return *this;
    }

    SetBuilder& stride(PositionType str)
    {
      m_stride = StridePolicyType(str);

      // Need to recompute size if range() function was used
      if(m_hasRange)
      {
        m_size = SizePolicyType(recomputeSize());
      }

      return *this;
    }

    SetBuilder& data(DataType bufPtr)
    {
      m_data = IndirectionPolicyType(bufPtr);
      return *this;
    }

    template <typename OtherParentType>
    SetBuilder& parent(OtherParentType* parSet)
    {
      m_parent = SubsettingPolicyType(parSet);
      return *this;
    }

    /// Alternate means of setting offset and size from contiguous range
    SetBuilder& range(PositionType lower, PositionType upper)
    {
      // Set by range rather than size and offset.
      m_offset = OffsetPolicyType(lower);

      m_hasRange = true;
      m_rangeLower = lower;
      m_rangeUpper = upper;

      // set size, account for a stride that has not yet been set
      m_size = SizePolicyType(recomputeSize());
      return *this;
    }

  private:
    /// Recomputes the size based on upper and lower range as well as stride
    PositionType recomputeSize() const
    {
      using axom::utilities::abs;
      using axom::utilities::ceil;

      if(m_hasRange)
      {
        const double str = m_stride.stride();
        const auto diff = (m_rangeUpper - m_rangeLower);

        // size is 0 if upper==lower, or signs of diff and stride differ
        return (diff == 0 || ((diff > 0) != (str > 0))) ? 0 : ceil(diff / str);
      }
      else
      {
        return m_size.size();
      }
    }

  private:
    SizePolicyType m_size;
    OffsetPolicyType m_offset;
    StridePolicyType m_stride;
    IndirectionPolicyType m_data;
    SubsettingPolicyType m_parent;

    // Some additional params to account for range() convenience function
    bool m_hasRange {false};
    PositionType m_rangeLower {0};  // only used when m_hasRange == true
    PositionType m_rangeUpper {0};  // only used when m_hasRange == true
  };

  /**
   * \class OrderedSetIterator
   * \brief An stl-compliant random iterator type for an ordered set
   *
   * Uses the set's policies for efficient iteration
   * \tparam T The result type of the iteration
   * \tparam Const Boolean to indicate if this is a const iterator
   *
   * \note Most operators are implemented via the \a IteratorBase class
   * \note The reference type is determined based on the IndirectionPolicy
   * and on the \a Const template parameter
   *
   * \note Use of a const template parameter with conditional member and pointer
   * operations based on ideas from https://stackoverflow.com/a/49425072
   */
  template <typename T, bool Const>
  class OrderedSetIterator
    : public IteratorBase<OrderedSetIterator<T, Const>, PositionType>
  {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = PositionType;

    using reference = typename std::conditional<
      Const,
      typename OrderedSet::IndirectionPolicyType::ConstIndirectionResult,
      typename OrderedSet::IndirectionPolicyType::IndirectionResult>::type;

    using pointer = typename std::conditional<Const, const T*, T*>::type;

    using IterBase = IteratorBase<OrderedSetIterator<T, Const>, PositionType>;

    using IndirectionType = typename OrderedSet::IndirectionPolicyType;
    using StrideType = typename OrderedSet::StridePolicyType;
    using OffsetType = typename OrderedSet::OffsetPolicyType;

    using IterBase::m_pos;

  public:
    /// \name Constructors, copying and assignment
    /// \{
    OrderedSetIterator() = default;

    OrderedSetIterator(PositionType pos, const OrderedSet& oSet)
      : IterBase(pos)
      , m_orderedSet(oSet)
    { }

    OrderedSetIterator(const OrderedSetIterator& it) = default;

    OrderedSetIterator& operator=(const OrderedSetIterator& it)
    {
      this->m_pos = it.m_pos;
      this->m_orderedSet = const_cast<OrderedSet&>(it.m_orderedSet);
      return *this;
    }
    /// \}

    /// \name Member and pointer operators
    /// \note We use the \a enable_if construct to implement both
    /// const and non-const iterators in the same implementation.
    /// \{

    /// Indirection operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, reference>::type operator*()
    {
      return m_orderedSet.IndirectionType::indirection(m_pos);
    }

    /// Indirection operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, reference>::type operator*() const
    {
      return m_orderedSet.IndirectionType::indirection(m_pos);
    }

    /// Structure dereference operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, pointer>::type operator->()
    {
      return &(m_orderedSet.IndirectionType::indirection(m_pos));
    }

    /// Structure dereference operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, pointer>::type operator->() const
    {
      return &(m_orderedSet.IndirectionType::indirection(m_pos));
    }

    /// Subscript operator for non-const iterator
    template <bool _Const = Const>
    typename std::enable_if<!_Const, reference>::type operator[](PositionType n)
    {
      return *(*this + n);
    }

    /// Subscript operator for const iterator
    template <bool _Const = Const>
    typename std::enable_if<_Const, reference>::type operator[](PositionType n) const
    {
      return *(*this + n);
    }

    /// \}

    /// \name Conversion operators
    /// \{

    /// Convert from iterator type to const_iterator type
    template <typename U>
    operator OrderedSetIterator<U, true>() const
    {
      return OrderedSetIterator<U, true>(this->m_pos, this->m_orderedSet);
    }
    /// \}

    PositionType index() const { return (m_pos - offset()) / stride(); }

  protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n * stride(); }

  private:
    inline const PositionType stride() const
    {
      return m_orderedSet.StrideType::stride();
    }

    inline const PositionType offset() const
    {
      return m_orderedSet.OffsetType::offset();
    }

  private:
    OrderedSet m_orderedSet;
  };

public:  // Functions related to iteration
  iterator begin() { return iterator(OffsetPolicyType::offset(), *this); }

  const_iterator begin() const
  {
    return const_iterator(OffsetPolicyType::offset(), *this);
  }

  iterator end()
  {
    return iterator(SizePolicyType::size() * StridePolicyType::stride() +
                      OffsetPolicyType::offset(),
                    *this);
  }

  const_iterator end() const
  {
    return const_iterator(SizePolicyType::size() * StridePolicyType::stride() +
                            OffsetPolicyType::offset(),
                          *this);
  }

  iterator_pair range() { return std::make_pair(begin(), end()); }

  const_iterator_pair range() const { return std::make_pair(begin(), end()); }

public:
  /**
   * \brief Given a position in the Set, return a position in the larger index
   *  space
   */
  inline typename IndirectionPolicy::ConstIndirectionResult operator[](
    PositionType pos) const
  {
    verifyPositionImpl(pos);
    return IndirectionPolicy::indirection(pos * StridePolicyType::stride() +
                                          OffsetPolicyType::offset());
  }

  inline typename IndirectionPolicy::IndirectionResult operator[](PositionType pos)
  {
    verifyPositionImpl(pos);
    return IndirectionPolicy::indirection(pos * StridePolicyType::stride() +
                                          OffsetPolicyType::offset());
  }
  inline ElementType at(PositionType pos) const { return operator[](pos); }

  inline PositionType size() const { return SizePolicyType::size(); }
  inline bool empty() const { return SizePolicyType::empty(); }

  bool isValid(bool verboseOutput = false) const;

  bool isSubset() const { return SubsettingPolicy::isSubset(); }

  /**
   * \brief checks whether the given position (index) is valid.
   *
   * An index pos is valid when \f$ 0 \le pos < size() \f$
   * \return true if the position is valid, false otherwise
   */
  inline bool isValidIndex(PositionType pos) const
  {
    return pos >= 0 && pos < size();
  }

  /**
   * \brief returns a PositionSet over the set's positions
   *
   * This can be used to simplify code to loop through the elements of a set.
   */
  PositionSet positions() const { return PositionSet(size()); }

private:
  inline void verifyPosition(PositionType pos) const
  {
    verifyPositionImpl(pos);
  }
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos)) const
  {
    SLIC_ASSERT_MSG(
      isValidIndex(pos),
      "SLAM::OrderedSet -- requested out-of-range element at position "
        << pos << ", but set only has " << size() << " elements.");
  }

private:
  /// NOTE: All data for OrderedSet is associated with parent policy classes
};

template <typename PosType,
          typename ElemType,
          typename SizePolicy,
          typename OffsetPolicy,
          typename StridePolicy,
          typename IndirectionPolicy,
          typename SubsettingPolicy>
bool OrderedSet<PosType,
                ElemType,
                SizePolicy,
                OffsetPolicy,
                StridePolicy,
                IndirectionPolicy,
                SubsettingPolicy>::isValid(bool verboseOutput) const
{
  bool bValid = SizePolicyType::isValid(verboseOutput) &&
    OffsetPolicyType::isValid(verboseOutput) &&
    StridePolicyType::isValid(verboseOutput) &&
    IndirectionPolicyType::isValid(size(),
                                   OffsetPolicy::offset(),
                                   StridePolicy::stride(),
                                   verboseOutput) &&
    SubsettingPolicyType::isValid(begin(), end(), verboseOutput);

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_ORDERED_SET_H_
