/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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
 * \file OrderedSet.hpp
 *
 * \brief Basic API for an ordered set of entities in a simulation
 * \note We are actually storing (ordered) multisets, since elements can be
 *  repeated an arbitrary number of times (e.g. for indirection sets)
 *
 */

#ifndef SLAM_ORDERED_SET_H_
#define SLAM_ORDERED_SET_H_

#include "axom/config.hpp"   // for AXOM_USE_CXX11


#include "axom/core/Types.hpp" // for nullptr
#include "axom/slic/interface/slic.hpp"

#include "axom/slam/Set.hpp"

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/OffsetPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/SubsettingPolicies.hpp"
#include "axom/slam/ModularInt.hpp"

#include "axom/slam/IteratorBase.hpp"

#ifdef AXOM_USE_CXX11
    #include <type_traits>
#endif

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
template<
  typename SizePolicy          = policies::RuntimeSize<Set::PositionType>,
  typename OffsetPolicy        = policies::ZeroOffset<Set::PositionType>,
  typename StridePolicy        = policies::StrideOne<Set::PositionType>,
  typename IndirectionPolicy   =
    policies::NoIndirection<Set::PositionType, Set::ElementType>,
  typename SubsettingPolicy    = policies::NoSubset >
struct OrderedSet : public Set,
                           SizePolicy,
                           OffsetPolicy,
                           StridePolicy,
                           IndirectionPolicy,
                           SubsettingPolicy
{
public:

  typedef Set::IndexType IndexType;
  typedef Set::PositionType PositionType;
  typedef IndexType ElementType;

  typedef SizePolicy SizePolicyType;
  typedef OffsetPolicy OffsetPolicyType;
  typedef StridePolicy StridePolicyType;
  typedef IndirectionPolicy IndirectionPolicyType;
  typedef SubsettingPolicy SubsettingPolicyType;

  typedef ModularInt<SizePolicy>  ModularIntType;

  struct SetBuilder;

#ifdef AXOM_USE_CXX11
  template<typename OrderedSetType> class OrderedSetIterator;

  typedef OrderedSetIterator<const OrderedSet>      const_iterator;
  typedef std::pair<const_iterator,const_iterator>  const_iterator_pair;

  typedef const_iterator iterator;
  typedef const_iterator_pair iterator_pair;
#endif // AXOM_USE_CXX11

public:
  OrderedSet(PositionType size    = SizePolicyType::DEFAULT_VALUE,
             PositionType offset  = OffsetPolicyType::DEFAULT_VALUE,
             PositionType stride  = StridePolicyType::DEFAULT_VALUE
                                    // Note: constructor does not yet take an
                                    // indirection type pointer...
                                    // const Set* parentSet = &s_nullSet
             )
    : SizePolicyType(size),
    OffsetPolicyType(offset),
    StridePolicyType(stride)
    //, SubsettingPolicyType(parentSet)
  {}

  OrderedSet(const SetBuilder & builder)
    : SizePolicyType(builder.m_size),
    OffsetPolicyType(builder.m_offset),
    StridePolicyType(builder.m_stride),
    IndirectionPolicyType(builder.m_data),
    SubsettingPolicyType(builder.m_parent)
  {}

  OrderedSet(const OrderedSet& oset)
    : SizePolicyType(oset),
    OffsetPolicyType(oset),
    StridePolicyType(oset),
    IndirectionPolicyType(oset),
    SubsettingPolicyType(oset)
  {}



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

    typedef typename IndirectionPolicyType::IndirectionBufferType DataType;
    typedef typename SubsettingPolicyType::ParentSetType ParentSetType;

    SetBuilder& size(PositionType sz)
    {
      m_size    = SizePolicyType(sz);
      return *this;
    }

    SetBuilder& offset(PositionType off)
    {
      m_offset  = OffsetPolicyType(off);
      return *this;
    }

    SetBuilder& stride(PositionType str)
    {
      m_stride  = StridePolicyType(str);
      return *this;
    }

    SetBuilder& data(DataType* bufPtr)
    {
      m_data   = IndirectionPolicyType(bufPtr);
      return *this;
    }

    SetBuilder& parent(ParentSetType* parSet)
    {
      m_parent = SubsettingPolicyType(parSet);
      return *this;
    }

    /** Alternate means of setting the offset and size from a contiguous range
       of values */
    SetBuilder& range(PositionType lower, PositionType upper)
    {
      // Set by range rather than size and offset.
      // Question: Should we ensure that only one of these options is called
      //   (e.g. size [+offset] or range, but not both)?
      m_offset = OffsetPolicyType(lower);
      m_size = SizePolicyType(upper - lower);
      return *this;
    }

private:
    SizePolicyType m_size;
    OffsetPolicyType m_offset;
    StridePolicyType m_stride;
    IndirectionPolicyType m_data;
    SubsettingPolicyType m_parent;
  };

#ifdef AXOM_USE_CXX11
  /**
   * \class OrderedSetIterator
   * \brief An iterator type for an ordered set
   *
   * Uses the set's policies for efficient iteration
   */
  template<typename OrderedSet>
  class OrderedSetIterator : public IteratorBase<
      OrderedSetIterator<OrderedSet>,
      typename OrderedSet::ElementType>
  {
public:

    typedef OrderedSetIterator<OrderedSet>              iter;
    typedef typename OrderedSet::ElementType ElementType;
    typedef typename OrderedSet::PositionType PositionType;

    typedef typename OrderedSet::IndirectionPolicyType IndirectionType;
    typedef typename OrderedSet::StridePolicyType StrideType;

    using IterBase = IteratorBase<OrderedSetIterator<OrderedSet>,
                                  typename OrderedSet::ElementType>;
    using IterBase::m_pos;

public:

    //OrderedSetIterator(PositionType pos) : IterBase(pos) {}
    OrderedSetIterator(PositionType pos, const OrderedSet& oSet)
      : IterBase(pos), m_orderedSet(oSet) {}


    const ElementType & operator*()    const {
      // Note: Since we return a reference to the pointed-to value, we need
      // different functions
      //       for OrderedSets with indirection buffers than with those that
      // have no indirection
      typedef policies::
        NoIndirection<PositionType,ElementType> NoIndirectionType;
      return indirection(
        HasIndirection<!std::is_same<IndirectionType,
                                     NoIndirectionType>::value >(), 0);
    }

    ElementType operator[](PositionType n) const
    {
      return *(this->operator+(n));
    }

protected:
    /** Implementation of advance() as required by IteratorBase */
    void advance(PositionType n) { m_pos += n * stride(); }

private:
    inline const PositionType stride() const
    {
      return m_orderedSet.StrideType::stride();
    }

    template<bool> class HasIndirection {};

    template<typename T>
    inline const ElementType& indirection(HasIndirection<true>, T) const
    {
      return m_orderedSet.IndirectionType::indirection(m_pos);
    }

    template<typename T>
    inline const ElementType& indirection(HasIndirection<false>, T) const
    {
      return m_pos;
    }

private:
    OrderedSet m_orderedSet;
  };

public:     // Functions related to iteration

  const_iterator      begin() const
  {
    return const_iterator( OffsetPolicyType::offset(), *this);
  }

  const_iterator      end()   const
  {
    return const_iterator( SizePolicyType::size() * StridePolicyType::stride()
                           + OffsetPolicyType::offset(), *this);
  }
  const_iterator_pair range() const { return std::make_pair(begin(), end()); }
#endif // AXOM_USE_CXX11

public:
  /**
   * \brief Given a position in the Set, return a position in the larger index
   *  space
   */
  inline typename IndirectionPolicy::IndirectionResult
  operator[](PositionType pos) const
  {
    verifyPositionImpl(pos);
    return IndirectionPolicy::indirection( pos * StridePolicyType::stride()
                                           + OffsetPolicyType::offset() );
  }

  inline ElementType  at(PositionType pos)         const
  {
    return operator[](pos);
  }


  inline PositionType size()  const { return SizePolicyType::size(); }
  inline bool         empty() const { return SizePolicyType::empty(); }

  bool isValid(bool verboseOutput = false) const;

  bool isSubset() const { return SubsettingPolicy::isSubset(); }

  /**
   * \brief checks whether the given position (index) is valid.
   *
   * An index pos is valid when \f$ 0 \le pos < size() \f$
   * \return true if the position is valid, false otherwise
   */
  bool                isValidIndex(PositionType pos) const
  {
    return pos >= 0 && pos < size();
  }

private:

  inline void verifyPosition(PositionType pos)       const
  {
    verifyPositionImpl(pos);
  }
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos))       const
  {
    SLIC_ASSERT_MSG(
      isValidIndex(pos),
      "SLAM::OrderedSet -- requested out-of-range element at position "
      << pos << ", but set only has " << size() << " elements." );
  }

private:
  /// NOTE: All data for OrderedSet is associated with parent policy classes
};

template<
  typename SizePolicy,
  typename OffsetPolicy,
  typename StridePolicy,
  typename IndirectionPolicy,
  typename SubsettingPolicy >
bool OrderedSet<SizePolicy,OffsetPolicy, StridePolicy, IndirectionPolicy,
                SubsettingPolicy>::isValid(bool verboseOutput) const
{
  bool bValid =  SizePolicyType::isValid(verboseOutput)
                && OffsetPolicyType::isValid(verboseOutput)
                && StridePolicyType::isValid(verboseOutput)
                && IndirectionPolicyType::isValid(size(),
                                                  OffsetPolicy::offset(),
                                                  StridePolicy::stride(),
                                                  verboseOutput)
#ifdef AXOM_USE_CXX11
                && SubsettingPolicyType::isValid(begin(), end(), verboseOutput)
#endif
  ;

  return bValid;
}


} // end namespace slam
} // end namespace axom

#endif //  SLAM_ORDERED_SET_H_
