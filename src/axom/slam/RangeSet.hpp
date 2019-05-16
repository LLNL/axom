// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file RangeSet.hpp
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef SLAM_RANGE_SET_H_
#define SLAM_RANGE_SET_H_

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{


/**
 * \class PositionSet
 * \brief Alias template for an OrderedSet whose elements belong
 * to a contiguous range \f$ \in [0,size) \f$
 *
 * \tparam P The PositionType
 * \tparam E The ElementType
 * \sa OrderedSet
 */
template<typename P = slam::PositionType, typename E = slam::ElementType>
class PositionSet : public OrderedSet<P,E>
{

  using OrderedSetType =  OrderedSet<P,E> ;

  static const PositionType DEFAULT_SIZE =
    OrderedSetType::SizePolicyType::DEFAULT_VALUE;

  static const PositionType DEFAULT_OFFSET =
    OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;

  static const PositionType DEFAULT_STRIDE =
    OrderedSetType::StridePolicyType::DEFAULT_VALUE;

public:
  using PositionType = typename OrderedSetType::PositionType;
  using ElementType = typename OrderedSetType::ElementType ;


public:
  PositionSet(PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}

  PositionSet(const typename OrderedSetType::SetBuilder & builder)
    : OrderedSetType(builder) {}
};

/**
 * \class GenericRangeSet
 * \brief Models a set whose elements belong to a contiguous range
 *  \f$ \in [lowerIndex,upperIndex) \f$
 *
 * \note The \a ElementType here needs to be computable as offsets
 *  (of \a PositionType) from the lowerIndex.
 *  Examples include: signed and unsigned integral types.
 *  This version of a range set still allows you to have
 *  different policies on striding, indirection and subsetting
 *  \sa OrderedSet, PositionSet, RangeSet
 */
template<
  typename P = slam::PositionType,
  typename E = slam::ElementType,
  typename StridingPolicy = policies::StrideOne<P>,
  typename IndirectionPolicy = policies::NoIndirection<P, E>,
  typename SubsettingPolicy = policies::NoSubset >
class GenericRangeSet : public OrderedSet<
    P,
    E,
    policies::RuntimeSize<P>,
    policies::RuntimeOffset<P>,
    StridingPolicy,
    IndirectionPolicy,
    SubsettingPolicy >
{
public:
  using PositionType = P;
  using ElementType = E;

private:
  using OrderedSetType =
          OrderedSet<P,E,
                     policies::RuntimeSize<P>,
                     policies::RuntimeOffset<P>,
                     StridingPolicy,
                     IndirectionPolicy,
                     SubsettingPolicy>;

public:
  GenericRangeSet(
    PositionType size = OrderedSetType::SizePolicyType::DEFAULT_VALUE)
    : OrderedSetType(size,
                     OrderedSetType::OffsetPolicyType::DEFAULT_VALUE,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE) {}

  GenericRangeSet(const typename OrderedSetType::SetBuilder & builder)
    : OrderedSetType(builder) {}

  GenericRangeSet(PositionType lowerIndex, PositionType upperIndex)
    : OrderedSetType(upperIndex - lowerIndex,
                     lowerIndex,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE) {}
};


/**
 * \class RangeSet
 * \brief A specialization of GenericRangeSet with stride 1 and no indirection
 *
 * \tparam P The PositionType
 * \tparam E The ElementType
 *  \sa GenericRangeSet, OrderedSet, PositionSet
 *
 * \note The \a ElementType needs to be computable as offsets (of
 * \a PositionType) from the lowerIndex.
 * Examples include: signed and unsigned integral types
 */
template<typename P = slam::PositionType, typename E = slam::ElementType>
class RangeSet : public GenericRangeSet<P,E>
{

private:
  using GenericRangeSetType = GenericRangeSet<P,E>;

public:
  using PositionType = typename GenericRangeSetType::PositionType ;
  using ElementType = typename GenericRangeSetType::ElementType ;
private:
  static const PositionType DEFAULT_SIZE =
    GenericRangeSetType::SizePolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_OFFSET =
    GenericRangeSetType::OffsetPolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_STRIDE =
    GenericRangeSetType::StridePolicyType::DEFAULT_VALUE;


public:
  RangeSet(PositionType size = DEFAULT_SIZE)
    : GenericRangeSetType(size) {}

  RangeSet(PositionType lowerIndex, PositionType upperIndex)
    : GenericRangeSetType(lowerIndex, upperIndex) {}

  RangeSet(const typename GenericRangeSetType::SetBuilder & builder) 
    : GenericRangeSetType(builder) {}

};


} // end namespace slam
} // end namespace axom

#endif //  SLAM_RANGE_SET_H_
