// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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
template <typename P = slam::DefaultPositionType, typename E = slam::DefaultElementType>
class PositionSet : public OrderedSet<P, E>
{
public:
  using PositionType = P;
  using ElementType = E;

private:
  using OrderedSetType = OrderedSet<P, E>;

  static const PositionType DEFAULT_SIZE;
  static const PositionType DEFAULT_OFFSET;
  static const PositionType DEFAULT_STRIDE;

public:
  PositionSet(PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE)
  { }

  PositionSet(const typename OrderedSetType::SetBuilder& builder)
    : OrderedSetType(builder)
  { }
};

template <typename P, typename E>
const P PositionSet<P, E>::DEFAULT_SIZE =
  PositionSet<P, E>::OrderedSetType::SizePolicyType::DEFAULT_VALUE;

template <typename P, typename E>
const P PositionSet<P, E>::DEFAULT_OFFSET =
  PositionSet<P, E>::OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;

template <typename P, typename E>
const P PositionSet<P, E>::DEFAULT_STRIDE =
  PositionSet<P, E>::OrderedSetType::StridePolicyType::DEFAULT_VALUE;

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
template <typename P = slam::DefaultPositionType,
          typename E = slam::DefaultElementType,
          typename StridingPolicy = policies::StrideOne<P>,
          typename IndirectionPolicy = policies::NoIndirection<P, E>,
          typename SubsettingPolicy = policies::NoSubset>
class GenericRangeSet : public OrderedSet<P,
                                          E,
                                          policies::RuntimeSize<P>,
                                          policies::RuntimeOffset<P>,
                                          StridingPolicy,
                                          IndirectionPolicy,
                                          SubsettingPolicy>
{
public:
  using PositionType = P;
  using ElementType = E;

private:
  using OrderedSetType = OrderedSet<P,
                                    E,
                                    policies::RuntimeSize<P>,
                                    policies::RuntimeOffset<P>,
                                    StridingPolicy,
                                    IndirectionPolicy,
                                    SubsettingPolicy>;

public:
  GenericRangeSet(PositionType size = OrderedSetType::SizePolicyType::DEFAULT_VALUE)
    : OrderedSetType(size,
                     OrderedSetType::OffsetPolicyType::DEFAULT_VALUE,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE)
  { }

  GenericRangeSet(const typename OrderedSetType::SetBuilder& builder)
    : OrderedSetType(builder)
  { }

  GenericRangeSet(PositionType lowerIndex, PositionType upperIndex)
    : OrderedSetType(upperIndex - lowerIndex,
                     lowerIndex,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE)
  { }
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
template <typename P = slam::DefaultPositionType, typename E = slam::DefaultElementType>
class RangeSet : public GenericRangeSet<P, E>
{
public:
  using PositionType = P;
  using ElementType = E;

private:
  using GenericRangeSetType = GenericRangeSet<P, E>;

private:
  static const PositionType DEFAULT_SIZE;
  static const PositionType DEFAULT_OFFSET;
  static const PositionType DEFAULT_STRIDE;

public:
  RangeSet(PositionType size = DEFAULT_SIZE) : GenericRangeSetType(size) { }

  RangeSet(PositionType lowerIndex, PositionType upperIndex)
    : GenericRangeSetType(lowerIndex, upperIndex)
  { }

  RangeSet(const typename GenericRangeSetType::SetBuilder& builder)
    : GenericRangeSetType(builder)
  { }
};

template <typename P, typename E>
const P RangeSet<P, E>::DEFAULT_SIZE =
  RangeSet<P, E>::GenericRangeSetType::SizePolicyType::DEFAULT_VALUE;

template <typename P, typename E>
const P RangeSet<P, E>::DEFAULT_OFFSET =
  RangeSet<P, E>::GenericRangeSetType::OffsetPolicyType::DEFAULT_VALUE;

template <typename P, typename E>
const P RangeSet<P, E>::DEFAULT_STRIDE =
  RangeSet<P, E>::GenericRangeSetType::StridePolicyType::DEFAULT_VALUE;

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_RANGE_SET_H_
