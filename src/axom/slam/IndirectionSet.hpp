// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IndirectionSet.hpp
 *
 * \brief A few indirection set classes
 */

#ifndef SLAM_INDIRECTION_SET_H_
#define SLAM_INDIRECTION_SET_H_

#include <cstddef>
#include <vector>

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{


/**
 * \class ArrayIndirectionSet
 *
 * \brief Concrete class (all template parameters fixed) for an array-based
 * indirection set
 */
template<
  typename PosType = slam::PositionType,
  typename ElemType = slam::ElementType >
class ArrayIndirectionSet : public OrderedSet<
    PosType,
    ElemType,
    policies::RuntimeSize<PosType>,
    policies::ZeroOffset<PosType>,
    policies::StrideOne<PosType>,
    policies::ArrayIndirection<PosType, ElemType> >
{
private:
  using OrderedSetType =
          OrderedSet<
            PosType,
            ElemType,
            policies::RuntimeSize<PosType>,
            policies::ZeroOffset<PosType>,
            policies::StrideOne<PosType>,
            policies::ArrayIndirection<PosType, ElemType> >;

public:
  using PositionType = PosType;
  using ElementType = ElemType;

  using ArraySetBuilder = typename OrderedSetType::SetBuilder;

private:
  static constexpr PositionType DEFAULT_SIZE =
    OrderedSetType::SizePolicyType::DEFAULT_VALUE;

  static constexpr PositionType DEFAULT_OFFSET =
    OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;

  static constexpr PositionType DEFAULT_STRIDE =
    OrderedSetType::StridePolicyType::DEFAULT_VALUE;

public:
  ArrayIndirectionSet (PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}

  ArrayIndirectionSet( const ArraySetBuilder& builder)
    : OrderedSetType(builder){}

  ~ArrayIndirectionSet () = default;
};

/**
 * \class VectorIndirectionSet
 *
 * \brief Concrete class (all template parameters fixed)
 * for an STL vector-based indirection set
 */
template<
  typename PosType = slam::PositionType,
  typename ElemType = slam::ElementType >
class VectorIndirectionSet : public OrderedSet<
    PosType,
    ElemType,
    policies::RuntimeSize<PosType>,
    policies::ZeroOffset<PosType>,
    policies::StrideOne<PosType>,
    policies::STLVectorIndirection<PosType, ElemType> >
  // add parent subset ?
{
private:
  using OrderedSetType =
          OrderedSet<
            PosType,
            ElemType,
            policies::RuntimeSize<PosType>,
            policies::ZeroOffset<PosType>,
            policies::StrideOne<PosType>,
            policies::STLVectorIndirection<PosType, ElemType>
            >;

  using IndirectionPolicyType = typename OrderedSetType::IndirectionPolicyType;

public:
  using PositionType = PosType;
  using ElementType = ElemType;

  using ArrType = typename IndirectionPolicyType::VectorType;
  using VectorSetBuilder = typename OrderedSetType::SetBuilder;

private:
  static constexpr PositionType DEFAULT_SIZE =
    OrderedSetType::SizePolicyType::DEFAULT_VALUE;

  static constexpr PositionType DEFAULT_OFFSET =
    OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;

  static constexpr PositionType DEFAULT_STRIDE =
    OrderedSetType::StridePolicyType::DEFAULT_VALUE;

public:
  VectorIndirectionSet (PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}

  VectorIndirectionSet( const VectorSetBuilder& builder)
    : OrderedSetType(builder){}

  ~VectorIndirectionSet () = default;
};


} // end namespace slam
} // end namespace axom

#endif //  SLAM_INDIRECTION_SET_H_
