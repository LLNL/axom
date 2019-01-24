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
 * \class
 *
 * \brief Concrete class (all template parameters fixed) for an array-based
 * indirection set
 */
class ArrayIndirectionSet : public OrderedSet<
    policies::RuntimeSize<Set::PositionType>,
    policies::ZeroOffset<Set::PositionType>,
    policies::StrideOne<Set::PositionType>,
    policies::ArrayIndirection<Set::PositionType, Set::ElementType> >
{
private:
  typedef OrderedSet<
      policies::RuntimeSize<Set::PositionType>,
      policies::ZeroOffset<Set::PositionType>,
      policies::StrideOne<Set::PositionType>,
      policies::ArrayIndirection<Set::PositionType, Set::ElementType>
      > OrderedSetType;

public:
  typedef OrderedSetType::PositionType PositionType;
  typedef OrderedSetType::IndexType IndexType;
  typedef OrderedSetType::ElementType ElementType;

  typedef OrderedSetType::SetBuilder ArraySetBuilder;

private:
  static const PositionType DEFAULT_SIZE =
    OrderedSetType::SizePolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_OFFSET =
    OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_STRIDE =
    OrderedSetType::StridePolicyType::DEFAULT_VALUE;

public:
  ArrayIndirectionSet (PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}

  ArrayIndirectionSet( const ArraySetBuilder& builder)
    : OrderedSetType(builder){}

  ~ArrayIndirectionSet () {}
};

/**
 * \class
 *
 * \brief Concrete class (all template parameters fixed)
 * for an STL vector-based indirection set
 */
class VectorIndirectionSet : public OrderedSet<
    policies::RuntimeSize<Set::PositionType>,
    policies::ZeroOffset<Set::PositionType>,
    policies::StrideOne<Set::PositionType>,
    policies::STLVectorIndirection<Set::PositionType, Set::ElementType> >
  // add parent subset ?
{
private:
  typedef OrderedSet<
      policies::RuntimeSize<Set::PositionType>,
      policies::ZeroOffset<Set::PositionType>,
      policies::StrideOne<Set::PositionType>,
      policies::STLVectorIndirection<Set::PositionType, Set::ElementType>
      > OrderedSetType;

  typedef OrderedSet::IndirectionPolicyType IndirectionPolicyType;

public:
  typedef OrderedSetType::PositionType PositionType;
  typedef OrderedSetType::IndexType IndexType;
  typedef OrderedSetType::ElementType ElementType;

  typedef IndirectionPolicyType::VectorType ArrType;
  typedef OrderedSetType::SetBuilder VectorSetBuilder;

private:
  static const PositionType DEFAULT_SIZE =
    OrderedSetType::SizePolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_OFFSET =
    OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;
  static const PositionType DEFAULT_STRIDE =
    OrderedSetType::StridePolicyType::DEFAULT_VALUE;

public:
  VectorIndirectionSet (PositionType size = DEFAULT_SIZE)
    : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}
  VectorIndirectionSet( const VectorSetBuilder& builder)
    : OrderedSetType(builder){}

  ~VectorIndirectionSet () {}
};


} // end namespace slam
} // end namespace axom

#endif //  SLAM_INDIRECTION_SET_H_
