/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_

/**
 * \file PolicyTraits.hpp
 *
 * A collection of utility traits classes for Slam policies.
 */

#include "slam/SizePolicies.hpp"
#include "slam/StridePolicies.hpp"

namespace axom {
namespace slam {
namespace policies {

  /**
   * \brief Definition of a type trait to adapt a StridePolicy into a SizePolicy
   */
  template<typename StridePolicyType, typename IntType, IntType VAL> struct StrideToSize;

  /**
   * \brief Specialization of StrideToSize trait for a RuntimeStride
   */
  template<>
  struct StrideToSize<
    RuntimeStride< Set::PositionType >,
    Set::PositionType,
    RuntimeStride< Set::PositionType >::DEFAULT_VALUE >
  {
    typedef RuntimeSize<typename Set::PositionType> SizeType;
  };

  /**
   * \brief Specialization of StrideToSize trait for a CompileTimeStride
   */
  template<Set::PositionType VAL>
  struct StrideToSize< CompileTimeStride<Set::PositionType, VAL>, Set::PositionType, VAL >
  {
    typedef CompileTimeSize<Set::PositionType, VAL> SizeType;
  };

  /**
   * \brief Specialization of StrideToSize trait for a StrideOne type
   */
  template<>
  struct StrideToSize< StrideOne<Set::PositionType>, Set::PositionType,  StrideOne<Set::PositionType>::DEFAULT_VALUE >
  {
    typedef CompileTimeSize<Set::PositionType, StrideOne<Set::PositionType  >::DEFAULT_VALUE > SizeType;
  };


} // end namespace policies
} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICY_TRAITS_H_
