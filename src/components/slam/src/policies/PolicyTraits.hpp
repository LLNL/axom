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
  template<typename StridePolicyType, typename IntType, int VAL = 1> struct StrideToSize;

  /**
   * \brief Specialization of StrideToSize trait for a RuntimeStride
   */
  template<typename IntType>
  struct StrideToSize < RuntimeStride<IntType>, IntType >
  {
    typedef RuntimeSize<IntType> SizeType;
  };

  /**
   * \brief Specialization of StrideToSize trait for a CompileTimeStride
   */
  template<typename IntType, int VAL>
  struct StrideToSize< CompileTimeStride<IntType, IntType(VAL)>, IntType, VAL >
  {
    typedef CompileTimeSize<IntType, IntType(VAL)> SizeType;
  };

  /**
   * \brief Specialization of StrideToSize trait for a StrideOne type
   */
  template<typename IntType>
  struct StrideToSize< StrideOne<IntType>, IntType >
  {
    typedef CompileTimeSize<IntType, StrideOne<IntType>::DEFAULT_VALUE > SizeType;
  };


} // end namespace policies
} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICY_TRAITS_H_
