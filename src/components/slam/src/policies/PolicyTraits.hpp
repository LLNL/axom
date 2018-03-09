/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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


#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_

/**
 * \file PolicyTraits.hpp
 *
 * A collection of utility traits classes for Slam policies.
 */

#include "slam/SizePolicies.hpp"
#include "slam/StridePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{

/**
 * \brief Definition of a type trait to adapt a StridePolicy into a SizePolicy
 */
template<typename StridePolicyType, typename IntType, int VAL = 1>
struct StrideToSize
{
  typedef CompileTimeSize<IntType, VAL> SizeType;
};

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
