// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_

/**
 * \file PolicyTraits.hpp
 *
 * A collection of utility traits classes for Slam policies.
 */

#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \brief Definition of a type trait to adapt a StridePolicy into a SizePolicy
 */
template <typename StridePolicyType, typename IntType, int VAL = 1>
struct StrideToSize
{
  using SizeType = CompileTimeSize<IntType, VAL>;
};

/**
 * \brief Specialization of StrideToSize trait for a RuntimeStride
 */
template <typename IntType>
struct StrideToSize<RuntimeStride<IntType>, IntType>
{
  using SizeType = RuntimeSize<IntType>;
};

/**
 * \brief Specialization of StrideToSize trait for a CompileTimeStride
 */
template <typename IntType, int VAL>
struct StrideToSize<CompileTimeStride<IntType, IntType(VAL)>, IntType, VAL>
{
  using SizeType = CompileTimeSize<IntType, IntType(VAL)>;
};

/**
 * \brief Specialization of StrideToSize trait for a StrideOne type
 */
template <typename IntType>
struct StrideToSize<StrideOne<IntType>, IntType>
{
  using SizeType = CompileTimeSize<IntType, StrideOne<IntType>::DEFAULT_VALUE>;
};

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICY_TRAITS_H_
