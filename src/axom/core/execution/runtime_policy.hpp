// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_
#define AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_

#include "axom/config.hpp"   /* for compile time defs. */
#include "axom/fmt/format.h" /* for axom::fmt */

#include <map>

// Helper preprocessor defines for using OPENMP, CUDA, and HIP policies.
// RAJA is required for using OpenMP, CUDA and HIP.
// UMPIRE is required for using CUDA and HIP.
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
    #define AXOM_RUNTIME_POLICY_USE_OPENMP
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
    #define AXOM_RUNTIME_POLICY_USE_CUDA
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
    #define AXOM_RUNTIME_POLICY_USE_HIP
  #endif
#endif

namespace axom
{
namespace core
{
namespace runtime_policy
{
/// Execution policies supported by Axom's configuration.
enum class Policy
{
  seq = 0
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  ,
  omp = 1
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  ,
  cuda = 2
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  ,
  hip = 3
#endif
};

//! @brief Mapping from policy name to policy enum.
// clang-format off
  static const std::map<std::string, Policy> s_nameToPolicy
  {
    {"seq", Policy::seq}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    , {"omp", Policy::omp}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    , {"cuda", Policy::cuda}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    , {"hip", Policy::hip}
#endif
  };

  //! @brief Mapping from policy enum to policy name.
  static const std::map<Policy, std::string> s_policyToName
  {
    {Policy::seq, "seq"}
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    , {Policy::omp, "omp"}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
    , {Policy::cuda, "cuda"}
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
    , {Policy::hip, "hip"}
#endif
  };
// clang-format on

/// Utility function to allow formating a Policy
static inline auto format_as(Policy pol) { return axom::fmt::underlying(pol); }

}  // end namespace runtime_policy
}  // end namespace core
}  // end namespace axom

#endif /* AXOM_CORE_EXECUTION_RUNTIME_POLICY_HPP_ */
