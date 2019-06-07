// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_PRIMAL_EXECUTIONSPACE_HPP_
#define AXOM_PRIMAL_EXECUTIONSPACE_HPP_

#include "axom/config.hpp"

// RAJA includes
#include "RAJA/RAJA.hpp"

// Umpire includes
#include "umpire/config.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/op/MemoryOperationRegistry.hpp"

#include "axom/core/memory_management.hpp"

/*!
 * \file
 *
 * \brief Defines the list of available execution spaces for spatial
 *  acceleration data-structures.
 *
 * \note The execution space is provided as a template argument to
 *
 */

namespace axom
{
namespace primal
{

/// \name Execution Spaces
/// @{

/*!
 * \brief Indicates parallel execution on the GPU with CUDA.
 */
template < int BLOCK_SIZE >
struct CUDA_EXEC{ };

/*!
 * \brief Indicates sequential execution on the CPU.
 */
struct SEQ_EXEC{ };

/*!
 * \brief Indicates parallel execution on the CPU using OpenMP.
 */
struct OMP_EXEC{ };

/// @}


/// \name Execution Space Traits
/// @{


/*!
 * \brief The execution_space is a traits class.
 *
 *  The execution_space traits class binds the execution space to:
 *  * corresponding RAJA execution policies and
 *  * the default memory allocator to use.
 *
 *  \tparam ExecSpace the execution space
 *
 */
template < typename ExecSpace >
struct execution_space
{
  using raja_exec   = void;
  using raja_reduce = void;
  using raja_atomic = void;

  static constexpr bool valid() noexcept { return false; };
  static constexpr char* name() noexcept { return (char*)"[UNDEFINED]"; };
  static constexpr int allocatorID() noexcept
  { return axom::INVALID_ALLOCATOR_ID; };
};

//------------------------------------------------------------------------------


#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA)

/*!
 * \brief execution_space traits specialization for CUDA_EXEC
 *
 * \tparam BLOCK_SIZE the number of CUDA threads to launch
 */
template < int BLOCK_SIZE >
struct execution_space< CUDA_EXEC< BLOCK_SIZE > >
{
  using raja_exec   = RAJA::cuda_exec< BLOCK_SIZE >;
  using raja_reduce = RAJA::cuda_reduce;
  using raja_atomic = RAJA::atomic::cuda_atomic;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[CUDA_EXEC]"; };
  static constexpr int allocatorID() noexcept
  { return umpire::resource::Unified; };
};

#endif

//------------------------------------------------------------------------------

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)

/*!
 * \brief execution_space traits specialization for OMP_EXEC
 */
template < >
struct execution_space< OMP_EXEC >
{
  using raja_exec   = RAJA::omp_parallel_for_exec;
  using raja_reduce = RAJA::omp_reduce;
  using raja_atomic = RAJA::atomic::omp_atomic;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[OMP_EXEC]"; };
  static constexpr int allocatorID() noexcept
  { return umpire::resource::Host; };
};

#endif

/*!
 * \brief execution_space traits specialization for SEQ_EXEC
 */
template < >
struct execution_space< SEQ_EXEC >
{
  using raja_exec                   = RAJA::loop_exec;
  using raja_reduce                 = RAJA::loop_reduce;
  using raja_atomic                 = RAJA::atomic::loop_atomic;
  static const int allocator_id = umpire::resource::Host;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[SEQ_EXEC]"; };
  static constexpr int allocatorID() noexcept
  { return umpire::resource::Host; };
};

/// @}
} /* namespace primal */
} /* namespace axom   */

#endif /* AXOM_PRIMAL_EXECUTIONSPACE_HPP_ */
