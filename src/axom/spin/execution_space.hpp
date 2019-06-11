// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_SPIN_EXECUTIONSPACE_HPP_
#define AXOM_SPIN_EXECUTIONSPACE_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

#if !defined(AXOM_USE_RAJA) || !defined(AXOM_USE_UMPIRE)
#error *** The execution_space traits class requires RAJA and Umpire ***
#endif


// RAJA includes
#include "RAJA/RAJA.hpp"

// Umpire includes
#include "umpire/Umpire.hpp"

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
namespace spin
{

/// \name Execution Spaces
/// @{

/*!
 * \brief Indicates sequential execution on the CPU.
 */
struct SEQ_EXEC{ };

/*!
 * \brief Indicates parallel execution on the CPU using OpenMP.
 */
#ifdef AXOM_USE_OPENMP
struct OMP_EXEC{ };
#endif


/*!
 * \brief Indicates parallel execution on the GPU with CUDA.
 * \tparam BLOCK_SIZE the number of CUDA threads in a block.
 */
#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
template < int BLOCK_SIZE >
struct CUDA_EXEC{ };
#endif

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
} /* namespace spin */
} /* namespace axom   */

#endif /* AXOM_SPIN_EXECUTIONSPACE_HPP_ */
