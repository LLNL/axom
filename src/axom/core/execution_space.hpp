// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_EXECUTIONSPACE_HPP_
#define AXOM_EXECUTIONSPACE_HPP_

#include "axom/config.hpp"
#include "axom/core/memory_management.hpp"

// RAJA includes
#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

// Umpire includes
#ifdef AXOM_USE_UMPIRE
#include "umpire/Umpire.hpp"
#endif

/*!
 * \file
 *
 * \brief Defines the list of available execution spaces for axom.
 *
 * The list of defined execution spaces are the following:
 *
 *  * <b>SEQ_EXEC<b> <br />
 *
 *    Indicates sequential execution on the CPU. Always defined.
 *
 *    When using this execution space, the data must reside on CPU/host memory.
 *
 *  * <b>OMP_EXEC<b> <br />
 *
 *    Indicates parallel execution with OpenMP on the CPU.
 *
 *    Defined when AXOM_USE_OPENMP and AXOM_USE_RAJA are defined. In addition,
 *    using this execution space requires linking to a RAJA that is configured
 *    with OpenMp, i.e., RAJA_ENABLE_OPENMP must be defined in the generated
 *    RAJA/config.hpp.
 *
 *    The default memory allocator when using this execution space is HOST.
 *
 *    When using this execution space, the data must reside on CPU/host memory.
 *
 *  * <b>CUDA_EXEC<BLOCKSIZE></b> <br />
 *
 *    Indicates parallel execution with CUDA on the GPU.
 *
 *    Defined when AXOM_USE_CUDA and AXOM_USE_RAJA are defined. In  addition,
 *    using this execution space requires linking to a RAJA that is configured
 *    with CUDA, i.e., RAJA_ENABLE_CUDA must be defined in the generated
 *    RAJA/config.hpp.
 *
 *    The CUDA_EXEC execution space strictly requires Umpire for memory
 *    management. The data must reside either on device memory which, is only
 *    accessible on the GPU, or on Unified memory which, is accessible on both
 *    the CPU or GPU. Note, when using Unified memory, the hardware transfers
 *    the data accordingly as needed. Consequently, using unified memory must
 *    be handled with extreme care, in order to avoid the latencies associated
 *    with data transfers between the CPU and GPU.
 *
 *    The default memory allocator when using this execution space is UNIFIED.
 *
 * The execution spaces bind the corresponding RAJA execution policies and
 * default memory space.
 */

namespace axom
{

/*!
 * \brief The execution_space is a traits class that binds the execution
 *  space to a corresponding RAJA execution policies and default memory
 *  allocator.
 *
 * \tparam ExecSpace the execution space
 *
 * \note This class is specialized for each execution space.
 *
 */
template < typename ExecSpace >
struct execution_space
{
  using loop_policy   = void;
  using loop2d_policy = void;
  using loop3d_policy = void;

  using reduce_policy = void;
  using atomic_policy = void;
  using sync_policy   = void;

  static constexpr bool valid() noexcept { return false; };
  static constexpr char* name() noexcept { return (char*)"[UNDEFINED]"; };
  static int allocatorID() noexcept { return axom::INVALID_ALLOCATOR_ID; };
};

/// \name Execution Spaces
/// @{

//-----------------------------------------------------------| SEQ_EXEC |-------

/// \name SEQ_EXEC
/// @{

/*!
 * \brief Indicates sequential execution on the CPU.
 */
struct SEQ_EXEC{ };

/*!
 * \brief execution_space traits specialization for SEQ_EXEC
 */
template < >
struct execution_space< SEQ_EXEC >
{

#ifdef AXOM_USE_RAJA
  using loop_policy   = RAJA::loop_exec;

  /* *INDENT-OFF* */
  using loop2d_policy =
      RAJA::KernelPolicy<
               RAJA::statement::For< 1, RAJA::loop_exec,   // j
                 RAJA::statement::For< 0, RAJA::loop_exec, // i
                   RAJA::statement::Lambda< 0 >
                 > // END i
               > // END j
            >; // END kernel

  using loop3d_policy =
      RAJA::KernelPolicy<
              RAJA::statement::For< 2, RAJA::loop_exec,       // k
                 RAJA::statement::For< 1, RAJA::loop_exec,    // j
                    RAJA::statement::For< 0, RAJA::loop_exec, // i
                       RAJA::statement::Lambda< 0 >
                    > // END i
                  > // END j
               > // END k
            >; // END kernel
  /* *INDENT-ON* */

  using reduce_policy = RAJA::loop_reduce;
  using atomic_policy = RAJA::loop_atomic;
#else
  using loop_policy   = void;
  using loop2d_policy = void;
  using loop3d_policy = void;

  using reduce_policy = void;
  using atomic_policy = void;
#endif

  using sync_policy = void;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[SEQ_EXEC]"; };
  static int allocatorID() noexcept
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getResourceAllocatorID( umpire::resource::Host );
#else
    return 0;
#endif
  };

};

/// @}

//-----------------------------------------------------------| OMP_EXEC |-------

/// \name OMP_EXEC
/// @{

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)

#if !defined(RAJA_ENABLE_OPENMP)
#error *** OMP_EXEC requires an OpenMP enabled RAJA ***
#endif

/*!
 * \brief Indicates parallel execution on the CPU using OpenMP.
 */
struct OMP_EXEC{ };

/*!
 * \brief execution_space traits specialization for OMP_EXEC
 */
template < >
struct execution_space< OMP_EXEC >
{
  using loop_policy = RAJA::omp_parallel_for_exec;

  /* *INDENT-OFF* */
  using loop2d_policy = RAJA::KernelPolicy<
      RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                 RAJA::ArgList< 1,0 >,
                                 RAJA::statement::Lambda< 0 > > >;

  using loop3d_policy = RAJA::KernelPolicy<
      RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                 RAJA::ArgList< 2,1,0 >,
                                 RAJA::statement::Lambda< 0 > > >;
  /* *INDENT-ON* */

  using reduce_policy = RAJA::omp_reduce;
  using atomic_policy = RAJA::omp_atomic;
  using sync_policy   = RAJA::omp_synchronize;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[OMP_EXEC]"; };

  static int allocatorID() noexcept
  {
#ifdef AXOM_USE_UMPIRE
    return axom::getResourceAllocatorID(umpire::resource::Host);
#else
    return 0;
#endif
  };

};

#endif

/// @}

//----------------------------------------------------------| CUDA_EXEC |-------

/// \name CUDA_EXEC
/// @{

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)

#if !defined(RAJA_ENABLE_CUDA)
#error *** CUDA_EXEC requires a CUDA enabled RAJA ***
#endif

#if !defined(UMPIRE_ENABLE_CUDA)
#error *** CUDA_EXEC requires a CUDA enabled UMPIRE ***
#endif

#if !defined(UMPIRE_ENABLE_UM)
#error *** CUDA_EXEC requires UMPIRE configured with UMPIRE_ENABLE_UM ***
#endif

/*!
 * \brief Indicates parallel execution on the GPU with CUDA.
 * \tparam BLOCK_SIZE the number of CUDA threads in a block.
 */
template < int BLOCK_SIZE >
struct CUDA_EXEC{ };

/*!
 * \brief execution_space traits specialization for CUDA_EXEC
 *
 * \tparam BLOCK_SIZE the number of CUDA threads to launch
 */
template < int BLOCK_SIZE >
struct execution_space< CUDA_EXEC< BLOCK_SIZE > >
{
  using loop_policy   = RAJA::cuda_exec< BLOCK_SIZE >;

  /* *INDENT-OFF* */
  using loop2d_policy =
      RAJA::KernelPolicy<
            RAJA::statement::CudaKernelFixed< 256,
              RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >;

  using loop3d_policy =
      RAJA::KernelPolicy<
        RAJA::statement::CudaKernelFixed< 256,
          RAJA::statement::For<2, RAJA::cuda_block_x_loop,
            RAJA::statement::For<1, RAJA::cuda_block_y_loop,
              RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                RAJA::statement::Lambda<0>
              >
            >
          >
        >
      >;
  /* *INDENT-ON* */

  using reduce_policy = RAJA::cuda_reduce;
  using atomic_policy = RAJA::cuda_atomic;
  using sync_policy   = RAJA::cuda_synchronize;

  static constexpr bool valid() noexcept { return true; };
  static constexpr char* name() noexcept { return (char*)"[CUDA_EXEC]"; };
  static int allocatorID() noexcept
  { return axom::getResourceAllocatorID(umpire::resource::Unified); };
};

#endif

/// @}

} /* namespace axom   */

#endif /* AXOM_SPIN_EXECUTIONSPACE_HPP_ */
