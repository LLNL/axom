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

#ifndef MINT_EXECUTION_POLICY_HPP_
#define MINT_EXECUTION_POLICY_HPP_

#include "axom/config.hpp"  // for AXOM_USE_RAJA

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace axom
{
namespace mint
{

/*!
 * \file
 *
 * \brief Defines the list of execution policies that Mint provides and an
 *  associated traits class.
 *
 *  \see interface.hpp
 */

//------------------------------------------------------------------------------
/// \name Execution Policy Traits
/// @{
//------------------------------------------------------------------------------

namespace policy
{

/// \name GPU Execution Policies
/// @{

/*!
 * \brief Parallel execution on the GPU.
 * \tparam BLOCK_SIZE the number of blocks to use
 *
 * \note This feature requires building Axom with RAJA.
 */
template < int BLOCK_SIZE >
struct parallel_gpu { };

/*!
 * \brief Asynchronous parallel execution on the GPU
 * \param BLOCK_SIZE the number of blocks to use.
 *
 * \note This feature requires building Axom with RAJA.
 *
 * \note When executing an asynchronous parallel execution policy with a
 *  loop traversal, the function returns immediately to the caller once
 *  the parallel execution is launched. It is up to the caller to handle
 *  appropriately the necessary synchronization.
 */
template < int BLOCK_SIZE >
struct parallel_gpu_async { };

/// @}

/// \name CPU Execution Policies
/// @{

/*!
 * \brief Parallel execution on the CPU.
 * \note This feature requires building Axom with RAJA.
 */
struct parallel_cpu { };

/*!
 * \brief Policy used for executing the traversal serially on the CPU.
 */
struct serial { };

/// @}

} /* namespace policy */

//------------------------------------------------------------------------------
/// @}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/// \name Execution Policy Traits
/// @{
//------------------------------------------------------------------------------

/*!
 * \brief Policy traits specialization
 */
template < typename ExecPolicy >
struct policy_traits
{
  using raja_exec_policy   = void;
  using raja_reduce_policy = void;
  using raja_sync_policy   = void;
  static constexpr bool valid() { return false; };
  static constexpr char* name() { return (char*)"[UNDEFINED]"; };
};

//------------------------------------------------------------------------------
// PARALLEL GPU POLICY TRAITS
//------------------------------------------------------------------------------
template < int BLOCK_SIZE >
struct policy_traits< policy::parallel_gpu< BLOCK_SIZE > >
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)
  using raja_exec_policy   = RAJA::cuda_exec< BLOCK_SIZE >;
  using raja_reduce_policy = RAJA::cuda_reduce;
  using raja_sync_policy   = RAJA::cuda_synchronize;

  /* *INDENT-OFF* */
  // TODO: use CudaCollapse policy when that is available
  using raja_2d_exec =
       RAJA::KernelPolicy<
          RAJA::statement::CudaKernelFixed< 256,
            RAJA::statement::For< 1, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
              RAJA::statement::For< 0, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
                RAJA::statement::Lambda< 0 >
              > // END i
            > // END j
          > // END CudaKernel
       >; // END kernel

  using raja_3d_exec =
      RAJA::KernelPolicy<
         RAJA::statement::CudaKernelFixed< 256,
           RAJA::statement::For< 2, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
             RAJA::statement::For< 1, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
               RAJA::statement::For< 0, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
                 RAJA::statement::Lambda< 0 >
             > // END i
            > // END j
          > // END k
        > // END CudaKernel
      >; // END kernel
  /* *INDENT-ON* */

#endif
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)"parallel_gpu"; };
  static constexpr int numBlocks( ) { return BLOCK_SIZE; };
};

//------------------------------------------------------------------------------
// PARALLEL GPU ASYNC POLICY TRAITS
//------------------------------------------------------------------------------
template < int BLOCK_SIZE >
struct policy_traits< policy::parallel_gpu_async< BLOCK_SIZE > >
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)
  using raja_exec_policy   = RAJA::cuda_exec_async< BLOCK_SIZE >;
  using raja_reduce_policy = RAJA::cuda_reduce;
  using raja_sync_policy   = RAJA::cuda_synchronize;

  /* *INDENT-OFF* */
  // TODO: use CudaCollapse policy when that is available
  using raja_2d_exec =
      RAJA::KernelPolicy<
         RAJA::statement::CudaKernelFixed< 256,
           RAJA::statement::For< 1, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
             RAJA::statement::For< 0, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
               RAJA::statement::Lambda< 0 >
             > // END i
           > // END j
         > // END CudaKernel
      >; // END kernel

  using raja_3d_exec =
      RAJA::KernelPolicy<
         RAJA::statement::CudaKernelFixed< 256,
           RAJA::statement::For< 2, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
             RAJA::statement::For< 1, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
               RAJA::statement::For< 0, RAJA::cuda_threadblock_exec<BLOCK_SIZE>,
                 RAJA::statement::Lambda< 0 >
             > // END i
            > // END j
          > // END k
         > // END CudaKernel
      >; // END kernel
  /* *INDENT-ON* */

#endif
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)"parallel_gpu_async"; };
  static constexpr int numBlocks( ) { return BLOCK_SIZE; };
};

//------------------------------------------------------------------------------
// PARALLEL CPU POLICY TRAITS
//------------------------------------------------------------------------------
template < >
struct policy_traits< policy::parallel_cpu >
{
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)
  using raja_exec_policy   = RAJA::omp_parallel_for_exec;
  using raja_reduce_policy = RAJA::omp_reduce;
  using raja_sync_policy   = RAJA::omp_synchronize;

  using raja_2d_exec =
          RAJA::KernelPolicy<
            RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                       RAJA::ArgList< 1,0 >,
                                       RAJA::statement::Lambda< 0 >
                                       > // END collapse
            >; // END kernel

  using raja_3d_exec =
          RAJA::KernelPolicy<
            RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                       RAJA::ArgList< 2,1,0 >,
                                       RAJA::statement::Lambda< 0 >
                                       > // END collapse
            >; // END kernel
#endif
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)"parallel_cpu"; };
};

//------------------------------------------------------------------------------
// SERIAL CPU POLICY TRAITS
//------------------------------------------------------------------------------
template < >
struct policy_traits< policy::serial >
{
#ifdef AXOM_USE_RAJA
  using raja_exec_policy   = RAJA::loop_exec;
  using raja_reduce_policy = RAJA::loop_reduce;
  using raja_sync_policy   = void;

  /* *INDENT-OFF* */
  using raja_2d_exec =
      RAJA::KernelPolicy<
         RAJA::statement::For< 1, RAJA::loop_exec,   // j
           RAJA::statement::For< 0, RAJA::loop_exec, // i
             RAJA::statement::Lambda< 0 >
           > // END i
         > // END j
      >; // END kernel

  using raja_3d_exec =
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

#endif
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)"serial"; };
};

//------------------------------------------------------------------------------
/// @}
//------------------------------------------------------------------------------

} /* namespace mint */
} /* namespace axom */


#endif /* MINT_EXECUTION_POLICY_HPP_ */
