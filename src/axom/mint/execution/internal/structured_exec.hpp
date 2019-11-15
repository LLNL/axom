// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MINT_STRUCTURED_EXEC_HPP_
#define AXOM_MINT_STRUCTURED_EXEC_HPP_

#include "axom/core/execution/execution_space.hpp"

// RAJA includes
#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif


namespace axom
{
namespace mint
{
namespace internal
{

template < typename ExecSpace >
struct structured_exec
{
  using loop2d_policy = void;
  using loop3d_policy = void;
};

//--------------------------------------------------------| SEQ_EXEC |----------
template < >
struct structured_exec< SEQ_EXEC >
{
#ifdef AXOM_USE_RAJA
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
#else
  using loop2d_policy = void;
  using loop3d_policy = void;
#endif
};

//--------------------------------------------------------| OMP_EXEC |----------
#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
template < >
struct structured_exec< OMP_EXEC >
{
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
};
#endif

//--------------------------------------------------------| CUDA_EXEC |---------
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA)

template < int BLOCK_SIZE >
struct structured_exec< CUDA_EXEC< BLOCK_SIZE, SYNCHRONOUS > >
{
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

  /* *INDENT-ON*  */
};

template < int BLOCK_SIZE >
struct structured_exec< CUDA_EXEC< BLOCK_SIZE, ASYNC > >
{
  /* *INDENT-OFF* */

  using loop2d_policy =
        RAJA::KernelPolicy<
              RAJA::statement::CudaKernelFixedAsync< 256,
                RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                    RAJA::statement::Lambda<0>
                  >
                >
              >
            >;

    using loop3d_policy =
        RAJA::KernelPolicy<
          RAJA::statement::CudaKernelFixedAsync< 256,
            RAJA::statement::For<2, RAJA::cuda_block_x_loop,
              RAJA::statement::For<1, RAJA::cuda_block_y_loop,
                RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >
        >;

  /* *INDENT-ON*  */
};

#endif

} /* namespace internal */

} /* namespace mint */

} /* namespace axom */

#endif /* AXOM_MINT_STRUCTURED_EXEC_HPP_ */
