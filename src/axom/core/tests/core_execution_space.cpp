// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

// spin includes
#include "axom/core/execution_space.hpp"

#ifdef AXOM_USE_UMPIRE
#include "umpire/Umpire.hpp"  // for Umpire
#endif

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"      // for RAJA
#endif

// for gtest macros
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>     // for strlen(),strcmp()
#include <type_traits> // for std::is_same

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template< typename ExecSpace >
void check_valid( )
{
  std::cout << "checking execution space:" <<
               axom::execution_space< ExecSpace >::name() << std::endl;

  EXPECT_TRUE( axom::execution_space< ExecSpace >::valid() );
  EXPECT_TRUE( strlen( axom::execution_space< ExecSpace >::name() ) > 0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace >
void check_invalid( )
{
  std::cout << "checking execution space:" <<
              axom::execution_space< ExecSpace >::name() << std::endl;

  EXPECT_FALSE( axom::execution_space< ExecSpace >::valid() );

  EXPECT_EQ( axom::execution_space< ExecSpace >::allocatorID(),
             axom::INVALID_ALLOCATOR_ID );

  EXPECT_EQ( strcmp(axom::execution_space< ExecSpace >::name(),"[UNDEFINED]"),
             0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace,
           typename LoopPolicy,
           typename Loop2dPolicy,
           typename Loop3dPolicy,
           typename ReducePolicy,
           typename AtomicPolicy,
           typename SyncPolicy >
void check_execution_mappings( int expectedAllocatorID )
{
  std::cout << "checking execution space: " <<
               axom::execution_space< ExecSpace >::name() << std::endl;

  using loop_pol   = typename axom::execution_space< ExecSpace >::loop_policy;
  using loop2d_pol = typename axom::execution_space< ExecSpace >::loop2d_policy;
  using loop3d_pol = typename axom::execution_space< ExecSpace >::loop3d_policy;
  using reduce_pol = typename axom::execution_space< ExecSpace >::reduce_policy;
  using atomic_pol = typename axom::execution_space< ExecSpace >::atomic_policy;
  using sync_pol   = typename axom::execution_space< ExecSpace >::sync_policy;

  bool valid_loop_policy   = std::is_same< loop_pol, LoopPolicy >::value;
  bool valid_loop2d_policy = std::is_same< loop2d_pol, Loop2dPolicy >::value;
  bool valid_loop3d_policy = std::is_same< loop3d_pol, Loop3dPolicy >::value;
  bool valid_reduce_policy = std::is_same< reduce_pol, ReducePolicy >::value;
  bool valid_atomic_policy = std::is_same< atomic_pol, AtomicPolicy >::value;
  bool valid_sync_policy   = std::is_same< sync_pol, SyncPolicy >::value;

  EXPECT_TRUE( valid_loop_policy );
  EXPECT_TRUE( valid_loop2d_policy );
  EXPECT_TRUE( valid_loop3d_policy );
  EXPECT_TRUE( valid_reduce_policy );
  EXPECT_TRUE( valid_atomic_policy );
  EXPECT_TRUE( valid_sync_policy );

  int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  EXPECT_EQ( expectedAllocatorID, allocatorID );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( core_execution_space, check_valid )
{
  check_valid< axom::SEQ_EXEC >( );

#ifdef AXOM_USE_OPENMP
  check_valid< axom::OMP_EXEC >( );
#endif

#ifdef AXOM_USE_CUDA
  check_valid< axom::CUDA_EXEC< 256 > >( );
  check_valid< axom::CUDA_EXEC< 256, axom::ASYNC > >( );
#endif
}

//------------------------------------------------------------------------------
TEST( core_execution_space, check_invalid )
{
  struct InvalidSpace { };
  check_invalid< InvalidSpace >( );
}

//------------------------------------------------------------------------------
TEST( core_execution_space, check_seq_exec )
{
  check_valid< axom::SEQ_EXEC >( );

  /* *INDENT-OFF* */
  using raja_loop2d_policy =
      RAJA::KernelPolicy<
               RAJA::statement::For< 1, RAJA::loop_exec,   // j
                 RAJA::statement::For< 0, RAJA::loop_exec, // i
                   RAJA::statement::Lambda< 0 >
                 > // END i
               > // END j
            >; // END kernel

  using raja_loop3d_policy =
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

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< axom::SEQ_EXEC,
                            RAJA::loop_exec,
                            raja_loop2d_policy,
                            raja_loop3d_policy,
                            RAJA::loop_reduce,
                            RAJA::loop_atomic,
                            void >( allocator_id );

}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP
TEST( core_execution_space, check_omp_exec )
{

  check_valid< axom::OMP_EXEC >( );

  /* *INDENT-OFF* */
  using raja_loop2d_policy = RAJA::KernelPolicy<
      RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                 RAJA::ArgList< 1,0 >,
                                 RAJA::statement::Lambda< 0 > > >;

  using raja_loop3d_policy = RAJA::KernelPolicy<
      RAJA::statement::Collapse< RAJA::omp_parallel_collapse_exec,
                                 RAJA::ArgList< 2,1,0 >,
                                 RAJA::statement::Lambda< 0 > > >;
  /* *INDENT-ON* */

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< axom::OMP_EXEC,
                            RAJA::omp_parallel_for_exec,
                            raja_loop2d_policy,
                            raja_loop3d_policy,
                            RAJA::omp_reduce,
                            RAJA::omp_atomic,
                            RAJA::omp_synchronize >( allocator_id );

}
#endif

//------------------------------------------------------------------------------
#ifdef AXOM_USE_CUDA

TEST( core_execution_space, check_cuda_exec )
{
  constexpr int BLOCK_SIZE = 256;

  check_valid< axom::CUDA_EXEC< BLOCK_SIZE > >( );

  /* *INDENT-OFF* */
  using raja_loop2d_policy =
      RAJA::KernelPolicy<
            RAJA::statement::CudaKernelFixed< 256,
              RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >;

  using raja_loop3d_policy =
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

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Unified );
  check_execution_mappings< axom::CUDA_EXEC< BLOCK_SIZE >,
                            RAJA::cuda_exec< BLOCK_SIZE >,
                            raja_loop2d_policy,
                            raja_loop3d_policy,
                            RAJA::cuda_reduce,
                            RAJA::cuda_atomic,
                            RAJA::cuda_synchronize >( allocator_id );

}

//------------------------------------------------------------------------------
TEST( core_execution_space, check_cuda_exec_async )
{
  constexpr int BLOCK_SIZE = 256;

  check_valid< axom::CUDA_EXEC< BLOCK_SIZE, axom::ASYNC > >( );

  /* *INDENT-OFF* */
  using raja_loop2d_policy =
      RAJA::KernelPolicy<
            RAJA::statement::CudaKernelFixedAsync< 256,
              RAJA::statement::For<1, RAJA::cuda_block_x_loop,
                RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                  RAJA::statement::Lambda<0>
                >
              >
            >
          >;

  using raja_loop3d_policy =
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
  /* *INDENT-ON* */

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Unified );
  check_execution_mappings< axom::CUDA_EXEC< BLOCK_SIZE, axom::ASYNC >,
                            RAJA::cuda_exec_async< BLOCK_SIZE >,
                            raja_loop2d_policy,
                            raja_loop3d_policy,
                            RAJA::cuda_reduce,
                            RAJA::cuda_atomic,
                            RAJA::cuda_synchronize >( allocator_id );

}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
