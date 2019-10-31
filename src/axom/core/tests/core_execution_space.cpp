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
           typename RajaExec,
           typename RajaReduce,
           typename RajaAtomic,
           typename RajaSync >
void check_execution_mappings( int expectedAllocatorID )
{
  std::cout << "checking execution space: " <<
               axom::execution_space< ExecSpace >::name() << std::endl;

  using exec   = typename axom::execution_space< ExecSpace >::loop_policy;
  using reduce = typename axom::execution_space< ExecSpace >::reduce_policy;
  using atomic = typename axom::execution_space< ExecSpace >::atomic_policy;
  using sync   = typename axom::execution_space< ExecSpace >::sync_policy;

  bool valid_raja_exec   = std::is_same< exec, RajaExec >::value;
  bool valid_raja_reduce = std::is_same< reduce, RajaReduce >::value;
  bool valid_raja_atomic = std::is_same< atomic, RajaAtomic >::value;
  bool valid_raja_sync   = std::is_same< sync, RajaSync >::value;

  EXPECT_TRUE( valid_raja_exec );
  EXPECT_TRUE( valid_raja_reduce );
  EXPECT_TRUE( valid_raja_atomic );
  EXPECT_TRUE( valid_raja_sync );

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

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< axom::SEQ_EXEC,
                            RAJA::loop_exec,
                            RAJA::loop_reduce,
                            RAJA::loop_atomic,
                            void >( allocator_id );

}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP
TEST( core_execution_space, check_omp_exec )
{

  check_valid< axom::OMP_EXEC >( );

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< axom::OMP_EXEC,
                            RAJA::omp_parallel_for_exec,
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

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Unified );
  check_execution_mappings< axom::CUDA_EXEC< BLOCK_SIZE >,
                            RAJA::cuda_exec< BLOCK_SIZE >,
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
