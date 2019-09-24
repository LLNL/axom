// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

// slic
#include "axom/slic.hpp"

// spin includes
#include "axom/spin/execution_space.hpp"

#include "umpire/Umpire.hpp"  // for Umpire
#include "RAJA/RAJA.hpp"      // for RAJA

// for gtest macros
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>     // for strlen(),strcmp()
#include <type_traits> // for std::is_same

// namespace aliases
namespace spin = axom::spin;

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template< typename ExecSpace >
void check_valid( )
{
  SLIC_INFO( "checking execution space:" <<
              spin::execution_space< ExecSpace >::name() );
  EXPECT_TRUE( spin::execution_space< ExecSpace >::valid() );
  EXPECT_TRUE( strlen( spin::execution_space< ExecSpace >::name() ) > 0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace >
void check_invalid( )
{
  SLIC_INFO( "checking execution space:" <<
              spin::execution_space< ExecSpace >::name() );

  EXPECT_FALSE( spin::execution_space< ExecSpace >::valid() );

  EXPECT_EQ( spin::execution_space< ExecSpace >::allocatorID(),
             axom::INVALID_ALLOCATOR_ID );

  EXPECT_EQ( strcmp(spin::execution_space< ExecSpace >::name(),"[UNDEFINED]"),
             0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace,
           typename RajaExec,
           typename RajaReduce,
           typename RajaAtomic >
void check_execution_mappings( int expectedAllocatorID )
{
  SLIC_INFO( "checking execution space: " <<
             spin::execution_space< ExecSpace >::name());

  using exec   = typename spin::execution_space< ExecSpace >::raja_exec;
  using reduce = typename spin::execution_space< ExecSpace >::raja_reduce;
  using atomic = typename spin::execution_space< ExecSpace >::raja_atomic;

  bool valid_raja_exec   = std::is_same< exec, RajaExec >::value;
  bool valid_raja_reduce = std::is_same< reduce, RajaReduce >::value;
  bool valid_raja_atomic = std::is_same< atomic, RajaAtomic >::value;

  EXPECT_TRUE( valid_raja_exec );
  EXPECT_TRUE( valid_raja_reduce );
  EXPECT_TRUE( valid_raja_atomic );

  int allocatorID = spin::execution_space< ExecSpace >::allocatorID();
  EXPECT_EQ( expectedAllocatorID, allocatorID );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( spin_execution_space, check_valid )
{
  check_valid< spin::SEQ_EXEC >( );

#ifdef AXOM_USE_OPENMP
  check_valid< spin::OMP_EXEC >( );
#endif

#ifdef AXOM_USE_CUDA
  check_valid< spin::CUDA_EXEC< 256 > >( );
#endif
}

//------------------------------------------------------------------------------
TEST( spin_execution_space, check_invalid )
{
  struct InvalidSpace { };
  check_invalid< InvalidSpace >( );
}

//------------------------------------------------------------------------------
TEST( spin_execution_space, check_seq_exec )
{
  check_valid< spin::SEQ_EXEC >( );

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< spin::SEQ_EXEC,
                            RAJA::loop_exec,
                            RAJA::loop_reduce,
                            RAJA::loop_atomic >( allocator_id );

}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP
TEST( spin_execution_space, check_omp_exec )
{

  check_valid< spin::OMP_EXEC >( );

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Host );
  check_execution_mappings< spin::OMP_EXEC,
                            RAJA::omp_parallel_for_exec,
                            RAJA::omp_reduce,
                            RAJA::omp_atomic >( allocator_id );

}
#endif

//------------------------------------------------------------------------------
#ifdef AXOM_USE_CUDA
TEST( spin_execution_space, check_cuda_exec )
{
  constexpr int BLOCK_SIZE = 256;

  check_valid< spin::CUDA_EXEC< BLOCK_SIZE > >( );

  int allocator_id = axom::getResourceAllocatorID( umpire::resource::Unified );
  check_execution_mappings< spin::CUDA_EXEC< BLOCK_SIZE >,
                            RAJA::cuda_exec< BLOCK_SIZE >,
                            RAJA::cuda_reduce,
                            RAJA::cuda_atomic >( allocator_id );

}
#endif

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
