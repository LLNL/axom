// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for compile time definitions

// slic
#include "axom/slic.hpp"

// primal includes
#include "axom/primal/spatial_acceleration/ExecutionSpace.hpp"

#include "umpire/Umpire.hpp"  // for Umpire
#include "RAJA/RAJA.hpp"      // for RAJA

// for gtest macros
#include "gtest/gtest.h"

// C/C++ includes
#include <cstring>     // for strlen(),strcmp()
#include <type_traits> // for std::is_same

// namespace aliases
namespace primal = axom::primal;

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template< typename ExecSpace >
void check_valid( )
{
  SLIC_INFO( "checking execution space:" <<
              primal::execution_space< ExecSpace >::name() );
  EXPECT_TRUE( primal::execution_space< ExecSpace >::valid() );
  EXPECT_TRUE( strlen( primal::execution_space< ExecSpace >::name() ) > 0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace >
void check_invalid( )
{
  SLIC_INFO( "checking execution space:" <<
              primal::execution_space< ExecSpace >::name() );

  EXPECT_FALSE( primal::execution_space< ExecSpace >::valid() );

  EXPECT_EQ( primal::execution_space< ExecSpace >::allocatorID(),
             axom::INVALID_ALLOCATOR_ID );

  EXPECT_EQ( strcmp(primal::execution_space< ExecSpace >::name(),"[UNDEFINED]"),
             0 );
}

//------------------------------------------------------------------------------
template < typename ExecSpace,
           typename RajaExec,
           typename RajaReduce,
           typename RajaAtomic,
           umpire::resource::MemoryResourceType UmpireAllocator >
void check_execution_mappings( )
{
  SLIC_INFO( "checking execution space: " <<
             primal::execution_space< ExecSpace >::name());

  using exec   = typename primal::execution_space< ExecSpace >::raja_exec;
  using reduce = typename primal::execution_space< ExecSpace >::raja_reduce;
  using atomic = typename primal::execution_space< ExecSpace >::raja_atomic;

  bool valid_raja_exec   = std::is_same< exec, RajaExec >::value;
  bool valid_raja_reduce = std::is_same< reduce, RajaReduce >::value;
  bool valid_raja_atomic = std::is_same< atomic, RajaAtomic >::value;

  EXPECT_TRUE( valid_raja_exec );
  EXPECT_TRUE( valid_raja_reduce );
  EXPECT_TRUE( valid_raja_atomic );

  int allocatorID = primal::execution_space< ExecSpace >::allocatorID();
  EXPECT_EQ( allocatorID, UmpireAllocator );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST( primal_execution_space, check_valid )
{
  check_valid< primal::SEQ_EXEC >( );

#ifdef AXOM_USE_OPENMP
  check_valid< primal::OMP_EXEC >( );
#endif

#ifdef AXOM_USE_CUDA
  check_valid< primal::CUDA_EXEC< 256 > >( );
#endif
}

//------------------------------------------------------------------------------
TEST( primal_execution_space, check_invalid )
{
  struct InvalidSpace { };
  check_invalid< InvalidSpace >( );
}

//------------------------------------------------------------------------------
TEST( primal_execution_space, check_seq_exec )
{
  check_valid< primal::SEQ_EXEC >( );

  check_execution_mappings< primal::SEQ_EXEC,
                            RAJA::loop_exec,
                            RAJA::loop_reduce,
                            RAJA::atomic::loop_atomic,
                            umpire::resource::Host >( );

}

//------------------------------------------------------------------------------
#ifdef AXOM_USE_OPENMP
TEST( primal_execution_space, check_omp_exec )
{

  check_valid< primal::OMP_EXEC >( );

  check_execution_mappings< primal::OMP_EXEC,
                            RAJA::omp_parallel_for_exec,
                            RAJA::omp_reduce,
                            RAJA::atomic::omp_atomic,
                            umpire::resource::Host >( );

}
#endif

//------------------------------------------------------------------------------
#ifdef AXOM_USE_CUDA
TEST( primal_execution_space, check_cuda_exec )
{
  constexpr int BLOCK_SIZE = 256;

  check_valid< primal::CUDA_EXEC< BLOCK_SIZE > >( );

  check_execution_mappings< primal::CUDA_EXEC< BLOCK_SIZE >,
                            RAJA::cuda_exec< BLOCK_SIZE >,
                            RAJA::cuda_reduce,
                            RAJA::atomic::cuda_atomic,
                            umpire::resource::Unified >( );

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
