/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

// Axom includes
#include "axom/config.hpp"                    // for compile-time definitions

// Mint includes
#include "axom/mint/config.hpp"               // mint compile-time definitions
#include "axom/mint/execution/policy.hpp"     // mint execution policies/traits
#include "axom/mint/execution/interface.hpp"  // mint::for_all()

// Slic includes
#include "axom/slic.hpp" // for SLIC macros

// gtest includes
#include "gtest/gtest.h" // for gtest

// namespace aliases
namespace mint   = axom::mint;
namespace policy = mint::policy;

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template < typename ExecPolicy >
void check_for_all( bool async=false )
{
  EXPECT_TRUE( mint::policy_traits< ExecPolicy >::valid() );
  SLIC_INFO( "check_for_all with [" <<
             mint::policy_traits< ExecPolicy >::name() << "]" );

  // STEP 0: set some constants
  constexpr int VALUE_1 = -42;
  constexpr int VALUE_2 =  42;
  constexpr int N       = 256;

  // STEP 0: allocate buffer
  int* a = axom::alloc< int >( N );

  // STEP 1: initialize to VALUE_1
  mint::for_all< ExecPolicy >( N,
    AXOM_LAMBDA(mint::IndexType idx)
    {
      a[ idx ] = VALUE_1;
    }
  );

  if ( async )
  {
    mint::synchronize< ExecPolicy >( );
  }

  // STEP 2: check array
  for ( int i=0 ; i < N ; ++i )
  {
    EXPECT_EQ( a[ i ], VALUE_1 );
  }

  // STEP 3: set all values to VALUE_2 with mint::for_all
  mint::for_all< ExecPolicy >( 0, N,
    AXOM_LAMBDA(mint::IndexType idx)
    {
      a[ idx ] = VALUE_2;
    }
  );

  if ( async )
  {
    mint::synchronize< ExecPolicy >( );
  }

  // STEP 4: check array
  for ( int i=0 ; i < N ; ++i )
  {
    EXPECT_EQ( a[ i ], VALUE_2 );
  }

  // STEP 2:
  axom::free( a );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_execution_forall_traversals, check_generic_loop )
{
  check_for_all< policy::serial >( );

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP) && \
  defined(RAJA_ENABLE_OPENMP)

  check_for_all< policy::parallel_cpu >( );
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && \
  defined(RAJA_ENABLE_CUDA)

  constexpr bool async = true;
  check_for_all< policy::parallel_gpu< 512 > >( );
  check_for_all< policy::parallel_gpu_async< 512 > >( async );
#endif

}

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
