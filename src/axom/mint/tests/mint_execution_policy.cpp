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

#include "axom/mint/execution/policy.hpp" // policy definitions and traits
#include "axom/slic/interface/slic.hpp"   // for SLIC macros

// gtest includes
#include "gtest/gtest.h"                  // for gtest macros

// C/C++ includes
#include <cstring>     // for strlen(),strcmp()
#include <type_traits> // for std::is_same

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace mint   = axom::mint;
namespace policy = mint::policy;

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

template < typename Pol,
           typename RajaExec,
           typename RajaReduce,
           typename RajaSync >
void check_policies( )
{
  SLIC_INFO("Check RAJA policies for: " << mint::policy_traits< Pol >::name() );

  using exec_policy   = typename mint::policy_traits< Pol >::raja_exec_policy;
  using reduce_policy = typename mint::policy_traits< Pol >::raja_reduce_policy;
  using sync_policy   = typename mint::policy_traits< Pol >::raja_sync_policy;

  bool valid_exec_policy   = std::is_same< exec_policy, RajaExec >::value;
  bool valid_reduce_policy = std::is_same< reduce_policy, RajaReduce >::value;
  bool valid_sync_policy   = std::is_same< sync_policy, RajaSync >::value;

  EXPECT_TRUE( valid_exec_policy );
  EXPECT_TRUE( valid_reduce_policy );
  EXPECT_TRUE( valid_sync_policy );
}

//------------------------------------------------------------------------------
template < typename MintPolicy >
void check_valid_mint_policy( )
{
  SLIC_INFO( "checking policy: " << mint::policy_traits< MintPolicy >::name() );
  EXPECT_TRUE( mint::policy_traits< MintPolicy >::valid() );
  EXPECT_TRUE( strlen( mint::policy_traits< MintPolicy >::name() ) > 0 );
}

//------------------------------------------------------------------------------
template < typename MintPolicy >
void check_invalid_mint_policy( )
{
  SLIC_INFO( "checking policy: " << mint::policy_traits< MintPolicy >::name() );
  EXPECT_FALSE( mint::policy_traits< MintPolicy >::valid() );
  EXPECT_EQ( strcmp(mint::policy_traits<MintPolicy>::name(),"[UNDEFINED]"), 0 );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_execution_policy, check_policy_traits )
{
  constexpr int BLOCK_SIZE = 1024;
  check_valid_mint_policy< policy::serial >( );
  check_valid_mint_policy< policy::parallel_cpu >( );
  check_valid_mint_policy< policy::parallel_gpu< BLOCK_SIZE > >( );
  check_valid_mint_policy< policy::parallel_gpu_async< BLOCK_SIZE > >( );

  struct InvalidPolicy { };
  check_invalid_mint_policy< InvalidPolicy >();
}

#ifdef AXOM_USE_RAJA

//------------------------------------------------------------------------------
TEST( mint_execution_policy, check_raja_mappings )
{
  check_policies< policy::serial, RAJA::loop_exec, RAJA::loop_reduce, void >( );

#if defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)
  check_policies< policy::parallel_cpu,
                  RAJA::omp_parallel_for_exec,
                  RAJA::omp_reduce,
                  RAJA::omp_synchronize >( );
#endif

#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)

  constexpr int BLOCK_SIZE = 1024;

  check_policies< policy::parallel_gpu< BLOCK_SIZE >,
                  RAJA::cuda_exec< BLOCK_SIZE >,
                  RAJA::cuda_reduce,
                  RAJA::cuda_synchronize >( );

  check_policies< policy::parallel_gpu_async< BLOCK_SIZE >,
                  RAJA::cuda_exec_async< BLOCK_SIZE >,
                  RAJA::cuda_reduce,
                  RAJA::cuda_synchronize >( );
#endif

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
