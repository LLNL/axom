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

#include "axom/mint/execution/xargs.hpp" // execution xargs types/traits
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

// gtest includes
#include "gtest/gtest.h"                 // for gtest macros

// C/C++ includes
#include <cstring>

namespace mint  = axom::mint;
namespace xargs = mint::xargs;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template < typename ArgType >
void check_valid( )
{
  SLIC_INFO( "checking [" << mint::xargs_traits< ArgType >::name() << "]" );
  EXPECT_TRUE( mint::xargs_traits< ArgType >::valid() );
  EXPECT_TRUE( strlen( mint::xargs_traits< ArgType >::name() ) > 0 );
}

//------------------------------------------------------------------------------
template < typename ArgType >
void check_invalid( )
{
  SLIC_INFO( "checking invalid policy!");
  EXPECT_FALSE( mint::xargs_traits< ArgType >::valid() );
  EXPECT_TRUE( strlen( mint::xargs_traits< ArgType >::name() ) == 0 );
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_execution_xargs, check_types )
{
  check_valid< xargs::index >( );
  check_valid< xargs::ij >( );
  check_valid< xargs::ijk >( );
  check_valid< xargs::x >( );
  check_valid< xargs::xy >( );
  check_valid< xargs::xyz >( );
  check_valid< xargs::nodeids >( );
  check_valid< xargs::coords >( );
  check_valid< xargs::faceids >( );
  check_valid< xargs::cellids >( );

  struct invalid_type { };
  check_invalid< invalid_type >( );
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
