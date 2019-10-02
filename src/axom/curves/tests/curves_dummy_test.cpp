// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Curves includes
#include "axom/curves/curves.hpp"            /* for dummy() */

// Slic includes
#include "axom/slic/interface/slic.hpp"      /* for slic macros */
#include "axom/slic/core/UnitTestLogger.hpp" /* for UnitTestLogger */

// gtest includes
#include "gtest/gtest.h"  /* for gtest macros */

// C/C++ includes

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( curves, dummy )
{
  EXPECT_TRUE( axom::curves::dummy() );
}

//------------------------------------------------------------------------------
using axom::slic::UnitTestLogger;

int main( int argc, char* argv[] )
{
  int result = 0;
  ::testing::InitGoogleTest( &argc, argv );
  UnitTestLogger logger;
  result = RUN_ALL_TESTS();
  return result;
}
