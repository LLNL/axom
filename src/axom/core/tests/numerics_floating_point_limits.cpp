// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom includes
#include "axom/core/Macros.hpp"
#include "axom/core/numerics/floating_point_limits.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <type_traits> // for std::is_floating_point

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template < typename T >
void check_type_limits( const std::string& typeName )
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

  SCOPED_TRACE( "Testing [" + typeName + "]" );

  EXPECT_DOUBLE_EQ( axom::numerics::floating_point_limits< T >::lowest(),
                    std::numeric_limits< T >::lowest() );

  EXPECT_DOUBLE_EQ( axom::numerics::floating_point_limits< T >::min(),
                    std::numeric_limits< T >::min() );

  EXPECT_DOUBLE_EQ( axom::numerics::floating_point_limits< T >::max(),
                    std::numeric_limits< T >::max() );
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( numerics_floating_point_limits, consistency_with_standard_numeric_limits )
{
  check_type_limits< float >( "float" );
  check_type_limits< double >( "double" );
  check_type_limits< long double >( "long double" );
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
