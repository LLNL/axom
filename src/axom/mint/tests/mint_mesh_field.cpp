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
#include "axom/Types.hpp"          // for axom Types
#include "mint/Field.hpp"          // for mint::Field
#include "mint/FieldTypes.hpp"     // for mint::field_traits
#include "mint/FieldVariable.hpp"  // for mint::FieldVariable
#include "mint/FEBasis.hpp"        // for FEBasisTypes

#include "slic/slic.hpp"           // for SLIC macros

// gtest includes
#include "gtest/gtest.h"

// namespace aliases
namespace mint   = axom::mint;
namespace common = axom::common;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
template < typename T >
void check_field_instantiation( )
{
  constexpr mint::IndexType NUM_TUPLES     = 4;
  constexpr mint::IndexType NUM_COMPONENTS = 1;


  mint::Field* f =
    new mint::FieldVariable< T >( "f", NUM_TUPLES, NUM_COMPONENTS );

  EXPECT_TRUE( f != AXOM_NULLPTR );
  EXPECT_EQ( f->getName(), "f" );
  EXPECT_EQ( f->getNumTuples(), NUM_TUPLES );
  EXPECT_EQ( f->getNumComponents(), NUM_COMPONENTS );
  EXPECT_EQ( f->getBasis(), MINT_UNDEFINED_BASIS );
  EXPECT_EQ( f->getType(), mint::field_traits< T >::type() );

  T* data = mint::Field::getDataPtr< T >( f );
  EXPECT_TRUE( data != AXOM_NULLPTR );

  delete f;
  f = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( mintMesh_field_DeathTest, invalid_operations )
{
  const char* IGNORE_OUTPUT = ".*";
  constexpr mint::IndexType NUM_TUPLES     = 4;
  constexpr mint::IndexType NUM_COMPONENTS = 1;

  mint::Field* f =
    new mint::FieldVariable< double >( "f", NUM_TUPLES, NUM_COMPONENTS );
  EXPECT_DEATH_IF_SUPPORTED( mint::Field::getDataPtr< int >( f ),
                             IGNORE_OUTPUT );

  delete f;
}

//------------------------------------------------------------------------------
TEST( mint_mesh_field, instantiate )
{
  check_field_instantiation< double >( );
  check_field_instantiation< float >( );
  check_field_instantiation< common::int32 >( );
  check_field_instantiation< common::int64 >( );
}

//------------------------------------------------------------------------------
TEST( mint_mesh_field, set_basis )
{
  constexpr mint::IndexType NUM_TUPLES     = 4;
  constexpr mint::IndexType NUM_COMPONENTS = 1;

  mint::Field* f =
    new mint::FieldVariable< double >( "f", NUM_TUPLES, NUM_COMPONENTS );
  EXPECT_TRUE( f != AXOM_NULLPTR );
  EXPECT_EQ( f->getBasis(), MINT_UNDEFINED_BASIS );

  f->setBasis( MINT_LAGRANGE_BASIS );
  EXPECT_EQ( f->getBasis(), MINT_LAGRANGE_BASIS );

  delete f;
  f = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
TEST( mint_mesh_field, get_dataptr )
{
  constexpr mint::IndexType NUM_TUPLES     = 4;
  constexpr mint::IndexType NUM_COMPONENTS = 1;

  double f[]  = { 1.0, 2.0, 3.0, 4.0 };
  mint::Field* field =
    new mint::FieldVariable< double >( "f", f, NUM_TUPLES, NUM_COMPONENTS );

  const double* data = mint::Field::getDataPtr< double >( field );
  EXPECT_TRUE( data != AXOM_NULLPTR );
  EXPECT_EQ( data, f );

  for ( int i=0 ; i < NUM_TUPLES ; ++i )
  {
    EXPECT_DOUBLE_EQ( data[ i ], f[ i ] );
  }

  delete field;
  field = AXOM_NULLPTR;

  EXPECT_TRUE( f != AXOM_NULLPTR );
  for ( int i=0 ; i < NUM_TUPLES ; ++i )
  {
    EXPECT_DOUBLE_EQ( f[ i ], static_cast< double >( i+1 ) );
  }

}

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
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
