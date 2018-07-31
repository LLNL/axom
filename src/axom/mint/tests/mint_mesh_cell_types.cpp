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
#include "axom/slic/interface/slic.hpp"      // for slic macros
#include "axom/mint/mesh/CellTypes.hpp" // for CellTypes

// gtest includes
#include "gtest/gtest.h"      // for gtest macros

// C/C++ includes
#include <cstring>  // for strlen()

// namespace aliases
namespace mint = axom::mint;


//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_mesh_cell_types, check_cell_types )
{
  const int num_nodes[ mint::NUM_CELL_TYPES ] = {
    1,    // VERTEX
    2,    // SEGMENT
    3,    // TRIANGLE
    4,    // QUAD
    4,    // TET
    8,    // HEX
    6,    // PRISM
    5,    // PYRAMID
    9,    // QUAD9
    27,   // HEX27
  };

  for ( int i=0 ; i < mint::NUM_CELL_TYPES ; ++i )
  {
    const int cell_id = cellTypeToInt( mint::cell_info[ i ].cell_type );
    SLIC_INFO( "[" << mint::cell_info[ i ].name << "] =>" <<
               "type=" << cell_id << " " <<
               "blueprint_name=(" << mint::cell_info[ i ].blueprint_name << ") " <<
               "vtk_type=" << mint::cell_info[ i ].vtk_type << " " <<
               "num_nodes=" << mint::cell_info[ i ].num_nodes
               );

    EXPECT_EQ( cell_id, i );
    EXPECT_TRUE( strlen( mint::cell_info[ i ].name ) > 0 );
    EXPECT_TRUE( strlen( mint::cell_info[ i ].blueprint_name ) > 0 );
    EXPECT_TRUE( mint::cell_info[ i ].vtk_type > 0 );
    EXPECT_EQ( mint::cell_info[ i ].num_nodes, num_nodes[ i ] );

  } // END for all types
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
