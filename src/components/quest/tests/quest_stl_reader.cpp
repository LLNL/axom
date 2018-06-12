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

#include "quest/STLReader.hpp"
#include "mint/UnstructuredMesh.hpp"

// gtest includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cstdio>   // for std::remove()
#include <string>   // for std::string
#include <fstream>  // for std::ofstream

// namespace aliases
namespace mint  = axom::mint;
namespace quest = axom::quest;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

/*!
 * \brief Generates an STL file consisting of a single triangle on the XY plane
 * \param [in] file the name of the file to generate.
 * \pre file.empty() == false
 */
void generate_stl_file( const std::string& file )
{
  EXPECT_FALSE( file.empty() );

  std::ofstream ofs( file.c_str( ) );
  EXPECT_TRUE( ofs.is_open() );

  ofs << "solid triangle" << std::endl;
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex 0.0 0.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 1.0 0.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 0.0 1.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;
  ofs << "endsolid triangle" << std::endl;

  ofs.close( );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( quest_stl_reader_DeathTest, read_to_invalid_mesh )
{
  const char* IGNORE_OUTPUT  = ".*";
  const std::string filename = "triangle.stl";

  // STEP 0: generate a temporary STL file for testing
  generate_stl_file( filename );

  // STEP 1: constructs mesh object to read in the mesh to
  mint::UnstructuredMesh< mint::SINGLE_SHAPE > trimesh( 2, mint::TRIANGLE );
  mint::UnstructuredMesh< mint::SINGLE_SHAPE > hexmesh( 3, mint::HEX );

  // STEP 2: read in the STL mesh data
  quest::STLReader reader;
  reader.setFileName( filename );
  int status = reader.read( );
  EXPECT_EQ( status, 0 );

  // STEP 3: death tests

  // read the STL mesh data to a 2D mint::Mesh should fail
  EXPECT_DEATH_IF_SUPPORTED( reader.getMesh( &trimesh ), IGNORE_OUTPUT );

  // read the STL mesh data to a mint::Mesh that has a different cell type
  EXPECT_DEATH_IF_SUPPORTED( reader.getMesh( &hexmesh), IGNORE_OUTPUT );

  // STEP 4: remove STL file
  std::remove( filename.c_str() );
}

//------------------------------------------------------------------------------
TEST( quest_stl_reader, read_missing_file )
{
  const std::string INVALID_FILE = "foo.stl";
  quest::STLReader reader;
  reader.setFileName( INVALID_FILE );
  int status = reader.read( );
  EXPECT_TRUE( status != 0 );
}

//------------------------------------------------------------------------------
TEST( quest_stl_reader, read_stl )
{
  const double x_expected[] = { 0.0, 1.0, 0.0 };
  const double y_expected[] = { 0.0, 0.0, 1.0 };
  const double z_expected[] = { 0.0, 0.0, 0.0 };

  const std::string filename = "triangle.stl";

  // STEP 0: generate a temporary STL file for testing
  generate_stl_file( filename );

  // STEP 1: create an STL reader and read-in the mesh data
  quest::STLReader reader;
  reader.setFileName( filename );
  int status = reader.read();
  EXPECT_EQ( status, 0 );

  // STEP 2: reading the STL mesh data into a mint::Mesh
  mint::UnstructuredMesh< mint::SINGLE_SHAPE > mesh( 3, mint::TRIANGLE );
  reader.getMesh( &mesh );

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ( mesh.getNumberOfCells(), 1 );
  EXPECT_EQ( mesh.getNumberOfNodes(), 3 );

  const double* x = mesh.getCoordinateArray( mint::X_COORDINATE );
  const double* y = mesh.getCoordinateArray( mint::Y_COORDINATE );
  const double* z = mesh.getCoordinateArray( mint::Z_COORDINATE );
  EXPECT_TRUE( x != AXOM_NULLPTR );
  EXPECT_TRUE( y != AXOM_NULLPTR );
  EXPECT_TRUE( z != AXOM_NULLPTR );

  mint::IndexType numNodes = mesh.getNumberOfNodes();
  for ( mint::IndexType inode=0; inode < numNodes; ++inode )
  {
    EXPECT_DOUBLE_EQ( x[ inode ], x_expected[ inode ] );
    EXPECT_DOUBLE_EQ( y[ inode ], y_expected[ inode ] );
    EXPECT_DOUBLE_EQ( z[ inode ], z_expected[ inode ] );
  } // END for all nodes

  // STEP 4: remove temporary STL file
  std::remove( filename.c_str() );
}

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  // finalized when exiting main scope

  return RUN_ALL_TESTS();
}
