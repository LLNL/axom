/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Dec 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

// ATK Toolkit includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "common/FileUtilities.hpp"
#include "common/Timer.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/BVHTree.hpp"
#include "quest/Point.hpp"
#include "quest/STLReader.hpp"
#include "quest/Triangle.hpp"
#include "quest/SignedDistance.hpp"

#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/UnstructuredMesh.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

using namespace asctoolkit;

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

struct {
  std::string fileName;
  int maxLevels;
  int maxObjects;
  int nx;
  int ny;
  int nz;
} Arguments;

//------------------------------------------------------------------------------
void showHelp()
{
  SLIC_INFO("Usage: ./quest_driver --file <myfile.stl> [options]" );
  SLIC_INFO("--file <file> specifies the STL file of the input surface mesh.");
  SLIC_INFO("--maxLevels <N> max levels for BVH decomposition.");
  SLIC_INFO("--maxObjects <N> max objects in a BVH bin (for refinement).");
  SLIC_INFO("--nx <N> number of cells in x-direction for sample grid.");
  SLIC_INFO("--ny <N> number of cells in y-direction for sample grid.");
  SLIC_INFO("--nz <N> number of cells in z-direction for sample grid.");
  SLIC_INFO("--help prints this help information.");
}

//------------------------------------------------------------------------------
void parse_args( int argc, char** argv )
{
  // Defaults
  Arguments.fileName   = "";
  Arguments.maxLevels  = 15;
  Arguments.maxObjects = 5;
  Arguments.nx         = 32;
  Arguments.ny         = 32;
  Arguments.nz         = 32;

  for ( int i=1; i < argc; ++i ) {

    if ( strcmp( argv[i],"--file")==0 ) {

      Arguments.fileName = std::string( argv[ ++i ] );

    } else if ( strcmp(argv[i],"--maxLevels")==0 ) {

      Arguments.maxLevels = std::atoi( argv[++i] );

    } else if ( strcmp(argv[i], "--maxObjects" )==0 ) {

      Arguments.maxObjects = std::atoi( argv[++i] );

    } else if ( strcmp(argv[i],"--nx")==0 ) {

      Arguments.nx = std::atoi( argv[++i] );

    } else if ( strcmp(argv[i],"--ny")==0 ) {

      Arguments.ny = std::atoi( argv[++i] );

    } else if ( strcmp(argv[i],"--nz")==0 ) {

      Arguments.nz = std::atoi( argv[++i] );

    } else if ( strcmp(argv[i],"--help")==0 ) {
      showHelp();
      exit(0);
    }

  } // END for

  if ( Arguments.fileName == "" ) {
    SLIC_ERROR( "NO STL input file provided. Provide one with --file" );
  }

}

//------------------------------------------------------------------------------
void write_point( const quest::Point< double, 3 >& pt,
                  const std::string& fileName )
{
  std::ofstream ofs;
  ofs.open( fileName.c_str() );
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Point " << fileName << "\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";
  ofs << "POINTS 1 double\n";
  ofs << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
  ofs << "CELLS 1 2\n";
  ofs << "1 0\n";
  ofs << "CELL_TYPES 1\n";
  ofs << "1\n";
  ofs.close();
}

//------------------------------------------------------------------------------
void write_vtk( mint::Mesh* mesh, const std::string& fileName )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  std::ofstream ofs;
  ofs.open( fileName.c_str() );

  // STEP 0: Write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " Unstructured Mesh (";
  ofs << mesh->getBlockId() << ", " << mesh->getPartitionId() << ")\n";

  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: Write mesh nodes
  const int num_nodes = mesh->getMeshNumberOfNodes();
  ofs << "POINTS " << num_nodes << " double\n";
  for ( int node=0; node < num_nodes; ++node ) {

    ofs << mesh->getMeshNodeCoordinate( node, 0 ) << " ";
    ofs << mesh->getMeshNodeCoordinate( node, 1 ) << " ";
    if ( mesh->getDimension()==3 ) {

        ofs << mesh->getMeshNodeCoordinate( node, 2 );

    } else {

        ofs << "0.0";
    }

    ofs << std::endl;

  } // END for all nodes

  // STEP 2: Write mesh cell connectivity
  // TODO: note currently this does not work with mixed cell types
  const int ncells   = mesh->getMeshNumberOfCells();
  const int maxnodes = mesh->getMeshNumberOfCellNodes( 0 );
  ofs << "CELLS " << ncells << " " << ncells*(maxnodes+1) << std::endl;

  std::vector< int > cell;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {

    const int nnodes = mesh->getMeshNumberOfCellNodes( cellIdx );
    cell.resize( nnodes );
    mesh->getMeshCell( cellIdx, &cell[0] );

    ofs << nnodes << " ";
    for ( int i=0; i < nnodes; ++i ) {
      ofs << cell[ i ] << " ";
    } // END for all nodes
    ofs << std::endl;

  } // END for all cells

  // STEP 3: Write cell types
  ofs << "CELL_TYPES " << ncells << std::endl;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {
    int ctype    = mesh->getMeshCellType( cellIdx );
    int vtk_type = mint::cell::vtk_types[ ctype ];
    ofs << vtk_type << std::endl;
  } // END for all cells

  // STEP 4: Write Cell Data
  ofs << "CELL_DATA " << ncells << std::endl;
  mint::FieldData* CD = mesh->getCellFieldData();
  for ( int f=0; f < CD->getNumberOfFields(); ++f ) {

      mint::Field* field = CD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == mint::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i ) {
             ofs << dataPtr[ i ] << std::endl;
          }
      }

  }

  // STEP 5: Write Point Data
  const int nnodes = mesh->getMeshNumberOfNodes();
  ofs << "POINT_DATA " << nnodes << std::endl;
  mint::FieldData* PD = mesh->getNodeFieldData();
  for ( int f=0; f < PD->getNumberOfFields(); ++f ) {

      mint::Field* field = PD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == mint::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      }
  }

  ofs.close();
}

//------------------------------------------------------------------------------
void write_triangles( mint::Mesh* mesh, const int* cells, int ncells,
                      const std::string& fileName )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  TriangleMesh* subset = new TriangleMesh(3);
  SLIC_ASSERT( subset->getMeshType() == mesh->getMeshType() );

  int cellIds[3];
  quest::Point< double, 3 > n1;
  quest::Point< double, 3 > n2;
  quest::Point< double, 3 > n3;

  int icount = 0;
  int new_cell[3];

  for ( int i=0; i < ncells; ++i ) {

     const int cellIdx = cells[ i ];
     mesh->getMeshCell( cellIdx, cellIds );

     mesh->getMeshNode( cellIds[0], n1.data() );
     mesh->getMeshNode( cellIds[1], n2.data() );
     mesh->getMeshNode( cellIds[2], n3.data() );

     subset->insertNode( n1[0], n1[1], n1[2] );
     new_cell[0] = icount; ++icount;
     subset->insertNode( n2[0], n2[1], n2[2] );
     new_cell[1] = icount; ++icount;
     subset->insertNode( n3[0], n3[1], n3[2] );
     new_cell[2] = icount; ++icount;

     subset->insertCell( new_cell, MINT_TRIANGLE,3 );
  }

  write_vtk( subset, fileName );

  delete subset;
}

//------------------------------------------------------------------------------
quest::BoundingBox< double,3 > compute_bounds( mint::Mesh* mesh)
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   using namespace quest;

   quest::BoundingBox< double,3 > meshBB;
   quest::Point< double,3 > pt;

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
   {
       mesh->getMeshNode( i, pt.data() );
       meshBB.addPoint( pt );
   } // END for all nodes

   SLIC_ASSERT( meshBB.isValid() );

   return meshBB;
}


//------------------------------------------------------------------------------
void distance_field( mint::Mesh* surface_mesh, mint::UniformMesh* umesh )
{
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
  SLIC_ASSERT( umesh != ATK_NULLPTR );

  SLIC_INFO("Max BVH levels: " << Arguments.maxLevels );
  SLIC_INFO("Max Object Threshold: " << Arguments.maxObjects );
  utilities::Timer timer1;
  timer1.start();

  quest::SignedDistance< 3 > signedDistance(
      surface_mesh, Arguments.maxObjects, Arguments.maxLevels);

  timer1.stop();
  SLIC_INFO("Constructed BVH in " << timer1.elapsed() << "s" );

#ifdef ATK_DEBUG
  // write the bucket tree to a file
  const quest::BVHTree< int, 3>* btree = signedDistance.getBVHTree();
  SLIC_ASSERT( btree != ATK_NULLPTR );

  btree->writeVtkFile( "bucket-tree.vtk" );

  // mark bucket IDs on surface mesh
  const int ncells = surface_mesh->getMeshNumberOfCells();
  mint::FieldData* CD = surface_mesh->getCellFieldData();
  CD->addField( new mint::FieldVariable<int>( "BucketID", ncells ) );
  int* bidx = CD->getField( "BucketID" )->getIntPtr();
  SLIC_ASSERT( bidx != ATK_NULLPTR );

  const int numObjects = btree->getNumberOfObjects();
  for ( int i=0; i < numObjects; ++i ) {

     const int idx = btree->getObjectBucketIndex( i );
     bidx[ i ] = idx;

  } // END for all objects

  write_vtk( surface_mesh, "partitioned_surface_mesh.vtk" );
#endif

  const int nnodes = umesh->getNumberOfNodes();
  mint::FieldData* PD = umesh->getNodeFieldData();
  SLIC_ASSERT( PD != ATK_NULLPTR );

  PD->addField( new mint::FieldVariable< double >("phi",nnodes) );
  PD->addField( new mint::FieldVariable< int >("nbuckets",nnodes) );
  PD->addField( new mint::FieldVariable< int >("ntriangles", nnodes) );

  double* phi     = PD->getField( "phi" )->getDoublePtr();
  int*  nbuckets  = PD->getField( "nbuckets" )->getIntPtr();
  int* ntriangles = PD->getField( "ntriangles" )->getIntPtr();

  SLIC_ASSERT( phi != ATK_NULLPTR );
  SLIC_ASSERT( nbuckets != ATK_NULLPTR );
  SLIC_ASSERT( ntriangles != ATK_NULLPTR );

  utilities::Timer timer2;
  timer2.start();

  for ( int inode=0; inode < nnodes; ++inode ) {

      quest::Point< double,3 > pt;
      umesh->getMeshNode( inode, pt.data() );

      std::vector< int > buckets;
      std::vector< int > triangles;
      std::vector< int > my_triangles;
      triangles.clear();
      buckets.clear();

      quest::Point< double,3 > closest_pt;
      phi[ inode ] = signedDistance.computeDistance( pt,
                                                     buckets,
                                                     triangles,
                                                     my_triangles,
                                                     closest_pt );

      nbuckets[ inode ]   = static_cast< int >( buckets.size() );
      ntriangles[ inode ] = static_cast< int >( triangles.size() );

#ifdef ATK_DEBUG
      std::ostringstream oss;
      oss << "BINS_" << inode << ".vtk";
      signedDistance.getBVHTree()->writeVtkFile(
            oss.str(), &buckets[0], buckets.size() );

      oss.str("");
      oss.clear();
      oss << "POINT_" << inode << ".vtk";
      write_point( pt, oss.str() );

      oss.str("");
      oss.clear();
      oss << "CLOSEST_POINT_" << inode << ".vtk";
      write_point( closest_pt, oss.str() );

      oss.str("");
      oss.clear( );
      oss << "TRIANGLES_" << inode << ".vtk";
      write_triangles( surface_mesh, &triangles[0],
                       ntriangles[inode], oss.str() );

      oss.str("");
      oss.clear();
      oss << "MY_TRIANGLES_" << inode << ".vtk";
      write_triangles( surface_mesh, &my_triangles[0], my_triangles.size(),
                       oss.str() );
#endif

  } // END for all nodes

  timer2.stop();
  SLIC_INFO("BVH Query Time: " << timer2.elapsed() << "s" );

}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingMsgLevel( asctoolkit::slic::message::Debug );

  // Create a more verbose message for this application (only level and message)
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  slic::GenericOutputStream* defaultStream =
          new slic::GenericOutputStream(&std::cout);
  slic::GenericOutputStream* compactStream =
          new slic::GenericOutputStream(&std::cout, slicFormatStr);
  slic::addStreamToMsgLevel(defaultStream, asctoolkit::slic::message::Fatal) ;
  slic::addStreamToMsgLevel(defaultStream, asctoolkit::slic::message::Error);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Warning);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Info);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Debug);

  // STEP 1: get file from user or use default
  parse_args( argc, argv );

  // STEP 2: read file
  SLIC_INFO( "Reading file: " << Arguments.fileName << "...");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( Arguments.fileName );
  reader->read();
  SLIC_INFO("done");

  // STEP 3: get surface mesh
  mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );
  SLIC_INFO("Mesh has "
          << surface_mesh->getMeshNumberOfNodes() << " nodes and "
          << surface_mesh->getMeshNumberOfCells() << " cells.");

  // STEP 4: Delete the reader
  delete reader;
  reader = ATK_NULLPTR;

  // STEP 5: write vtk file
  write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 6: compute bounds
  quest::BoundingBox< double,3 > meshBB = compute_bounds( surface_mesh);
  SLIC_INFO("Mesh bounding box: " << meshBB );


  quest::BoundingBox< double,3 > queryBounds = meshBB;
  queryBounds.expand( 10 );

  double h[3];
  const quest::Vector< double,3 >& bbDiff = queryBounds.range();
  h[0] = bbDiff[0] / Arguments.nx;
  h[1] = bbDiff[1] / Arguments.ny;
  h[2] = bbDiff[2] / Arguments.nz;
  SLIC_INFO("grid dimensions:" << Arguments.nx << ", "
                               << Arguments.ny << ", "
                               << Arguments.nz );
  SLIC_INFO("grid cell size:(" << h[0] << "," << h[1] << "," << h[2] << ")\n" );

  int node_ext[6];
  node_ext[0] = 0;
  node_ext[1] = Arguments.nx;
  node_ext[2] = 0;
  node_ext[3] = Arguments.ny;
  node_ext[4] = 0;
  node_ext[5] = Arguments.nz;

  // STEP 8: Construct uniform mesh
  mint::UniformMesh* umesh =
          new mint::UniformMesh(3, queryBounds.getMin().data(), h, node_ext);

  // STEP 9: Compute the distance field on the uniform mesh
  SLIC_INFO( "computing distance field..." );
  utilities::Timer timer;
  timer.start();
  distance_field( surface_mesh, umesh );
  timer.stop();
  SLIC_INFO( "Time to compute distance field: " << timer.elapsed() );

  // STEP 10: write the uniform mesh
  write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 11: clean up
  delete surface_mesh;
  surface_mesh = ATK_NULLPTR;

  delete umesh;
  umesh = ATK_NULLPTR;

  // STEP 12: Finalize SLIC environment
  slic::finalize();
  return 0;
}
