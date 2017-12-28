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

// axom includes
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "axom_utils/FileUtilities.hpp"
#include "axom_utils/Timer.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/BVHTree.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"

#include "quest/STLReader.hpp"
#include "quest/SignedDistance.hpp"

#include "mint/config.hpp"
#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

using namespace axom;

using axom::primal::Point;
using axom::primal::Triangle;
using axom::primal::BoundingBox;
using axom::primal::BVHTree;
using axom::primal::Vector;

using quest::SignedDistance;

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

static struct
{
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
  SLIC_INFO("Usage: ./quest_distance_driver_ex --file <myfile.stl> [options]" );
  SLIC_INFO("--file <file> specifies the STL file of the input surface mesh.");
  SLIC_INFO("--maxLevels <N> max levels for BVH decomposition.");
  SLIC_INFO("--maxObjects <N> max objects in a BVH bin (for refinement).");
  SLIC_INFO("--nx <N> number of cells in x-direction for sample grid.");
  SLIC_INFO("--ny <N> number of cells in y-direction for sample grid.");
  SLIC_INFO("--nz <N> number of cells in z-direction for sample grid.");
  SLIC_INFO("--help prints this help information.");
}

//------------------------------------------------------------------------------
void parse_args( int argc, char * * argv )
{
  // Defaults
  Arguments.fileName   = "";
  Arguments.maxLevels  = 15;
  Arguments.maxObjects = 5;
  Arguments.nx         = 32;
  Arguments.ny         = 32;
  Arguments.nz         = 32;

  for ( int i=1 ; i < argc ; ++i )
  {

    if ( strcmp( argv[i],"--file")==0 )
    {

      Arguments.fileName = std::string( argv[ ++i ] );

    }
    else if ( strcmp(argv[i],"--maxLevels")==0 )
    {

      Arguments.maxLevels = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i], "--maxObjects" )==0 )
    {

      Arguments.maxObjects = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i],"--nx")==0 )
    {

      Arguments.nx = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i],"--ny")==0 )
    {

      Arguments.ny = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i],"--nz")==0 )
    {

      Arguments.nz = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i],"--help")==0 )
    {
      showHelp();
      exit(0);
    }

  } // END for

  if ( Arguments.fileName == "" )
  {
    SLIC_ERROR( "NO STL input file provided. Provide one with --file" );
  }

}

//------------------------------------------------------------------------------
void write_point( const Point< double, 3 >& pt,
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
void write_triangles( axom::mint::Mesh * mesh,
                      const axom::mint::IndexType * cells,
                      axom::mint::IndexType ncells,
                      const std::string& fileName )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  TriangleMesh * subset = new TriangleMesh(3);
  SLIC_ASSERT( subset->getMeshType() == mesh->getMeshType() );

  axom::mint::IndexType cellIds[3];
  Point< double, 3 > n1;
  Point< double, 3 > n2;
  Point< double, 3 > n3;

  int icount = 0;
  axom::mint::IndexType new_cell[3];

  for ( axom::mint::IndexType i=0 ; i < ncells ; ++i )
  {

    const axom::mint::IndexType cellIdx = cells[ i ];
    mesh->getMeshCell( cellIdx, cellIds );

    mesh->getMeshNode( cellIds[0], n1.data() );
    mesh->getMeshNode( cellIds[1], n2.data() );
    mesh->getMeshNode( cellIds[2], n3.data() );

    subset->addNode( n1[0], n1[1], n1[2] );
    new_cell[0] = icount; ++icount;
    subset->addNode( n2[0], n2[1], n2[2] );
    new_cell[1] = icount; ++icount;
    subset->addNode( n3[0], n3[1], n3[2] );
    new_cell[2] = icount; ++icount;

    subset->addCell( new_cell, MINT_TRIANGLE );
  }

  axom::mint::write_vtk( subset, fileName );

  delete subset;
}

//------------------------------------------------------------------------------
BoundingBox< double,3 > compute_bounds( axom::mint::Mesh * mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  using namespace axom::quest;

  BoundingBox< double,3 > meshBB;
  Point< double,3 > pt;

  for ( int i=0 ; i < mesh->getMeshNumberOfNodes() ; ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  }  // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}


//------------------------------------------------------------------------------
void distance_field( axom::mint::Mesh * surface_mesh,
                     axom::mint::UniformMesh * umesh )
{
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  SLIC_INFO("Max BVH levels: " << Arguments.maxLevels );
  SLIC_INFO("Max Object Threshold: " << Arguments.maxObjects );
  utilities::Timer timer1;
  timer1.start();

  SignedDistance< 3 > signedDistance( surface_mesh,
                                      Arguments.maxObjects,
                                      Arguments.maxLevels );

  timer1.stop();
  SLIC_INFO("Constructed BVH in " << timer1.elapsed() << "s" );

#ifdef AXOM_DEBUG
  // write the bucket tree to a file
  const BVHTree< int, 3> * btree = signedDistance.getBVHTree();
  SLIC_ASSERT( btree != AXOM_NULLPTR );

  btree->writeVtkFile( "bucket-tree.vtk" );

  // mark bucket IDs on surface mesh
//  int * bidx = surface_mesh->addCellField< int >( "BucketID", 1 )->getIntPtr();
  mint::FieldData& CD = surface_mesh->getCellFieldData( );
  const int ncells = surface_mesh->getMeshNumberOfCells( );
  CD.addField( new mint::FieldVariable< int >( "BucketID", ncells ) );
  int* bidx = CD.getField( "BucketID" )->getIntPtr( );
  SLIC_ASSERT( bidx != AXOM_NULLPTR );

  const int numObjects = btree->getNumberOfObjects();
  for ( int i=0 ; i < numObjects ; ++i )
  {

    const int idx = btree->getObjectBucketIndex( i );
    bidx[ i ] = idx;

  } // END for all objects

  axom::mint::write_vtk( surface_mesh, "partitioned_surface_mesh.vtk" );
#endif

  const axom::mint::IndexType nnodes = umesh->getNumberOfNodes();

  mint::FieldData& PD = umesh->getNodeFieldData();
  PD.addField( new mint::FieldVariable< double >("phi", nnodes) );
  PD.addField( new mint::FieldVariable< int >( "nbuckets", nnodes ) );
  PD.addField( new mint::FieldVariable< int >( "ntriangles", nnodes ) );

  double * phi     = PD.getField("phi")->getDoublePtr();
  int * nbuckets   = PD.getField("nbuckets")->getIntPtr();
  int * ntriangles = PD.getField("ntriangles")->getIntPtr();

  SLIC_ASSERT( phi != AXOM_NULLPTR );
  SLIC_ASSERT( nbuckets != AXOM_NULLPTR );
  SLIC_ASSERT( ntriangles != AXOM_NULLPTR );

  utilities::Timer timer2;
  timer2.start();

  for ( axom::mint::IndexType inode=0 ; inode < nnodes ; ++inode )
  {

    Point< double,3 > pt;
    umesh->getMeshNode( inode, pt.data() );

    std::vector< int > buckets;
    std::vector< axom::mint::IndexType > triangles;
    std::vector< axom::mint::IndexType > my_triangles;
    triangles.clear();
    buckets.clear();

    Point< double,3 > closest_pt;
    phi[ inode ] = signedDistance.computeDistance( pt,
                                                   buckets,
                                                   triangles,
                                                   my_triangles,
                                                   closest_pt );

    nbuckets[ inode ]   = static_cast< int >( buckets.size() );
    ntriangles[ inode ] = static_cast< int >( triangles.size() );

#ifdef AXOM_DEBUG
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
int main( int argc, char * * argv )
{
  // STEP 0: Initialize SLIC Environment
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug );

  // Create a more verbose message for this application (only level and message)
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  axom::slic::GenericOutputStream * defaultStream =
    new axom::slic::GenericOutputStream(&std::cout);
  axom::slic::GenericOutputStream * compactStream =
    new axom::slic::GenericOutputStream(&std::cout, slicFormatStr);
  axom::slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);

  // STEP 1: get file from user or use default
  parse_args( argc, argv );

  // STEP 2: read file
  SLIC_INFO( "Reading file: " << Arguments.fileName << "...");
  quest::STLReader * reader = new quest::STLReader();
  reader->setFileName( Arguments.fileName );
  reader->read();
  SLIC_INFO("done");

  // STEP 3: get surface mesh
  axom::mint::Mesh * surface_mesh = new TriangleMesh( 3 );
  reader->getMesh( static_cast<TriangleMesh *>( surface_mesh ) );
  SLIC_INFO("Mesh has "
            << surface_mesh->getMeshNumberOfNodes() << " nodes and "
            << surface_mesh->getMeshNumberOfCells() << " cells.");

  // STEP 4: Delete the reader
  delete reader;
  reader = AXOM_NULLPTR;

  // STEP 5: write vtk file
  axom::mint::write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 6: compute bounds
  BoundingBox< double,3 > meshBB = compute_bounds( surface_mesh);
  SLIC_INFO("Mesh bounding box: " << meshBB );


  BoundingBox< double,3 > queryBounds = meshBB;
  queryBounds.expand( 10 );

  double h[3];
  const Vector< double,3 >& bbDiff = queryBounds.range();
  h[0] = bbDiff[0] / Arguments.nx;
  h[1] = bbDiff[1] / Arguments.ny;
  h[2] = bbDiff[2] / Arguments.nz;
  SLIC_INFO("grid dimensions:" << Arguments.nx << ", "
                               << Arguments.ny << ", "
                               << Arguments.nz );
  SLIC_INFO("grid cell size:(" << h[0] << "," << h[1] << "," << h[2] << ")\n" );

  axom::mint::int64 node_ext[6];
  node_ext[0] = 0;
  node_ext[1] = Arguments.nx;
  node_ext[2] = 0;
  node_ext[3] = Arguments.ny;
  node_ext[4] = 0;
  node_ext[5] = Arguments.nz;

  // STEP 8: Construct uniform mesh
  axom::mint::UniformMesh * umesh =
    new axom::mint::UniformMesh(3, queryBounds.getMin().data(), h, node_ext);

  // STEP 9: Compute the distance field on the uniform mesh
  SLIC_INFO( "computing distance field..." );
  utilities::Timer timer;
  timer.start();
  distance_field( surface_mesh, umesh );
  timer.stop();
  SLIC_INFO( "Time to compute distance field: " << timer.elapsed() );

  // STEP 10: write the uniform mesh
  axom::mint::write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 11: clean up
  delete surface_mesh;
  surface_mesh = AXOM_NULLPTR;

  delete umesh;
  umesh = AXOM_NULLPTR;

  // STEP 12: Finalize SLIC environment
  axom::slic::finalize();
  return 0;
}
