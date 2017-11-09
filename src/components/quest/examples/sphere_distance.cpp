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
#include "axom/config.hpp"

#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "axom_utils/FileUtilities.hpp"
#include "axom_utils/Timer.hpp"

#include "primal/BVHTree.hpp"
#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Sphere.hpp"
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"

#include "primal/orientation.hpp"
#include "primal/squared_distance.hpp"

#include "quest/STLReader.hpp"
#include "quest/SignedDistance.hpp"

#include "mint/DataTypes.hpp"
#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

#include "slam/Utilities.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

using namespace axom;

using axom::primal::Point;
using axom::primal::Triangle;
using axom::primal::BoundingBox;
using axom::primal::BVHTree;
using axom::primal::Vector;
using axom::primal::Sphere;


static const int NDIMS = 3;
typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

static struct
{
  std::string fileName;
  bool runN2Algorithm;
  double sphere_radius;
  double inflate_factor;
  Point< double,3 > sphere_center;
  int nx;
  int ny;
  int nz;
} Arguments;

//------------------------------------------------------------------------------
// Function Prototypes
//------------------------------------------------------------------------------
void init();
void parse_args( int argc, char** argv );
void read_stl_mesh( TriangleMesh* stl_mesh );
BoundingBox< double,NDIMS > compute_bounds( axom::mint::Mesh* mesh );
void get_uniform_mesh( TriangleMesh* surface_mesh,
                       axom::mint::UniformMesh*& umesh);
BoundingBox< double,NDIMS > getCellBoundingBox(
  int cellIdx, axom::mint::Mesh* surface_mesh );
void computeUsingBucketTree( axom::mint::Mesh* surface_mesh,
                             axom::mint::UniformMesh* umesh );
void n2( axom::mint::Mesh* surface_mesh, axom::mint::UniformMesh* umesh );
void expected_phi(axom::mint::UniformMesh* umesh);
void compute_norms( axom::mint::UniformMesh* umesh,
                    double& l1, double& l2, double& linf );
void showHelp();
void finalize();

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: Initialize SLIC Environment
  init();
  parse_args( argc, argv );

  // STEP 1: read file
  TriangleMesh* surface_mesh = new TriangleMesh( 3 );
  read_stl_mesh( surface_mesh );
  axom::mint::write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 2: get uniform mesh
  axom::mint::UniformMesh* umesh = AXOM_NULLPTR;
  get_uniform_mesh( surface_mesh, umesh );
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  utilities::Timer timer;

  // STEP 3: Run the n^2 algorithm.
  double n2time = 0.0;
  if ( Arguments.runN2Algorithm )
  {

    SLIC_INFO( "Running n^2 algorithm..." );

    timer.start();
    n2(surface_mesh,umesh);
    timer.stop();
    n2time = timer.elapsed();

    SLIC_INFO( "N^2 algorithm run in " << n2time << "s " );
  }

  // STEP 4: Run the BVH algorithm
  SLIC_INFO( "Running BVHTree algorithm..." );
  timer.reset();
  timer.start();
  computeUsingBucketTree(surface_mesh,umesh);
  timer.stop();
  double btreetime = timer.elapsed();
  SLIC_INFO( "BVH algorithm run in " << btreetime << "s " );
  if ( Arguments.runN2Algorithm )
  {
    SLIC_INFO( "Speedup with respect to n^2 algorithm: " <<
               (n2time/btreetime) );
  }

  // STEP 11: compute expected phi and error
  SLIC_INFO( "Compute expected phi & error...." );
  expected_phi( umesh );
  double l1=0.0, l2=0.0, linf=0.0;
  compute_norms( umesh, l1, l2, linf );
  SLIC_INFO( "done." );

  SLIC_INFO( "l1="<< l1 << " l2="<< l2 << " linf=" << linf );

  axom::mint::write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 12: clean up
  delete surface_mesh;
  surface_mesh = AXOM_NULLPTR;

  delete umesh;
  umesh = AXOM_NULLPTR;

  // STEP 13: Finalize SLIC environment
  finalize();
  return 0;
}

//------------------------------------------------------------------------------
void get_uniform_mesh( TriangleMesh* surface_mesh,
                       axom::mint::UniformMesh*& umesh)
{
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );

  BoundingBox< double,NDIMS > meshBounds = compute_bounds(surface_mesh);
  meshBounds.expand(Arguments.inflate_factor);

  double h[3];
  h[0] = ( meshBounds.getMax()[0]-meshBounds.getMin()[0] ) / Arguments.nx;
  h[1] = ( meshBounds.getMax()[1]-meshBounds.getMin()[1] ) / Arguments.ny;
  h[2] = ( meshBounds.getMax()[2]-meshBounds.getMin()[2] ) / Arguments.nz;

  mint::globalIndex ext[6];
  ext[0] = 0;
  ext[1] = Arguments.nx;
  ext[2] = 0;
  ext[3] = Arguments.ny;
  ext[4] = 0;
  ext[5] = Arguments.nz;

  umesh = new axom::mint::UniformMesh(3,meshBounds.getMin().data(),h,ext);
}

//------------------------------------------------------------------------------
void parse_args( int argc, char** argv )
{
  // Defaults
  Arguments.fileName          = "";
  Arguments.runN2Algorithm    = false;
  Arguments.sphere_radius     = 0.5;
  Arguments.sphere_center[0]  = 0.0;
  Arguments.sphere_center[1]  = 0.0;
  Arguments.sphere_center[2]  = 0.0;
  Arguments.inflate_factor    = 2.0;
  Arguments.nx                = 32;
  Arguments.ny                = 32;
  Arguments.nz                = 32;

  for ( int i=1 ; i < argc ; ++i )
  {

    if ( strcmp(argv[i],"--file")==0 )
    {

      Arguments.fileName = std::string( argv[ ++i ] );

    }
    else if ( strcmp(argv[i],"--n2")==0 )
    {

      Arguments.runN2Algorithm = true;

    }
    else if ( strcmp(argv[i],"--radius")==0 )
    {

      Arguments.sphere_radius = std::atof( argv[++i] );

    }
    else if ( strcmp(argv[i],"--center")==0 )
    {

      Arguments.sphere_center[0] = std::atof( argv[++i] );
      Arguments.sphere_center[1] = std::atof( argv[++i] );
      Arguments.sphere_center[2] = std::atof( argv[++i] );

    }
    else if ( strcmp(argv[i],"--ndims")==0 )
    {

      Arguments.nx = std::atoi( argv[++i] );
      Arguments.ny = std::atoi( argv[++i] );
      Arguments.nz = std::atoi( argv[++i] );

    }
    else if ( strcmp(argv[i],"--inflate")==0 )
    {

      Arguments.inflate_factor = std::atof( argv[++i] );

    }
    else if ( strcmp(argv[i],"--help")==0 )
    {
      showHelp();
      exit(0);
    }

  } // END for all arguments

  if ( Arguments.fileName == "" )
  {
    SLIC_ERROR( "No STL input file provided. Provide one with --file" );
  }
}

//------------------------------------------------------------------------------
void showHelp()
{
  SLIC_INFO("Usage: ./quest_sphere_distance_ex --file <sphere.stl> [options]" );
  SLIC_INFO("--file <file> sets the filename of the STL file. Required.");
  SLIC_INFO("--n2 flag that indicates whether to run the N^2 algorithm.");
  SLIC_INFO("--radius <radius> the radius of the sphere. Default 0.5.");
  SLIC_INFO("--center <x y z> the center of the sphere.Default (0.0,0.0,0.0)");
  SLIC_INFO("--inflate <n> factor to inflate bounding box. Default is 2.0");
  SLIC_INFO("--help prints help information.");
}

//------------------------------------------------------------------------------
void read_stl_mesh( TriangleMesh* stl_mesh )
{
  SLIC_INFO("Reading file: " << Arguments.fileName << "...");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( Arguments.fileName );
  reader->read();
  SLIC_INFO("done");

  reader->getMesh( stl_mesh  );

  delete reader;
  reader = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
void init()
{
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug );

  // Create a more verbose message for this application (only level and message)
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  axom::slic::GenericOutputStream* defaultStream =
    new axom::slic::GenericOutputStream(&std::cout);
  axom::slic::GenericOutputStream* compactStream =
    new axom::slic::GenericOutputStream(&std::cout, slicFormatStr);
  axom::slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);

}

//------------------------------------------------------------------------------
void finalize()
{
  axom::slic::finalize();
}

//------------------------------------------------------------------------------
void compute_norms( axom::mint::UniformMesh* umesh,
                    double& l1, double& l2, double& linf )
{
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  const int nnodes = umesh->getNumberOfNodes();

  // STEP 0: grab field pointers
  mint::FieldData& PD = umesh->getNodeFieldData();
  double* phi_computed  = PD.getField( "phi" )->getDoublePtr();
  double* phi_expected  = PD.getField( "expected_phi" )->getDoublePtr();

  SLIC_ASSERT( phi_computed != AXOM_NULLPTR );
  SLIC_ASSERT( phi_expected != AXOM_NULLPTR );

  // STEP 1: add field to store erroraddCell
  double* error = umesh->addNodeField< double >( "error", 1 )->getDoublePtr();
  SLIC_ASSERT( error != AXOM_NULLPTR );

  // STEP 2: loop over nodes and calculate norms
  l1 = l2 = 0.0;
  linf = std::numeric_limits< double >::min( );

  for ( int inode=0 ; inode < nnodes ; ++inode )
  {

    const double dx    = phi_computed[ inode ] - phi_expected[ inode ];
    const double absdx = std::fabs( dx );

    l1 += absdx;
    l2 += dx*dx;

    if ( absdx > linf )
    {
      linf = absdx;
    }

    error[ inode ] = absdx;
  }  // END for all nodes

  l2 = std::sqrt( l2 );
}

//------------------------------------------------------------------------------
void expected_phi(axom::mint::UniformMesh* umesh)
{
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  // STEP 0: Construct sphere with user-supplied center and radius.
  SLIC_INFO("computing expected signed distance field...");
  SLIC_INFO("sphere radius: " << Arguments.sphere_radius );
  SLIC_INFO("sphere center: " << Arguments.sphere_center );

  Sphere< double, 3 > sphere( Arguments.sphere_radius );

  // STEP 1: Add node field to stored exact distance field.
  const int nnodes = umesh->getNumberOfNodes();
  double* phi = umesh->addNodeField< double >( "expected_phi", 1)->getDoublePtr();
  SLIC_ASSERT( phi != AXOM_NULLPTR );

  // STEP 2: loop over uniform mesh nodes and compute distance field
  SLIC_INFO( "Calculating analytic distance field...");
#ifdef AXOM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0 ; i < nnodes ; ++i )
  {

    double pnt[3];
    umesh->getNode( i, pnt );

    phi[ i ]  = sphere.computeSignedDistance( pnt );
  }

  SLIC_INFO("done.");
}

//------------------------------------------------------------------------------
void n2( axom::mint::Mesh* surface_mesh, axom::mint::UniformMesh* umesh )
{
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  // STEP 1: Setup node-centered signed distance field on uniform mesh
  const int nnodes = umesh->getNumberOfNodes();
  double* phi = umesh->addNodeField< double >("n2_phi", 1)->getDoublePtr();
  SLIC_ASSERT( phi != AXOM_NULLPTR );

  // STEP 2: loop over uniform mesh nodes and compute distance field
  SLIC_INFO("Calculating distance field...");

#ifdef AXOM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0 ; i < nnodes ; ++i )
  {

    // get target node
    Point< double,NDIMS > pnt;
    umesh->getNode( i, pnt.data() );
    Point< double,NDIMS > Q = pnt;

    double unsignedMinDistSQ = std::numeric_limits< double >::max();
    int sign = 0;

    const axom::mint::localIndex ncells = surface_mesh->getMeshNumberOfCells();
    for (axom::mint::localIndex j=0; j < ncells; ++j ) 
    {
      // find minimum distance from query point to triangle
      axom::mint::localIndex closest_cell[ 3 ];
      surface_mesh->getMeshCell( j, closest_cell );

      Point< double,NDIMS > a,b,c;
      surface_mesh->getMeshNode( closest_cell[0], a.data() );
      surface_mesh->getMeshNode( closest_cell[1], b.data() );
      surface_mesh->getMeshNode( closest_cell[2], c.data() );
      Triangle< double,3 > T( a,b,c);
      const double sqDist = axom::primal::squared_distance( Q, T );
      if ( sqDist < unsignedMinDistSQ)
      {

        bool negSide =
          axom::primal::orientation(Q,T)==axom::primal::ON_NEGATIVE_SIDE;

        sign = negSide ? -1 : 1;
        unsignedMinDistSQ = sqDist;
      }     // END if
    }   // END for all cells on the surface mesh

    phi[i] = sign * std::sqrt( unsignedMinDistSQ );

  }  // END for all nodes on the uniform mesh
  SLIC_INFO("done." );
}

//------------------------------------------------------------------------------
void computeUsingBucketTree( axom::mint::Mesh* surface_mesh,
                             axom::mint::UniformMesh* umesh )
{
  // Sanity Checks
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );
  SLIC_ASSERT( umesh != AXOM_NULLPTR );

  quest::SignedDistance< NDIMS > signedDistance( surface_mesh, 25, 32 );

  const int nnodes = umesh->getNumberOfNodes();
  double* phi = umesh->addNodeField< double >("phi", 1)->getDoublePtr();
  SLIC_ASSERT( phi != AXOM_NULLPTR );

  for ( int inode=0 ; inode < nnodes ; ++inode )
  {

    Point< double,NDIMS > pt;
    umesh->getMeshNode( inode, pt.data() );

    phi[ inode ] = signedDistance.computeDistance( pt );

  } // END for all nodes

#ifdef AXOM_DEBUG
  // write the bucket tree to a file
  const BVHTree< int,NDIMS >* btree = signedDistance.getBVHTree();
  SLIC_ASSERT( btree != AXOM_NULLPTR );

  btree->writeVtkFile( "bucket-tree.vtk" );

  // mark bucket IDs on surface mesh
  int* bidx = umesh->addCellField< int >("BucketID", 1)->getIntPtr();
  SLIC_ASSERT( bidx != AXOM_NULLPTR );

  const int numObjects = btree->getNumberOfObjects();
  for ( int i=0 ; i < numObjects ; ++i )
  {

    const int idx = btree->getObjectBucketIndex( i );
    bidx[ i ] = idx;

  } // END for all objects

  axom::mint::write_vtk( surface_mesh, "partitioned_surface_mesh.vtk" );
#endif
}

//------------------------------------------------------------------------------
BoundingBox< double,NDIMS > getCellBoundingBox( axom::mint::localIndex cellIdx,
                                                axom::mint::Mesh* surface_mesh )
{
  // Sanity checks
  SLIC_ASSERT( surface_mesh != AXOM_NULLPTR );
  SLIC_ASSERT( cellIdx >= 0 && cellIdx < surface_mesh->getMeshNumberOfCells());

  axom::mint::localIndex cell[3];
  surface_mesh->getMeshCell( cellIdx, cell );

  BoundingBox< double,3 > bb;
  Point< double,3 > pt;

  for ( int i=0 ; i < 3 ; ++i )
  {
    surface_mesh->getMeshNode( cell[i], pt.data() );
    bb.addPoint( pt );
  }  // END for all cell nodes

  return bb;
}

//------------------------------------------------------------------------------
BoundingBox< double,NDIMS > compute_bounds( axom::mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  BoundingBox< double,NDIMS > meshBB;
  Point< double,NDIMS > pt;

  for ( int i=0 ; i < mesh->getMeshNumberOfNodes() ; ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  }  // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}
