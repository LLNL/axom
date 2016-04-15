/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "quest/quest.hpp"

// Quest includes
#include "common/CommonTypes.hpp"
#include "quest/STLReader.hpp"
#include "quest/SignedDistance.hpp"
#include "slic/slic.hpp"

#ifdef USE_MPI
#include "quest/PSTLReader.hpp"
#endif

namespace quest
{

// NOTE: supporting just one region/surface for now
static int NDIMS = 3;
static SignedDistance< 3 >* region = ATK_NULLPTR;
static meshtk::Mesh* surface_mesh  = ATK_NULLPTR;

typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;

//------------------------------------------------------------------------------
#ifdef USE_MPI
void initialize( MPI_Comm comm, const std::string& fileName,
                 int ndims, int maxElements, int maxLevels )
{
  SLIC_ASSERT( comm != MPI_COMM_NULL );
  SLIC_ASSERT( region == ATK_NULLPTR );

  NDIMS = ndims;
  SLIC_ASSERT( NDIMS==2 || NDIMS==3 );

  quest::PSTLReader* reader = new quest::PSTLReader( comm );
  reader->setFileName( fileName );
  reader->read();

  surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  region = new SignedDistance< 3 >( surface_mesh, maxElements, maxLevels );
}
#else
//------------------------------------------------------------------------------
void initialize( const std::string& fileName,
                 int ndims, int maxElements, int maxLevels )
{
  SLIC_ASSERT( region == ATK_NULLPTR );

  NDIMS = ndims;
  SLIC_ASSERT( NDIMS==2 || NDIMS==3 );

  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( fileName );
  reader->read();

  surface_mesh = new TriangleMesh( 3 );
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );

  reader->getMesh( static_cast< TriangleMesh* >( surface_mesh ) );
  delete reader;

  region = new SignedDistance< 3 >( surface_mesh, maxElements, maxLevels );
}
#endif

//------------------------------------------------------------------------------
double distance( double x, double y, double z )
{
  SLIC_ASSERT( region != ATK_NULLPTR );

  // TODO: assume 3-D for now
  Point< double,3 > pt;
  pt[0] = x;
  pt[1] = y;
  pt[2] = z;

  return( region->computeDistance( pt ) );
}

//------------------------------------------------------------------------------
void distance( const double* xyz, double* dist, int npoints )
{
  SLIC_ASSERT( xyz != ATK_NULLPTR );
  SLIC_ASSERT( dist != ATK_NULLPTR );
  SLIC_ASSERT( region != ATK_NULLPTR );

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0; i < npoints; ++i ) {

     // TODO: assume 3-D for now
     Point< double,3 > pt;
     pt[0] = xyz[i*3];
     pt[1] = xyz[i*3+1];
     pt[2] = xyz[i*3+2];

     dist[ i ] = region->computeDistance( pt );
  }

}

//------------------------------------------------------------------------------
int inside( double x, double y, double z )
{
  SLIC_ASSERT( region != ATK_NULLPTR );

  // TOOD: assume 3-D for now
  Point< double,3 > pt;
  pt[0] = x;
  pt[1] = y;
  pt[2] = z;

  int sign = -1;
  const quest::SignedDistance<3>::BVHTreeType* tree = region->getBVHTree();
  SLIC_ASSERT( tree != ATK_NULLPTR );

  if ( !tree->contains( pt ) ) {

    sign = 1;

  } else {

    sign = ( region->computeDistance( pt ) < 0.0f )? -1 : 1;

  }

  return( sign );
}

//------------------------------------------------------------------------------
void inside( const double* xyz, int* in, int npoints )
{
  SLIC_ASSERT( xyz != ATK_NULLPTR );
  SLIC_ASSERT( in != ATK_NULLPTR );
  SLIC_ASSERT( region != ATK_NULLPTR );

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0; i < npoints; ++i ) {
     // TODO: assume 3-D for now
     in[ i ] = quest::inside( xyz[i*3], xyz[i*3+1], xyz[i*3+2] );
  }
}

//------------------------------------------------------------------------------
void finalize()
{
  if ( region != ATK_NULLPTR ) {

     delete region;
     region = ATK_NULLPTR;

  } // END if

  if ( surface_mesh != ATK_NULLPTR ) {

     delete surface_mesh;
     surface_mesh = ATK_NULLPTR;
  }

}

} /* end namespace quest */

