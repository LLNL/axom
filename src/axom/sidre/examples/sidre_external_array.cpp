// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core.hpp"  // for Axom types and macros
#include "axom/sidre.hpp" // for sidre
#include "axom/slic.hpp"  // for logging with slic

// MPI includes
#include <mpi.h>

// C/C++ includes
#include <iostream>

// aliases
namespace sidre = axom::sidre;
namespace slic  = axom::slic;

//------------------------------------------------------------------------------
void sidre_write( MPI_Comm comm,
                  const std::string& file,
                  int* data,
                  axom::IndexType numTuples,
                  axom::IndexType numComponents )
{
  SLIC_ASSERT( comm != MPI_COMM_NULL );
  SLIC_ASSERT( data != nullptr );
  SLIC_ASSERT( !file.empty( ) );

  int nranks = -1;
  MPI_Comm_size( comm, &nranks );

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();
  sidre::View*  view = root->createView( "data" );

  sidre::IndexType shape[2];
  shape[ 0 ] = numTuples;
  shape[ 1 ] = numComponents;

  view->setExternalDataPtr( sidre::INT32_ID, 2, shape, data );

// DEBUG
  SLIC_INFO( "Here is the data that is begin dumped:" );
  root->print();
  std::cout << std::endl;
// DEBUG


  // STEP 2: save the array data in to a file
  sidre::IOManager sidre_io( comm );
  sidre_io.write( root, nranks, file, "sidre_hdf5" );
}

//------------------------------------------------------------------------------
void sidre_read( MPI_Comm comm,
                 const std::string& file,
                 int*& data,
                 axom::IndexType& numTuples,
                 axom::IndexType& numComponents )
{
  SLIC_ASSERT( comm != MPI_COMM_NULL );
  SLIC_ASSERT( data == nullptr );
  SLIC_ASSERT( !file.empty( ) );

  int nranks = -1;
  MPI_Comm_size( comm, &nranks );

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  sidre::IOManager sidre_io( comm );
  sidre_io.read( root, file );

// DEBUG
  SLIC_INFO( "Here is the data that was read back:" );
  root->print();
  std::cout << std::endl;
// DEBUG

  SLIC_ASSERT( root->hasChildView("data") );
  sidre::View* view = root->getView( "data" );

  sidre::IndexType shape[2];
  view->getShape( 2, shape );
  numTuples     = shape[ 0 ];
  numComponents = shape[ 1 ];

  axom::IndexType nelems = view->getNumElements();
  SLIC_ASSERT( nelems==(numTuples*numComponents) );

  data = axom::allocate< int >(nelems);
  SLIC_ASSERT( data != nullptr );

  memcpy( data, view->getVoidPtr(), nelems*sizeof(int) );
}

//------------------------------------------------------------------------------
int main ( int argc, char** argv )
{
  MPI_Init( &argc, &argv );
  MPI_Comm problem_comm = MPI_COMM_WORLD;

  slic::UnitTestLogger logger;

  // STEP 0: create some data
  constexpr axom::IndexType NUM_NODES = 10;
  constexpr axom::IndexType DIMENSION = 4;
  constexpr axom::IndexType NSIZE     = NUM_NODES * DIMENSION;

  int* data = axom::allocate< int >( NSIZE );
  SLIC_ASSERT( data != nullptr );

  for ( int i=0; i < NSIZE; ++i )
  {
    data[ i ] = (i+1) * 10;
  }

  // STEP 1: dump the data to a file using sidre
  SLIC_INFO( "Writting data..." );
  sidre_write( problem_comm, "mesh", data, NUM_NODES, DIMENSION );
  SLIC_INFO( "[DONE]" );

  // STEP 2: read the data from a file using sidre
  int* data2              = nullptr;
  axom::IndexType ntuples = -1;
  axom::IndexType ncomp   = -1;

  SLIC_INFO( "Reading data..." );
  sidre_read( problem_comm, "mesh.root", data2, ntuples, ncomp );
  SLIC_INFO( "[DONE]" );

  // STEP 3: check the data
  SLIC_ASSERT( ntuples == NUM_NODES );
  SLIC_ASSERT( ncomp == DIMENSION );
  for ( int i=0; i < NSIZE; ++i )
  {
    SLIC_ASSERT( data[ i ] == data2[ i ] );
  }

  // STEP 4: deallocate
  axom::deallocate( data );
  axom::deallocate( data2 );

  MPI_Finalize();
  return 0;
}


