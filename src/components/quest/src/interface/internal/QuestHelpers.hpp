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

#ifndef QUEST_HELPERS_HPP_
#define QUEST_HELPERS_HPP_

// Axom includes
#include "axom/config.hpp"           // for compile-time definitions

// Mint includes
#include "mint/Mesh.hpp"             // for mint::Mesh
#include "mint/UnstructuredMesh.hpp" // for mint::UnstructuredMesh

// Slic includes
#include "slic/slic.hpp"  // for SLIC macros
#include "slic/GenericOutputStream.hpp"
#if defined(AXOM_USE_MPI) && defined(AXOM_USE_LUMBERJACK)
  #include "slic/LumberjackStream.hpp"
#elif defined(AXOM_USE_MPI) && !defined(AXOM_USE_LUMBERJACK)
  #include "slic/SynchronizedStream.hpp"
#endif

// Quest includes
#include "quest/mpicomm_wrapper.hpp"
#include "quest/STLReader.hpp"
#ifdef AXOM_USE_MPI
#include "quest/PSTLReader.hpp"
#endif

// C/C++ includes
#include <string> // for C++ string

/*!
 * \file
 *
 * \brief Helper methods that can be used across the different Quest queries.
 */
namespace axom
{
namespace quest
{
namespace internal
{

constexpr int READ_FAILED  = -1;
constexpr int READ_SUCCESS = 0;

/// \name Mesh I/O methods
/// @{

/*!
 * \brief Reads in the surface mesh from the specified file.
 *
 * \param [in] file the file consisting of the surface
 * \param [out] m user-supplied pointer to point to the mesh object.
 * \param [in] comm the MPI communicator, only applicable when MPI is available.
 *
 * \note This method currently expects the surface mesh to be given in STL
 *  format.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function.
 *
 * \return status set to zero on success, or to a non-zero value otherwise.
 *
 * \pre m == AXOM_NULLPTR
 * \pre !file.empty()
 *
 * \post m != AXOM_NULLPTR
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::TRIANGLE
 *
 * \see STLReader
 * \see PSTLReader
 */
int read_mesh( const std::string& file,
               mint::Mesh*& m,
               MPI_Comm comm=MPI_COMM_SELF )
{
  // NOTE: STL meshes are always 3D
  constexpr int DIMENSION = 3;
  using TriangleMesh      = mint::UnstructuredMesh< mint::SINGLE_SHAPE >;

  // STEP 0: check input mesh pointer
  if ( m != AXOM_NULLPTR )
  {
    SLIC_WARNING( "supplied mesh pointer is not null!" );
    return READ_FAILED;
  }

  // STEP 1: allocate output mesh object
  m = new TriangleMesh( DIMENSION, mint::TRIANGLE );

  // STEP 2: allocate reader
  quest::STLReader* reader = AXOM_NULLPTR;
#ifdef AXOM_USE_MPI
  reader = new quest::PSTLReader( comm );
#else
  static_cast< void >( comm );        // to silence compiler warnings
  reader = new quest::STLReader();
#endif

  // STEP 3: read the mesh from the STL file
  reader->setFileName( file );
  int rc = reader->read( );
  if ( rc == READ_SUCCESS )
  {
    reader->getMesh( static_cast< TriangleMesh* >( m )  );
  }
  else
  {
    SLIC_WARNING( "reading STL file failed, setting mesh to NULL" );
    delete m;
    m = AXOM_NULLPTR;
  }

  // STEP 4: delete the reader
  delete reader;
  reader = AXOM_NULLPTR;

  return rc;
}

/// @}

/// \name Mesh Helper Methods
/// @{

/*!
 * \brief Computes the bounds of the given mesh.
 *
 * \param [in] mesh pointer to the mesh whose bounds will be computed.
 * \param [out] lo buffer to store the lower bound mesh coordinates
 * \param [out] hi buffer to store the upper bound mesh coordinates
 *
 * \pre mesh != AXOM_NULLPTR
 * \pre lo != AXOM_NULLPTR
 * \pre hi != AXOM_NULLPTR
 * \pre hi & lo must point to buffers that are at least N long, where N
 *  corresponds to the mesh dimension.
 */
void compute_mesh_bounds( const mint::Mesh* mesh, double* lo, double* hi )
{

  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  SLIC_ASSERT( lo != AXOM_NULLPTR );
  SLIC_ASSERT( hi != AXOM_NULLPTR );

  const int ndims = mesh->getDimension();

  // STEP 0: initialize lo,hi
  for ( int i=0 ; i < ndims ; ++i )
  {
    lo[ i ] = std::numeric_limits< double >::max();
    hi[ i ] = std::numeric_limits< double >::lowest();
  } // END for all dimensions

  // STEP 1: compute lo,hi
  double pt[ 3 ];
  const mint::IndexType numNodes = mesh->getNumberOfNodes();
  for ( mint::IndexType inode=0 ; inode < numNodes ; ++inode )
  {
    mesh->getNode( inode, pt );
    for ( int i=0 ; i < ndims ; ++i )
    {
      lo[ i ] = ( pt[ i ] < lo[ i ] ) ? pt[ i ] : lo[ i ];
      hi[ i ] = ( pt[ i ] > hi[ i ] ) ? pt[ i ] : hi[ i ];
    } // END for all dimensions

  } // END for all nodes

}
/// @}

/// \name Logger Initialize/Finalize Methods
/// @{

/*!
 * \brief Helper method to initialize the Slic logger if needed.
 *
 * \param [in,out] isInitialized indicates if Slic is already initialized.
 * \param [out] mustFinalize inidicates if the caller would be responsible
 *  for finalizing the Slic logger.
 * \param [in] verbose flag to control the verbosity
 * \param [in] comm the MPI communicator (applicable when compiled with MPI)
 *
 * \note If Slic is not already initialized, this method will initialize the
 *  Slic Logging environment and set the `isInitialized` flag to true.
 *
 * \note The 'verbose' flag is only applicable when the Slic logging environment
 *  is not already initialized by the calling application. In that case, when
 *  'verbose' is true, all messages will get logged to the console, including,
 *  Info and debug messages. Otherwise, if 'false', only errors will be printed
 *  out.
 *
 *  \see logger_finalize
 */
void logger_init( bool& isInitialized,
                  bool& mustFinalize,
                  bool verbose,
                  MPI_Comm comm )
{
  if ( isInitialized )
  {
    // Query has already initialized the logger
    return;
  }

  if ( slic::isInitialized() )
  {
    // logger is initialized by an application, the application will finalize
    isInitialized = true;
    mustFinalize  = false;
    return;
  }

  // The SignedDistance Query must initialize the Slic logger and is then
  // also responsible for finalizing it when done
  isInitialized = true;
  mustFinalize  = true;
  slic::initialize();

  slic::LogStream* ls = AXOM_NULLPTR;
  std::string msgfmt  = "[<LEVEL>]: <MESSAGE>\n";

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_LUMBERJACK)
  msgfmt.insert(0,"[<RANK>]", 8 );
  constexpr int RLIMIT = 8;
  ls = new slic::LumberjackStream( &std::cout, comm, RLIMIT, msgfmt );
#elif defined(AXOM_USE_MPI) && !defined(AXOM_USE_LUMBERJACK)
  msgfmt.insert( 0, "[<RANK>]", 8 );
  ls = new slic::SynchronizedStream( &std::cout, comm, msgfmt );
#else
  static_cast< void >( comm );        // to silence compiler warnings
  ls = new slic::GenericOutputStream( &std::cout, msgfmt );
#endif

  slic::addStreamToAllMsgLevels( ls );
  slic::setLoggingMsgLevel(
    ( verbose ) ? slic::message::Info : slic::message::Error );
}

/*!
 * \brief Finalizes the Slic logger (if needed)
 *
 * \param [in] mustFinalize flag that indicates whether the query is responsible
 *  for finalizing the Slic logger.
 *
 * \see logger_init
 */
void logger_finalize( bool mustFinalize )
{

  if ( mustFinalize )
  {
    slic::finalize();
  }

}
/// @}

} /* end namespace internal */
} /* end namespace quest    */
} /* end namespace axom     */

#endif /* QUEST_HELPERS_HPP_ */
