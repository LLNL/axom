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

// Axom utils includes
#include "axom_utils/Utilities.hpp"       // for processAbort()
#include "axom_utils/Timer.hpp"

// Mint includes
#include "mint/config.hpp"                // for compile-time mint
#include "mint/UniformMesh.hpp"           // for mint::UniformMesh
#include "mint/vtk_utils.hpp"             // for mint::write_vtk()

// Quest includes
#include "quest/signed_distance.hpp"      // for Quest's signed distance query

// Slic includes
#include "slic/slic.hpp"                  // for SLIC macros
#include "slic/GenericOutputStream.hpp"   // for GenericOutputStream

#ifdef AXOM_USE_MPI
  #include <mpi.h>                         // for MPI
  #include "slic/SynchronizedStream.hpp"   // for SynchronizedStream
#else
using MPI_Comm = int;
#endif

// C/C++ includes
#include <cstring>    // for strcmp()

/*!
 * \file
 *
 * \brief Simple example that illustrates the use of Quest's C-style Signed
 *  Distance interface from within a C++ application.
 */

// namespace aliases
namespace mint      = axom::mint;
namespace quest     = axom::quest;
namespace slic      = axom::slic;
namespace utilities = axom::utilities;

// Function prototypes
void initialize_logger( );
void finalize_logger( );
void parse_args( int argc, char** argv );
void generate_uniform_box_mesh( mint::UniformMesh*& mesh );
void show_help( );



//------------------------------------------------------------------------------
// GLOBALS
//------------------------------------------------------------------------------
MPI_Comm global_comm;

/*!
 * \brief Holds command-line arguments
 */
static struct
{
  std::string fileName;
  int ndims;
  int maxLevels;
  int maxOccupancy;
  int box_dims[3];
  double box_min[3];
  double box_max[3];
  bool specified_box_min;
  bool specified_box_max;
  bool dump_vtk;
} Arguments;

//------------------------------------------------------------------------------
// PROGRAM MAIN
//------------------------------------------------------------------------------
int main ( int argc, char** argv )
{
  // STEP 0: initialize MPI
#ifdef AXOM_USE_MPI
  MPI_Init( &argc, &argv );
  global_comm = MPI_COMM_WORLD;
#endif

  utilities::Timer timer;

  // STEP 1: initialize the logger
  initialize_logger( );

  // STEP 2: parse command line arguments
  parse_args( argc, argv );

  // STEP 3: initialize the signed distance interface
  SLIC_INFO( "initializing signed distance function..." );
  SLIC_INFO( "input file: " << Arguments.fileName );
  SLIC_INFO( "max_levels=" << Arguments.maxLevels );
  SLIC_INFO( "max_occupancy=" << Arguments.maxOccupancy );
  slic::flushStreams();

  timer.start();
  quest::signed_distance_set_max_levels( Arguments.maxLevels );
  quest::signed_distance_set_max_occupancy( Arguments.maxOccupancy );
  int rc = quest::signed_distance_init( Arguments.fileName, global_comm );
  timer.stop();

  SLIC_ERROR_IF( (rc != 0), "Signed Distance query initialization failed!" );
  SLIC_INFO( "time to initialize: " << timer.elapsed() << "s"  );
  slic::flushStreams();

  // STEP 5: Generate computational mesh
  mint::UniformMesh* mesh = AXOM_NULLPTR;
  generate_uniform_box_mesh( mesh );
  SLIC_ERROR_IF( mesh==AXOM_NULLPTR, "box mesh is null!" );
  double* phi = mesh->createField< double >( "phi", mint::NODE_CENTERED );

  // STEP 6: evaluate the signed distance field on the given mesh
  SLIC_INFO( "evaluating signed distance field on specified box mesh..." );
  slic::flushStreams();

  const mint::IndexType nnodes = mesh->getNumberOfNodes( );

  timer.reset();
  timer.start();
  for ( mint::IndexType inode=0; inode < nnodes; ++inode )
  {
    double pt[ 3 ];
    mesh->getNode( inode, pt );
    phi[ inode ] = quest::signed_distance_evaluate( pt[0], pt[1], pt[2] );
  } // END for all nodes
  timer.stop();
  SLIC_INFO( "time to evaluate: " << timer.elapsed() << "s" );
  slic::flushStreams();

  // STEP 7: vtk output
  if ( Arguments.dump_vtk )
  {
    SLIC_INFO( "writing vtk output" );
    slic::flushStreams();
    mint::write_vtk( mesh, "uniform_mesh.vtk" );
  }

  // STEP 8: finalize
  delete mesh;
  mesh = AXOM_NULLPTR;

  quest::signed_distance_finalize( );
  finalize_logger( );
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}

//------------------------------------------------------------------------------
//  FUNCTION PROTOTYPE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void generate_uniform_box_mesh( mint::UniformMesh*& mesh )
{
  SLIC_ASSERT( mesh == AXOM_NULLPTR );
  SLIC_ASSERT( quest::signed_distance_initialized() );

  double mesh_box_min[ 3 ];
  double mesh_box_max[ 3 ];

  double* lo = AXOM_NULLPTR;
  double* hi = AXOM_NULLPTR;

  if ( Arguments.specified_box_max && Arguments.specified_box_min )
  {
    SLIC_INFO( "using specified box bounds" );
    lo = Arguments.box_min;
    hi = Arguments.box_max;
  }
  else
  {
    SLIC_INFO( "computing mesh bounds..." );
    quest::signed_distance_get_mesh_bounds( mesh_box_min, mesh_box_max );
    lo = mesh_box_min;
    hi = mesh_box_max;
  }


  SLIC_INFO( "box min: [" << lo[0] << " " << lo[1] << " " << lo[2] << "]" );
  SLIC_INFO( "box max: [" << hi[0] << " " << hi[1] << " " << hi[2] << "]" );
  SLIC_INFO( "constructing Uniform Mesh: [" << Arguments.box_dims[0] <<
             " " << Arguments.box_dims[ 1 ] <<
             " " << Arguments.box_dims[ 2 ] << "]" );

  mint::IndexType Ni = static_cast< mint::IndexType >( Arguments.box_dims[0] );
  mint::IndexType Nj = static_cast< mint::IndexType >( Arguments.box_dims[1] );
  mint::IndexType Nk = static_cast< mint::IndexType >( Arguments.box_dims[2] );
  mesh = new mint::UniformMesh( lo, hi, Ni, Nj, Nk );

  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  slic::flushStreams();
}

//------------------------------------------------------------------------------
void parse_args( int argc, char** argv )
{
  // Set defaults
  Arguments.ndims         = 3;
  Arguments.fileName      = "";
  Arguments.maxLevels     = 15;
  Arguments.maxOccupancy  = 5;
  Arguments.box_dims[ 0 ] =
  Arguments.box_dims[ 1 ] =
  Arguments.box_dims[ 2 ] = 32;
  Arguments.specified_box_max = false;
  Arguments.specified_box_min = false;
  Arguments.dump_vtk          = true;

  for ( int i=1 ; i < argc ; ++i )
  {
    if ( strcmp( argv[i], "--file" )==0 || strcmp( argv[i], "-f"  )==0 )
    {
      Arguments.fileName = std::string( argv[ ++i] );
    }
    else if ( strcmp(argv[i], "--dimension")==0 )
    {
      Arguments.ndims = std::atoi( argv[++i] );
    }
    else if ( strcmp( argv[i], "--maxLevels" )==0 )
    {
      Arguments.maxLevels = std::atoi( argv[++i] );
    }
    else if ( strcmp( argv[i], "--maxOccupancy" )==0 )
    {
      Arguments.maxOccupancy = std::atoi( argv[++i] );
    }
    else if ( strcmp( argv[i], "--box-dims" )==0 )
    {
      for ( int j=0 ; j < Arguments.ndims ; ++j )
      {
        Arguments.box_dims[ j ] = std::atoi( argv[++i] );
      } // END for all j
    }
    else if ( strcmp(argv[i], "--box-min")==0 )
    {
      Arguments.specified_box_min = true;
      for ( int j=0 ; j < Arguments.ndims ; ++j )
      {
        Arguments.box_min[ j ] = std::atof( argv[++i] );
      } // END for all j
    }
    else if ( strcmp( argv[i], "--box-max")==0 )
    {
      Arguments.specified_box_max = true;
      for ( int j=0 ; j < Arguments.ndims ; ++j )
      {
        Arguments.box_max[ j ] = std::atof( argv[++i] );
      } // END for all j
    }
    else if ( strcmp( argv[i], "--no-vtk")==0 )
    {
      Arguments.dump_vtk = false;
    }
    else if ( strcmp( argv[i], "--help" )==0 )
    {
      show_help();
      utilities::processAbort();
    }
    else
    {
      SLIC_WARNING( "skipping undefined argument [" << argv[i] << "]..." );
    }

  } // END for all i

  SLIC_ERROR_IF( (Arguments.ndims != 3),
                 "The signed distance is currently only supported in 3-D" );
  SLIC_ERROR_IF( Arguments.fileName.empty(),
                 "Must provide an STL input file. Provide one with --file" );
  SLIC_ERROR_IF( Arguments.specified_box_max == !Arguments.specified_box_min,
                 "Both min/max bounds must be specified.");
  slic::flushStreams();
}

//------------------------------------------------------------------------------
void show_help( )
{
  SLIC_INFO( "Usage:./quest_signed_distance_interface_ex -f <file> [options]");
  SLIC_INFO( "-f, --file <file> specifies the input mesh file" );
  SLIC_INFO( "--dimension <DIM> the problem dimension" );
  SLIC_INFO( "--maxLevels <N> max levels of Subdivision for the BVH" );
  SLIC_INFO( "--maxOccupancy <N> max number of item per BVH bin." );
  SLIC_INFO( "--box-dims <N0> <N1> <N2> the dimensions of the box mesh");
  SLIC_INFO( "--box-min <X0> <Y0> <Z0> the lower corner of the box mesh" );
  SLIC_INFO( "--box-max <XN> <YN> <ZN> the upper cordner of the box mesh" );
  SLIC_INFO( "--no-vtk disables VTK output." );
  SLIC_INFO( "--help prints this help information" );
  slic::flushStreams();
}

//------------------------------------------------------------------------------
void initialize_logger( )
{
  // initialize logger
  slic::initialize( );
  slic::setLoggingMsgLevel( slic::message::Info );

  // setup the logstreams
  std::string fmt              = "";
  slic::LogStream* logStream   = AXOM_NULLPTR;

#ifdef AXOM_USE_MPI
  fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::SynchronizedStream( &std::cout, global_comm, fmt );
#else
  fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream   = new slic::GenericOutputStream( &std::cout, fmt );
#endif

  // register stream objects with the logger
  slic::addStreamToAllMsgLevels( logStream );
}

//------------------------------------------------------------------------------
void finalize_logger( )
{
  slic::flushStreams();
  slic::finalize( );
}
