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

// Mint includes
#include "mint/config.hpp"                // for compile-time mint
#include "mint/UniformMesh.hpp"           // for mint::UniformMesh

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

  // STEP 1: initialize the logger
  initialize_logger( );

  // STEP 2: parse command line arguments
  parse_args( argc, argv );

  // STEP 3: initialize the signed distance interface
  quest::signed_distance_set_max_levels( Arguments.maxLevels );
  quest::signed_distance_set_max_occupancy( Arguments.maxOccupancy );
  int rc = quest::signed_distance_init( Arguments.fileName, global_comm );
  SLIC_ERROR_IF( (rc != 0), "Signed Distance query initialization failed!" );

  // STEP 5: Generate computational mesh
  mint::UniformMesh* mesh = AXOM_NULLPTR;

  // TODO: implement this2

  // STEP 6: evaluate the signed distance field on the given mesh
  // TODO: implement this

  // STEP 7: finalize
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

  for ( int i=1; i < argc; ++i )
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
    else if ( strcmp( argv[i], "--box_dims" )==0 )
    {
      for ( int j=0; j < Arguments.ndims; ++j )
      {
        Arguments.box_dims[ j ] = std::atoi( argv[++i] );
      } // END for all j
    }
    else if ( strcmp(argv[i], "--box-min")==0 )
    {
      Arguments.specified_box_min = true;
      for ( int j=0; j < Arguments.ndims; ++j )
      {
         Arguments.box_min[ j ] = std::atof( argv[++i] );
      } // END for all j
    }
    else if ( strcmp( argv[i], "--box-max")==0 )
    {
      Arguments.specified_box_max = true;
      for ( int j=0; j < Arguments.ndims; ++j )
      {
        Arguments.box_max[ j ] = std::atof( argv[++i] );
      } // END for all j
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
      "The Signed Distance query is currently only supported in 3-D" );
  SLIC_ERROR_IF( Arguments.fileName.empty(),
      "Must provide an STL input file. Provide one with --file" );
  SLIC_ERROR_IF( Arguments.specified_box_max == !Arguments.specified_box_min,
  "To explicitly define the box domain you must specify both min/max bounds.");
}

//------------------------------------------------------------------------------
void show_help( )
{
  SLIC_INFO( "Usage:./quest_signed_distance_interface_ex -f <file> [options]");
  SLIC_INFO( "-f, --file <file> specifies the input mesh file" );
  SLIC_INFO( "--dimension <DIM> the problem dimension" );
  SLIC_INFO( "--maxLevels <N> max levels of Subdivision for the BVH" );
  SLIC_INFO( "--maxOccupancy <N> max number of item per BVH bin." );
  SLIC_INFO( "--box_dims <N0> <N1> <N2> the dimensions of the box mesh");
  SLIC_INFO( "--box_min <X0> <Y0> <Z0> the lower corner of the box mesh" );
  SLIC_INFO( "--box_max <XN> <YN> <ZN> the upper cordner of the box mesh" );
  SLIC_INFO( "--help prints this help information" );
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
  logStream   = new slic::SynchronizedStream( &std::cout, global_comm, fmt );
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
  slic::finalize( );
}
