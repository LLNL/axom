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

#include "quest/signed_distance.hpp"

#include "quest/QuestHelpers.hpp"
#include "quest/SignedDistance.hpp"

#include "mint/Mesh.hpp"

#include "slic/slic.hpp"

namespace axom
{
namespace quest
{

//------------------------------------------------------------------------------
// INTERNAL DATA STRUCTURES
//------------------------------------------------------------------------------
namespace
{
constexpr int INIT_FAILED = -1;
constexpr int INIT_SUCCESS = 0;

using SignedDistance3D = SignedDistance< 3 >;
using SignedDistance2D = SignedDistance< 2 >;


/*!
 * \brief Holds the options for the SignedDistance query.
 */
static struct parameters_t
{
  int dimension;       /*!< the dimension, 2 or 3 */
  int geometry_type;   /*!< the geometry type, default is WATERTIGHT */
  int max_levels;      /*!< max levels of subdivision for the BVH */
  int max_occupancy;   /*!< max occupancy per BVH bin */
  bool verbose;        /*!< logger verbosity */

  /*!
   * \brief Default Constructor. Sets default values for the parameters.
   */
  parameters_t( ) :
    dimension( 3 ),
    geometry_type( WATERTIGHT ),
    max_levels( 12 ),
    max_occupancy( 5 ),
    verbose( false )
  { }

} Parameters;

// TODO: note the SignedDistance query is currently only supported in 3-D
static SignedDistance3D* s_query    = AXOM_NULLPTR;
static mint::Mesh* s_surface_mesh   = AXOM_NULLPTR;
static bool s_must_delete_mesh      = false;
static bool s_must_finalize_logger  = false;
static bool s_logger_is_initialized = false;

} // end anonymous namespace

//------------------------------------------------------------------------------
// SIGNED DISTANCE QUERY INTERFACE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int signed_distance_init( const std::string& file, MPI_Comm comm )
{
  // STEP 0: initialize logger
  internal::logger_init( s_logger_is_initialized,
                         s_must_finalize_logger,
                         Parameters.verbose,
                         comm );

  SLIC_ASSERT( s_query == AXOM_NULLPTR );

  if ( Parameters.dimension != 3 )
  {
    SLIC_WARNING("the SignedDistance Query is currently only supported in 3D");
    return INIT_FAILED;
  }

  // STEP 0: read the STL mesh
  int rc = internal::read_mesh( file, s_surface_mesh, comm );
  if ( rc != 0 )
  {
    SLIC_WARNING( "reading mesh from [" << file << "] failed!" );
    return INIT_FAILED;
  }

  // STEP 1: initialized the signed distance query
  s_must_delete_mesh = true;
  rc = signed_distance_init( s_surface_mesh, comm );
  return rc;
}

//------------------------------------------------------------------------------
int signed_distance_init( const mint::Mesh* m, MPI_Comm comm )
{
  internal::logger_init( s_logger_is_initialized,
                         s_must_finalize_logger,
                         Parameters.verbose,
                         comm );

  SLIC_ERROR_IF( signed_distance_initialized(),
                 "signed distance query has already been initialized!" );
  SLIC_ERROR_IF( m->getDimension() != 3,
                 "signed distance query currently only support 3-D meshes" );
  SLIC_ERROR_IF(
    m->getMeshType() != mint::UNSTRUCTURED_MESH,
    "signed distance query currently only supports unstructured meshes" );
  SLIC_ERROR_IF(
    m->hasMixedCellTypes()==true,
    "signed distance query does not support meshes with mixed shape topology" );
  SLIC_ERROR_IF(
    m->getCellType() != mint::TRIANGLE,
    "signed distance currently only support 3D triangular surface meshes" );

  if ( s_surface_mesh != m )
  {
    SLIC_ASSERT( s_surface_mesh == AXOM_NULLPTR );
    s_surface_mesh     = const_cast< mint::Mesh* >( m );
    s_must_delete_mesh = false;
  }

  s_query = new SignedDistance3D( s_surface_mesh,
                                  Parameters.max_occupancy,
                                  Parameters.max_levels );

  return INIT_SUCCESS;
}

//------------------------------------------------------------------------------
bool signed_distance_initialized( )
{
  return ( s_query != AXOM_NULLPTR );
}

//------------------------------------------------------------------------------
void signed_distance_get_mesh_bounds( double* lo, double* hi )
{
  SLIC_ERROR_IF( !signed_distance_initialized(),
                 "signed distance query must be initialized prior to" <<
                 "calling get_mesh_bounds()" );
  SLIC_ERROR_IF( lo==AXOM_NULLPTR, "supplied buffer is null" );
  SLIC_ERROR_IF( hi==AXOM_NULLPTR, "supplied buffer is null" );

  internal::compute_mesh_bounds( s_surface_mesh, lo, hi );
}

//------------------------------------------------------------------------------
void signed_distance_set_dimension( int dim )
{
  SLIC_ERROR_IF( dim != 3, "The signed distance query only support 3D" );
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.dimension = dim;
}

//------------------------------------------------------------------------------
void signed_distance_set_geometry( int type )
{
  SLIC_ERROR_IF( type != WATERTIGHT, "invalid geometry type" );
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.geometry_type = type;
}

//------------------------------------------------------------------------------
void signed_distance_set_max_levels( int maxLevels )
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.max_levels = maxLevels;
}

//------------------------------------------------------------------------------
void signed_distance_set_max_occupancy( int threshold )
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.max_occupancy = threshold;
}

//------------------------------------------------------------------------------
void signed_distance_set_verbose( bool status )
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.verbose = status;
}

//------------------------------------------------------------------------------
double signed_distance_evaluate( double x, double y, double z )
{
  SLIC_ERROR_IF(
    !signed_distance_initialized(),
    "signed distance query must be initialized prior to calling evaluate()!" );

  double phi = s_query->computeDistance( x, y, z );
  return ( phi );
}

//------------------------------------------------------------------------------
void signed_distance_evaluate( const double* x,
                               const double* y,
                               const double* z,
                               int npoints,
                               double* phi )
{
  SLIC_ERROR_IF(
    !signed_distance_initialized(),
    "signed distance query must be initialized prior to calling evaluate()!" );
  SLIC_ERROR_IF( x == AXOM_NULLPTR, "x-coords array is null" );
  SLIC_ERROR_IF( y == AXOM_NULLPTR, "y-coords array is null" );
  SLIC_ERROR_IF( z == AXOM_NULLPTR, "z-coords array is null" );
  SLIC_ERROR_IF ( phi == AXOM_NULLPTR, "output phi array is null" );

#ifdef AXOM_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for ( int i=0 ; i < npoints ; ++i )
  {
    phi[ i ] = s_query->computeDistance( x[i], y[i], z[i] );
  }
}

//------------------------------------------------------------------------------
void signed_distance_finalize( )
{
  if ( s_query != AXOM_NULLPTR )
  {
    delete s_query;
    s_query = AXOM_NULLPTR;
  }

  if ( s_surface_mesh != AXOM_NULLPTR  && s_must_delete_mesh )
  {
    delete s_surface_mesh;
  }

  s_surface_mesh = AXOM_NULLPTR;

  SLIC_ASSERT( !signed_distance_initialized() );
  internal::logger_finalize( s_must_finalize_logger);
}

} // end namespace quest
} // end namespace axom
