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

#include "slic/slic.hpp"            // for slic macros

#include "axom_utils/Utilities.hpp" // for utilities::max()

#include "mint/config.hpp"          // for mint::IndexType
#include "mint/MeshCoordinates.hpp" // for mint::MeshCoordinates

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"          // for sidre::Group, sidre::View
#endif

#include "gtest/gtest.h" // for gtest macros

// namespace aliases
namespace mint      = axom::mint;
namespace utilities = axom::utilities;
namespace sidre     = axom::sidre;

// constants used in tests
constexpr mint::IndexType DEFAULT_CAPACITY    = 100;
constexpr mint::IndexType ZERO_NUM_NODES      = 0;
constexpr mint::IndexType SMALL_NUM_NODES     = 4;
constexpr mint::IndexType LARGE_NUM_NODES     = 256;
constexpr mint::IndexType IGNORE_CAPACITY     = -1;
constexpr mint::IndexType SMALL_NODE_CAPACITY = 5;
constexpr mint::IndexType LARGE_NODE_CAPACITY = 256;

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace detail
{

/*!
 * \brief Checks that the values of two arrays are the same.
 * \param [in] actual pointer to the buffer of values to check
 * \param [in] expected pointer to the buffer consisting the expected values
 * \param [in] N the length of the buffer.
 *
 * \pre actual != AXOM_NULLPTR
 * \pre expected != AXOM_NULLPTR
 * \pre N > 0
 *
 */
void check_array_values( const double* actual,
                         const double* expected,
                         mint::IndexType N )
{
  SLIC_ASSERT( actual != AXOM_NULLPTR );
  SLIC_ASSERT( expected != AXOM_NULLPTR );
  SLIC_ASSERT( N > 0 );

  for ( mint::IndexType i=0; i < N; ++i )
  {
    EXPECT_DOUBLE_EQ( actual[ i ], expected[ i ] );
  }
}

#ifdef MINT_USE_SIDRE
/*!
 * \brief Create a valid sidre datastore to use in tests
 *
 * \param [in,out] ds pointer to the datastore
 * \param [in] dimension the requested dimension
 *
 * \pre 1 <= dimension <= 3
 * \note Caller must delete the datastore instance and group.
 */
void create_sidre_data( sidre::DataStore& ds, int dimension )
{
  SLIC_ASSERT( (dimension >= 1) && (dimension <= 3) );

  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double* ptrs[ 3 ]           = { x, y, z };

  sidre::Group* gp = ds.getRoot();
  SLIC_ASSERT( gp != AXOM_NULLPTR );

  gp->createView( "type" )->setString( "explicit" );
  sidre::Group* values = gp->createGroup( "values" );

  const char* coord_names[ 3 ] = { "x", "y", "z" };

  for ( int idim=0; idim < dimension; ++idim )
  {
    const char* name = coord_names[ idim ];
    sidre::View* coord_view = values->createView( std::string( name ) );

    // NOTE: even though the array goes out-of-scope here, the data
    // remains persistent in sidre
    mint::Array< double > coord_array (
        coord_view, SMALL_NUM_NODES, 1, SMALL_NUM_NODES );

    coord_array.set( ptrs[ idim ], SMALL_NUM_NODES, 0 );
  } // END for all dimensions
}
#endif

/*!
 * \brief Tests that the coordinate arrays are not AXOM_NULLPTR
 * \param [in] coords const pointer to the MeshCoordinates object.
 */
void check_coordinate_arrays( const mint::MeshCoordinates* coords )
{
  SLIC_ASSERT( coords != AXOM_NULLPTR );

  const int ndims = coords->dimension();
  for ( int i=0; i < ndims; ++i )
  {
    const double* coordsptr = coords->getCoordinateArray( i );
    EXPECT_TRUE( coordsptr != AXOM_NULLPTR );
  }
}

/*!
 * \brief Tests construction of a MeshCoordinates object of specified dimension.
 * \param [in] dimension the mesh dimension
 */
void check_constructor( int dimension )
{
  EXPECT_TRUE( dimension >= 1 && dimension <= 3 );

  // STEP 0: construct MeshCoordinates object
  mint::MeshCoordinates coords( dimension );

  // STEP 1: check post-conditions
  EXPECT_EQ( dimension, coords.dimension() );
  EXPECT_EQ( 0, coords.numNodes() );
  EXPECT_TRUE( coords.capacity() > 0 );
  EXPECT_EQ( coords.capacity(), DEFAULT_CAPACITY );
  EXPECT_TRUE( coords.empty() );

  // STEP 2: check coordinate arrays
  check_coordinate_arrays( &coords );
}

/*!
 * \brief Tests construction of a MeshCoordinates instance of specified
 *  dimension, number of nodes and optionally initial max capacity.
 *
 * \param [in] dimension the mesh dimension
 * \param [in] numNodes the number of nodes in the mesh
 * \param [in] capacity max initial capacity (optional)
 */
void check_constructor( int dimension,
                       mint::IndexType numNodes,
                       mint::IndexType capacity=IGNORE_CAPACITY )
{
  EXPECT_TRUE( dimension >= 1 && dimension <= 3 );

  // STEP 0: construct the MeshCoordinates object
  mint::MeshCoordinates *coords = AXOM_NULLPTR;
  if ( capacity == IGNORE_CAPACITY )
  {
    coords = new mint::MeshCoordinates( dimension, numNodes );
  }
  else
  {
    coords = new mint::MeshCoordinates( dimension, numNodes, capacity );
  }
  EXPECT_TRUE( coords != AXOM_NULLPTR );

  // STEP 1: check post-conditions
  EXPECT_EQ( dimension, coords->dimension() );
  EXPECT_EQ( numNodes, coords->numNodes() );
  EXPECT_TRUE( coords->numNodes() <= coords->capacity() );

  if ( numNodes==ZERO_NUM_NODES )
  {
    EXPECT_TRUE( coords->empty() );
  }

  // STEP 2: check actual capacity
  const mint::IndexType actual_capacity = coords->capacity();

  const double ratio = mint::Array< double >::DEFAULT_RESIZE_RATIO;
  const mint::IndexType expected_computed_capacity =
      utilities::max( DEFAULT_CAPACITY,
                      static_cast< mint::IndexType >( numNodes*ratio+0.5 ) );

  if ( capacity==IGNORE_CAPACITY )
  {
    EXPECT_EQ( actual_capacity, expected_computed_capacity );
  }
  else
  {
    EXPECT_EQ( actual_capacity, capacity );
  }

  // STEP 3: check coordinate arrays
  check_coordinate_arrays( coords );

  // STEP 4: delete MeshCoordinates object
  delete coords;
  coords = AXOM_NULLPTR;
}

/*!
 * \brief Checks the append() method on the MeshCoordinates object
 * \param [in,out] mc the MeshCoordinates object to test with
 * \pre mc != AXOM_NULLPTR
 */
void check_append( mint::MeshCoordinates* mc )
{
  EXPECT_TRUE( mc != AXOM_NULLPTR );

  constexpr int NUM_APPENDS   = 2;
  constexpr double TEST_VALUE = 7;

  // Construct a test node to append
  const int ndims = mc->dimension();
  mint::Array< double > new_node( ndims, 1, ndims );
  for ( int i=0; i < ndims; ++i )
  {
    new_node( i ) = TEST_VALUE + i;
  }

  for ( int iter=0; iter < NUM_APPENDS; ++iter )
  {
    // Append a new node
    mint::IndexType currentNumNodes = mc->numNodes();
    mint::IndexType idx = mc->append( new_node.getData() );
    EXPECT_EQ( idx, currentNumNodes );
    EXPECT_EQ( mc->numNodes(), currentNumNodes+1 );

    // Ensure the data on the new node is what we expect
    for ( int i=0; i < ndims; ++i )
    {
      EXPECT_DOUBLE_EQ( new_node(i), mc->getCoordinate( idx, i ) );
    }

    // Ensure invariant holds after the append
    EXPECT_TRUE( mc->numNodes() <= mc->capacity() );

    // shrink the buffer so that the next append will trigger a realloc
    mc->shrink();
    currentNumNodes = mc->numNodes();
    mint::IndexType currentCapacity = mc->capacity();
    EXPECT_EQ( currentNumNodes, currentCapacity);
  }

}

/*!
 * \brief Test set/get a node from the given MeshCoordinates object.
 * \param [in] mc the mesh coordinates object to test
 * \pre mc != AXOM_NULLPTR
 */
void check_set_and_get( mint::MeshCoordinates* mc )
{
  EXPECT_TRUE( mc != AXOM_NULLPTR );

  constexpr double TEST_VALUE         = 7;
  constexpr mint::IndexType targetIdx = 3;

  const int nnodes = mc->numNodes( );
  const int ndims  = mc->dimension( );

  EXPECT_TRUE( targetIdx < nnodes );

  mint::Array< double > node( ndims, 1, ndims );
  for ( int i=0; i < ndims; ++i )
  {
    node( i ) = TEST_VALUE + 0.5 ;
  }

  mc->set( targetIdx, node.getData() );

  for ( int i=0; i < ndims; ++i )
  {
    EXPECT_DOUBLE_EQ( mc->getCoordinate( targetIdx, i ), node( i ) );
  }

}

} /* end detail namespace*/

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST( mint_mesh_coordinates, append )
{
  constexpr int NDIMS = 3;

  for ( int dim=1; dim <= NDIMS; ++dim )
  {
    // Construct an empty MeshCoordinates object
    mint::MeshCoordinates mc1( dim );
    detail::check_coordinate_arrays( &mc1 );
    detail::check_append( &mc1 );

#ifdef MINT_USE_SIDRE
    // Construct a MeshCoordinates object from sidre
    sidre::DataStore ds;
    detail::create_sidre_data( ds, dim );
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    mint::MeshCoordinates mc2( coords_group );
    detail::check_coordinate_arrays( &mc2 );
    detail::check_append( &mc2 );
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, set_and_get )
{
  constexpr int NDIMS = 3;

  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };

  for ( int dim=1; dim <= NDIMS; ++dim )
  {

    // MeshCoordinates object with native store
    mint::MeshCoordinates mc1( dim, SMALL_NUM_NODES );
    detail::check_set_and_get( &mc1 );

#ifdef MINT_USE_SIDRE
    // MeshCoordinates object from sidre store
    sidre::DataStore ds;
    mint::MeshCoordinates mc2(
                          ds.getRoot(), dim, SMALL_NUM_NODES, SMALL_NUM_NODES );
    detail::check_set_and_get( &mc2 );
#endif

    // MeshCoordinates object tied to external buffers
    if ( dim==1 )
    {
      mint::MeshCoordinates mc3( SMALL_NUM_NODES, x);
      detail::check_set_and_get( &mc3 );
    }
    else if ( dim==2 )
    {
      mint::MeshCoordinates mc3( SMALL_NUM_NODES, x, y );
      detail::check_set_and_get( &mc3 );
    }
    else
    {
      EXPECT_EQ( dim, 3 );
      mint::MeshCoordinates mc3( SMALL_NUM_NODES, x, y, z );
      detail::check_set_and_get( &mc3 );
    }

  }

}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, reserve )
{
  constexpr int NDIMS = 3;

  for ( int dim=1; dim <= NDIMS; ++dim )
  {
    // Construct an empty MeshCoordinates object
    mint::MeshCoordinates mc1( dim );
    mint::IndexType numNodes1 = mc1.numNodes();

    mc1.reserve( LARGE_NODE_CAPACITY );

    EXPECT_EQ( mc1.numNodes(), numNodes1);
    EXPECT_EQ( mc1.capacity(), LARGE_NODE_CAPACITY );

#ifdef MINT_USE_SIDRE
    // Construct a MeshCoordinates object from sidre
    sidre::DataStore ds;
    detail::create_sidre_data( ds, dim );
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    mint::MeshCoordinates mc2( coords_group );
    mint::IndexType numNodes2 = mc2.numNodes();

    mc2.reserve( LARGE_NODE_CAPACITY );

    EXPECT_EQ( mc2.numNodes(), numNodes2 );
    EXPECT_EQ( mc2.capacity(), LARGE_NODE_CAPACITY );
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coorindates, resize )
{
  constexpr int NDIMS = 3;

  for ( int dim=1; dim <= NDIMS; ++dim )
  {
    // Construct an empty MeshCoordinates object
    mint::MeshCoordinates mc1( dim );
    mc1.resize( LARGE_NUM_NODES );
    EXPECT_EQ( mc1.numNodes(), LARGE_NUM_NODES );
    EXPECT_TRUE( mc1.numNodes() <= mc1.capacity() );

#ifdef MINT_USE_SIDRE
    // Construct a MeshCoordinates object from sidre
    sidre::DataStore ds;
    detail::create_sidre_data( ds, dim );
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    mint::MeshCoordinates mc2( coords_group );
    mc2.resize( LARGE_NUM_NODES );
    EXPECT_EQ( mc2.numNodes(), LARGE_NUM_NODES );
    EXPECT_TRUE( mc2.numNodes() <= mc2.capacity() );
#endif
  }

}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, shrink )
{
  constexpr int NDIMS = 3;

  for ( int dim=1; dim <= NDIMS; ++dim )
  {
    // test shrink() when constructed using native store
    mint::MeshCoordinates mc1( dim, SMALL_NUM_NODES );
    EXPECT_EQ( mc1.numNodes(), SMALL_NUM_NODES );
    EXPECT_TRUE( mc1.numNodes() < mc1.capacity() );

    mc1.shrink();

    EXPECT_EQ( mc1.numNodes(), mc1.capacity() );

#ifdef MINT_USE_SIDRE
    // test shrink() when constructed using sidre
    sidre::DataStore ds;
    detail::create_sidre_data( ds, dim );
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    mint::MeshCoordinates mc2( coords_group );

    mc2.shrink();

    EXPECT_EQ( mc2.numNodes(), mc2.capacity() );
#endif
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, change_resize_ratio )
{
  constexpr int NDIMS               = 3;
  constexpr double DEFAULT_RESIZE_RATIO = 
                                    mint::Array< double >::DEFAULT_RESIZE_RATIO;
  constexpr double NEW_RESIZE_RATIO = 2.5;

  mint::MeshCoordinates mc( NDIMS );
  EXPECT_DOUBLE_EQ( mc.getResizeRatio(), DEFAULT_RESIZE_RATIO );

  mc.setResizeRatio( NEW_RESIZE_RATIO );
  EXPECT_DOUBLE_EQ( mc.getResizeRatio(), NEW_RESIZE_RATIO );
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, native_constructors )
{
  for ( int dim=1; dim <= 3; ++dim )
  {
    detail::check_constructor( dim );

    detail::check_constructor( dim, ZERO_NUM_NODES );
    detail::check_constructor( dim, ZERO_NUM_NODES, SMALL_NODE_CAPACITY );
    detail::check_constructor( dim, ZERO_NUM_NODES, LARGE_NODE_CAPACITY );

    detail::check_constructor( dim, SMALL_NUM_NODES );
    detail::check_constructor( dim, SMALL_NUM_NODES, SMALL_NODE_CAPACITY );
    detail::check_constructor( dim, SMALL_NUM_NODES, LARGE_NODE_CAPACITY );

    detail::check_constructor( dim, LARGE_NUM_NODES );
    detail::check_constructor( dim, LARGE_NUM_NODES, LARGE_NODE_CAPACITY );
  }
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, external_constructor )
{
  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };

  // Test 1-D
  {
    constexpr int EXPECTED_DIMENSION = 1;
    mint::MeshCoordinates coords( SMALL_NUM_NODES, x );
    detail::check_coordinate_arrays( &coords );
    EXPECT_EQ( coords.dimension(), EXPECTED_DIMENSION );
    EXPECT_EQ( coords.numNodes(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.capacity(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.getCoordinateArray( mint::X_COORDINATE ), x );
  }
  for ( int i = 0; i < SMALL_NUM_NODES; ++i )
  {
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
  }

  // Test 2-D
  {
    constexpr int EXPECTED_DIMENSION = 2;
    mint::MeshCoordinates coords( SMALL_NUM_NODES, x, y );
    detail::check_coordinate_arrays( &coords );
    EXPECT_EQ( coords.dimension(), EXPECTED_DIMENSION );
    EXPECT_EQ( coords.numNodes(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.capacity(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.getCoordinateArray( mint::X_COORDINATE ), x );
    EXPECT_EQ( coords.getCoordinateArray( mint::Y_COORDINATE ), y );
  }
  for ( int i = 0; i < SMALL_NUM_NODES; ++i )
  {
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
  }

  // Test 3-D
  {
    constexpr int EXPECTED_DIMENSION = 3;
    mint::MeshCoordinates coords( SMALL_NUM_NODES, x, y, z );
    detail::check_coordinate_arrays( &coords );
    EXPECT_EQ( coords.dimension(), EXPECTED_DIMENSION );
    EXPECT_EQ( coords.numNodes(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.capacity(), SMALL_NUM_NODES );
    EXPECT_EQ( coords.getCoordinateArray( mint::X_COORDINATE ), x );
    EXPECT_EQ( coords.getCoordinateArray( mint::Y_COORDINATE ), y );
    EXPECT_EQ( coords.getCoordinateArray( mint::Z_COORDINATE ), z );
  }
  for ( int i = 0; i < SMALL_NUM_NODES; ++i )
  {
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
    EXPECT_DOUBLE_EQ( x[i], i + 1 );
  }

}

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

TEST( mint_mesh_coordinates, sidre_pull_constructor )
{
  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double* expected_data[ 3 ]  = { x, y, z };
  const char* coord_names[3]  = { "x", "y", "z" };

  for ( int dim=1; dim <= 3; ++dim )
  {
    // STEP 0: create a test datastore with a coordinates groupp
    sidre::DataStore ds;
    detail::create_sidre_data( ds, dim );
    sidre::Group* coords_group = ds.getRoot();
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    // STEP 1: create a MeshCoordinates instance within a scope from
    // the sidre group.
    // BEGIN SCOPE
    {
      mint::MeshCoordinates coords( coords_group );
      EXPECT_EQ( coords.dimension(), dim );
      detail::check_coordinate_arrays( &coords );

      EXPECT_FALSE( coords.empty() );
      EXPECT_EQ( coords.numNodes(), SMALL_NUM_NODES );
      EXPECT_TRUE( coords.numNodes() <= coords.capacity() );

      detail::check_coordinate_arrays( &coords );

      for ( int j=0; j < dim; ++j )
      {
        detail::check_array_values( coords.getCoordinateArray( j ),
                                    expected_data[ j ],
                                    SMALL_NUM_NODES );
      }

    }
    // END SCOPE

    // STEP 2: ensure data is persistent in the data-store after the
    // MeshCoordinates object goes out-of-scope
    EXPECT_TRUE( coords_group != AXOM_NULLPTR );

    // Ensure that the data remains persistent in the data-store
    sidre::Group* values_group = coords_group->getGroup( "values" );
    EXPECT_TRUE( values_group != AXOM_NULLPTR );
    for ( int j=0; j < dim; ++j )
    {
      sidre::View* coords_view =
          values_group->getView( std::string( coord_names[ j ] ) );

      double* coord_data =
          static_cast< double* >( coords_view->getVoidPtr() );

      detail::check_array_values( coord_data,
                                  expected_data[ j ],
                                  SMALL_NUM_NODES );
    }

  } // END for all dimensions

}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates, sidre_push_constructor )
{
  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double* data[ 3 ]           = { x, y, z };
  const char* coord_names[3]  = { "x", "y", "z" };

  for ( int dim=1; dim <= 3; ++dim )
  {
    sidre::DataStore ds;
    sidre::Group* coords_group = ds.getRoot();
    SLIC_ASSERT( coords_group != AXOM_NULLPTR );

    EXPECT_TRUE( coords_group->getNumGroups()==0 );
    EXPECT_TRUE( coords_group->getNumViews()==0 );

    // BEGIN SCOPE
    {
      mint::MeshCoordinates mesh_coords( coords_group, dim, SMALL_NUM_NODES );
      detail::check_coordinate_arrays( &mesh_coords );

      // ensure the sidre::Group is populated accordingly
      EXPECT_TRUE( coords_group->getNumGroups()==1 );
      EXPECT_TRUE( coords_group->getNumViews()==1 );
      EXPECT_TRUE( coords_group->hasView( "type" ) );
      EXPECT_TRUE( coords_group->hasChildGroup( "values" ) );

      sidre::Group* values_group = coords_group->getGroup( "values" );
      EXPECT_TRUE( values_group != AXOM_NULLPTR );
      EXPECT_TRUE( values_group->getNumViews()==static_cast< size_t> (dim ) );

      for ( int j=0; j < dim; ++j )
      {
        EXPECT_TRUE( values_group->hasChildView( coord_names[ j ] ) );

        sidre::View* coord_view = values_group->getView( coord_names[j] );
        EXPECT_TRUE( coord_view != AXOM_NULLPTR );

        EXPECT_EQ( mesh_coords.getCoordinateArray( j ),
                   coord_view->getVoidPtr() );
      }

      EXPECT_EQ( mesh_coords.dimension(), dim );
      EXPECT_EQ( mesh_coords.numNodes(), SMALL_NUM_NODES );
      EXPECT_TRUE( mesh_coords.numNodes() <= mesh_coords.capacity() );
      EXPECT_EQ( mesh_coords.capacity(), 2*SMALL_NUM_NODES );

      // populate the coordinates, writes to the corresponding sidre views
      mint::Array< double > xx( dim, 1, dim );
      for ( int inode=0; inode < SMALL_NUM_NODES; ++inode )
      {
        for ( int j=0; j < dim; ++j )
        {
          xx( j ) = data[ j ][ inode ];
        }

        mesh_coords.set( inode, xx.getData() );
      }

    }
    // END SCOPE

    // Ensure the data is persistent in sidre
    EXPECT_TRUE( coords_group->getNumGroups()==1 );
    EXPECT_TRUE( coords_group->getNumViews()==1 );
    EXPECT_TRUE( coords_group->hasView( "type" ) );
    EXPECT_TRUE( coords_group->hasChildGroup( "values" ) );

    sidre::Group* values_group = coords_group->getGroup( "values" );
    EXPECT_TRUE( values_group != AXOM_NULLPTR );
    EXPECT_TRUE( values_group->getNumViews()==static_cast< size_t >(dim) );

    for ( int j=0; j < dim; ++j )
    {
      EXPECT_TRUE( values_group->hasChildView( coord_names[ j ] ) );

      sidre::View* coord_view = values_group->getView( coord_names[j] );
      EXPECT_TRUE( coord_view != AXOM_NULLPTR );

      double* actual_data = static_cast< double* >( coord_view->getVoidPtr() );
      detail::check_array_values( actual_data, data[ j ], SMALL_NUM_NODES );
    }

  } // END for all dimensions

}

#endif /* MINT_USE_SIDRE */

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates_DeathTest, invalid_operations )
{
  // NOTE: this test ensures that when constructing a mint::MeshCoordinates
  // object with

  const char* IGNORE_OUTPUT   = ".*";
  constexpr mint::IndexType N = 4;
  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double y[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  double z[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };

  mint::MeshCoordinates coords( N, x, y, z );

  EXPECT_DEATH_IF_SUPPORTED( coords.append( x ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( coords.shrink(), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( coords.resize( 10 ), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( coords.reserve( 10 ), IGNORE_OUTPUT );
}

//------------------------------------------------------------------------------
TEST( mint_mesh_coordinates_DeathTest, invalid_construction )
{
  const char* IGNORE_OUTPUT = ".*";

  // STEP 0: test construction with invalid dimension
  EXPECT_DEATH_IF_SUPPORTED( mint::MeshCoordinates(0), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED( mint::MeshCoordinates(4), IGNORE_OUTPUT );
  EXPECT_DEATH_IF_SUPPORTED(
      mint::MeshCoordinates(0,SMALL_NUM_NODES,SMALL_NODE_CAPACITY),
      IGNORE_OUTPUT );

  // STEP 1: test construction with invalid numNodes,max_capacity settings
  EXPECT_DEATH_IF_SUPPORTED(
      mint::MeshCoordinates(2, LARGE_NUM_NODES, SMALL_NODE_CAPACITY),
      IGNORE_OUTPUT );

  // STEP 2: test invalid construction with null external buffers
  EXPECT_DEATH_IF_SUPPORTED( 
      mint::MeshCoordinates( 10, AXOM_NULLPTR ), IGNORE_OUTPUT );

  // STEP 3: test invalid construction from external buffers with 0 nodes
  double x[ SMALL_NUM_NODES ] = { 1.0, 2.0, 3.0, 4.0 };
  EXPECT_DEATH_IF_SUPPORTED(
      mint::MeshCoordinates(ZERO_NUM_NODES, x), IGNORE_OUTPUT );


#ifdef MINT_USE_SIDRE

  // STEP 4: test construction with a null Sidre group
  EXPECT_DEATH_IF_SUPPORTED(
      mint::MeshCoordinates( AXOM_NULLPTR ), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(
      mint::MeshCoordinates( AXOM_NULLPTR, 2, 10, 10 ),
      IGNORE_OUTPUT );

  // STEP 5: test sidre pull-constructor that does not conform to blueprint
  sidre::DataStore ds;
  EXPECT_DEATH_IF_SUPPORTED( mint::MeshCoordinates( ds.getRoot() ),
                             IGNORE_OUTPUT );

  // STEP 6: test sidre push-constructor with with a non-empty group
  detail::create_sidre_data( ds, 2 );
  EXPECT_DEATH_IF_SUPPORTED( mint::MeshCoordinates( ds.getRoot(),2,5,5),
                             IGNORE_OUTPUT );

#endif /* MINT_USE_SIDRE */

}

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
