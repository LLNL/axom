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

#include "gtest/gtest.h"

// Axom includes
#include "axom/Types.hpp" // for AXOM_NULLPTR

#include "axom_utils/Matrix.hpp"       // for Matrix class definition
#include "axom_utils/Determinants.hpp" // for numerics::determinant()

// Mint includes
#include "mint/CellType.hpp"
#include "mint/DataTypes.hpp"
#include "mint/FEBasis.hpp"
#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/FiniteElement.hpp"
#include "mint/Lagrange.hpp"
#include "mint/Mesh.hpp"
#include "mint/ShapeFunction.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"

// Slic includes
#include "slic/slic.hpp"

// C/C++ includes
#include <cmath>   // for sin(), cos(), sqrt(), etc.

//#define MINT_FEM_DEBUG

using namespace axom;

//------------------------------------------------------------------------------
//  INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

/*!
 * \brief Helper method to calculate the centroid of Finite Element instance.
 * \param [in] fe pointer to a finite element instance.
 * \param [out] centroid buffer to store the centroid
 */
void compute_centroid( mint::FiniteElement * fe, double * centroid )
{
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_TRUE(  centroid != AXOM_NULLPTR );

  const int ndims  = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  numerics::Matrix< double > physical_nodes(
    ndims, nnodes, fe->getPhysicalNodes(), true );

  if ( fe->getCellType() == MINT_PYRAMID )
  {

    // Hard-code expected centroid from VTK, since averaging the nodes
    // doesn't yield the centroid for a pyramid.
    //
    // NOTE: these values are hard-coded for the specific pyramid in this test.

    centroid[ 0 ] = 4.00001;
    centroid[ 1 ] = 4.28284;
    centroid[ 2 ] = 2.0;

  }
  else
  {

    for ( int i=0 ; i < ndims ; ++i )
    {

      typedef numerics::Matrix< double >::IndexType IndexType;
      IndexType p = 0;
      IndexType N = 0;
      const double * xp_i = physical_nodes.getRow(i,p,N);

      centroid[ i ] = 0.0;
      for ( int j=0 ; j < N ; j+=p )
      {
        centroid[ i ] += xp_i[ j ];
      } // END for all nodes

      centroid[ i ] /= nnodes;

    } // END for all dimensions

  } // END else

}

/*!
 * \brief Helper method to construct a single element (unstructured) mesh of
 *  the given cell type.
 *
 * \note Constructs the mesh element by scaling and rotating the reference
 *  element.
 *
 * \tparam BasisType the FEM basis of the element, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the cell type of the element, e.g., MINT_QUAD
 *
 * \post mesh != AXOM_NULLPTR
 * \post mesh->getNumberOfCells()==1
 * \post mesh->getNumberOfNodes()==N
 *
 * \note Ownership of the returned pointer is propagated to the caller. The
 *  calling method is therefore responsible for managing and properly
 *  deallocating the returned mesh object.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see get_fe_mesh()
 */
template < int BasisType, int CellType >
mint::UnstructuredMesh< CellType > * single_element_mesh( )
{
  typedef typename mint::UnstructuredMesh< CellType > MeshType;
  typedef typename mint::FEBasis< BasisType,CellType > FEMType;
  typedef typename FEMType::ShapeFunctionType ShapeFunctionType;

  const double SCALE = 10.0;
  const double ANGLE = 0.785398; // 45 degrees
  const double SINT  = sin( ANGLE );
  const double COST  = cos( ANGLE );

  const int ndims = ShapeFunctionType::dimension();
  EXPECT_TRUE( ndims==2 || ndims==3 );

  const int ndofs = ShapeFunctionType::numDofs();

  mint::localIndex * cell = new mint::localIndex[ ndofs ];

  double * center = new double[ ndims ];
  ShapeFunctionType::center( center );

  double * nodes  = new double[ ndofs*ndims ];
  ShapeFunctionType::coords( nodes );

  MeshType * m = new MeshType( ndims, ndofs, 1 );

  double centroid[]  = { 0.0, 0.0, 0.0}; // used to compute the pyramid apex

  for ( int i=0 ; i < ndofs ; ++i )
  {
    cell[ i ]    = i;
    double * node = &nodes[ i*ndims ];

    // scale & rotate cell from the reference space to get a  test cell
    const double dx   = node[ 0 ] - center[ 0 ];
    const double dy   = node[ 1 ] - center[ 1 ];

    node[ 0 ] = SCALE * ( (dx*COST - dy*SINT) + center[ 0 ] );
    node[ 1 ] = SCALE * ( (dx*SINT + dy*COST) + center[ 1 ] );
    if ( ndims==3 )
    {
      node[ 2 ] *= SCALE;
    }

    if ( CellType==MINT_PYRAMID && i < 4 )
    {
      centroid[ 0 ] += node[ 0 ];
      centroid[ 1 ] += node[ 1 ];
      centroid[ 2 ] += node[ 2 ];
    }

    if ( CellType==MINT_PYRAMID && i==4 )
    {
      // generate right pyramid, ensure the apex is prependicular to the
      // base of the pyramid to facilitate testing.
      node[ 0 ] = 0.25*centroid[ 0 ];
      node[ 1 ] = 0.25*centroid[ 1 ];
      node[ 2 ] = 0.25*centroid[ 2 ] + SCALE;
    }

    m->addNode( node );
  }

  m->addCell( cell, CellType );

#ifdef MINT_FEM_DEBUG
  std::string vtkFile = std::string( mint::cell::name[ CellType ] ) + ".vtk";
  mint::write_vtk( m, vtkFile );
#endif

  // clean up
  delete [] center;
  delete [] nodes;
  delete [] cell;

  return ( m );
}

/*!
 * \brief Constructs a finite element mesh consisting of a single element.
 *
 * \param [out] m pointer to the mesh object
 * \param [out] fe finite element object bound to the mesh
 *
 * \note Ownership of the mesh and finite element objects is propagated to the
 *  caller. The calling method is therefore responsible for managing and
 *  properly deallocating these objects.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see single_element_mesh()
 */
template < int BasisType, int CellType >
void get_fe_mesh( mint::UnstructuredMesh< CellType > *& m,
                  mint::FiniteElement *& fe )
{
  EXPECT_TRUE(  m==AXOM_NULLPTR );
  EXPECT_TRUE(  fe==AXOM_NULLPTR );

  // STEP 0: construct a mesh with a single element
  m = single_element_mesh< BasisType, CellType >( );

  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 1: construct FE instance.
  fe = new mint::FiniteElement( m, 0 );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_TRUE(  fe->getBasisType()==MINT_UNDEFINED_BASIS );

  mint::bind_basis< BasisType, CellType >( *fe );
  EXPECT_FALSE( fe->getBasisType()==MINT_UNDEFINED_BASIS );
}

/*!
 * \brief Ensures the FiniteElement instance is properly constructed by
 *  checking the underlying ShapeFunction.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != AXOM_NULLPTR
 */
template < typename ShapeFunctionType >
void check_reference_element( mint::FiniteElement * fe )
{
  EXPECT_TRUE( fe != AXOM_NULLPTR);
  EXPECT_FALSE( fe->getBasisType()==MINT_UNDEFINED_BASIS );

  // STEP 0: check reference element attributes, eg, dofs, dimension, etc.
  ShapeFunctionType sf;
  EXPECT_EQ(  sf.numDofs(),   fe->getNumDofs() );
  EXPECT_EQ(  sf.dimension(), fe->getReferenceDimension() );

  const int fe_dim = fe->getPhysicalDimension();
  EXPECT_TRUE(  (fe_dim >= 1) && (fe_dim <= 3) );

  const int ref_dim = fe->getReferenceDimension();
  EXPECT_TRUE(  ref_dim <= fe_dim );

  // STEP 1: check reference coordinates
  const int ndofs = sf.numDofs();
  const int ndims = sf.dimension();

  // populate a matrix with shape function coordinates (expected)
  numerics::Matrix< double > shape_coords( ndims, ndofs );
  sf.coords( shape_coords.data() );

  // populate a matrix with the element's reference coordinates to test
  numerics::Matrix< double > ref_coords(
    ndims, ndofs, fe->getReferenceNodes(), true );

  for ( int i=0 ; i < ndofs ; ++i )
  {

    const double * xr_expected = shape_coords.getColumn( i );
    const double * xr          = ref_coords.getColumn( i );

    for ( int j=0 ; j < ndims ; ++j )
    {
      EXPECT_DOUBLE_EQ( xr_expected[ j ], xr[ j ] );
    }

  } // END for all dofs

  // STEP 2: check reference center
  double * shape_center = new double[ ndims ];
  sf.center( shape_center );

  for ( int j=0 ; j < ref_dim ; ++j )
  {
    const double expected = shape_center[ j ];
    const double actual   = fe->getReferenceCenter( )[ j ];
    EXPECT_DOUBLE_EQ( expected, actual );
  }

  // STEP 3: clean up
  delete [] shape_center;
}

/*!
 * \brief Checks the mapping from reference space to physical space.
 *
 *  The test loops over all the nodes of the reference element, maps them to
 *  physical space and makes sure they match with the corresponding physical
 *  coordinates of the element.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != AXOM_NULLPTR
 *
 * \see check_forward_map
 */
void test_forward_map( mint::FiniteElement * fe, double TOL=1.e-9 )
{
  EXPECT_TRUE( fe != AXOM_NULLPTR );
  EXPECT_FALSE( fe->getBasisType()==MINT_UNDEFINED_BASIS );
  EXPECT_TRUE(  fe->getReferenceDimension()==fe->getPhysicalDimension() );
  EXPECT_TRUE(  fe->getNumNodes()==fe->getNumDofs() );

  int numdofs = fe->getNumDofs();
  int nnodes  = fe->getNumNodes();
  int ndims   = fe->getReferenceDimension();

  // STEP 0: get the reference coordinates
  numerics::Matrix< double > reference_coords(
    ndims, numdofs, fe->getReferenceNodes(), true );

  // STEP 1: get the physical coordinates of the element
  numerics::Matrix< double > element_coords(
    ndims, nnodes, fe->getPhysicalNodes(), true );

  // STEP 2: buffer to store the mapped physical coordinates
  numerics::Matrix< double > physical_coords( ndims, nnodes );

  // STEP 3: map each reference coordinate to physical
  for ( int i=0 ; i < numdofs ; ++i )
  {

    // forward map
    const double * xr = reference_coords.getColumn( i );
    double * xp       = physical_coords.getColumn( i );
    fe->computePhysicalCoords( xr, xp );

    double * expected_xp = element_coords.getColumn( i );

    // check mapping
    for ( int j=0 ; j < ndims ; ++j )
    {
      EXPECT_NEAR( expected_xp[ j ], xp[ j ], TOL );
    }

  } // END for all dofs

  // STEP 3: map center
  const double * xr_center = fe->getReferenceCenter();
  double * xp_center       = new double[ ndims ];
  fe->computePhysicalCoords( xr_center, xp_center );

  // STEP 4: compute expected centroid by averaging coordinates
  double * expected_centroid = new double[ ndims ];
  compute_centroid( fe, expected_centroid );

  // STEP 5: check mapping of centroid
  for ( int i=0 ; i < ndims ; ++i )
  {
    EXPECT_NEAR( expected_centroid[ i ], xp_center[ i ], TOL );
  }

  // STEP 6: cleanup
  delete [] xp_center;
  delete [] expected_centroid;
}

/*!
 * \brief Checks the mapping from physical space to reference space.
 *
 *  The test loops over all the nodes of the element in physical space, maps
 *  them to the reference space and makes sure they match with the corresponding
 *  reference coordinates of the element.
 *
 * \param [in] fe pointer to the finite element instance to test
 *
 * \pre fe != AXOM_NULLPTR
 *
 * \see check_inverse_map
 */
void test_inverse_map( mint::FiniteElement * fe, double TOL=1.e-9 )
{
  EXPECT_TRUE( fe != AXOM_NULLPTR );
  EXPECT_FALSE( fe->getBasisType()==MINT_UNDEFINED_BASIS );
  EXPECT_TRUE(  fe->getReferenceDimension()==fe->getPhysicalDimension() );
  EXPECT_TRUE(  fe->getNumNodes()==fe->getNumDofs() );

  // STEP 0: get Matrix of physical nodes; nodal coordinates stored in columns.
  const int ndims  = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  numerics::Matrix< double > physical_nodes(
    ndims, nnodes, fe->getPhysicalNodes(), true );

  // STEP 1: get reference coordinates
  numerics::Matrix< double > reference_nodes(
    ndims, nnodes, fe->getReferenceNodes(), true );

  // STEP 2: allocate buffer to store the computed reference coordinates
  double * xr = new double[ ndims ];

  // STEP 3: loop over physical nodes and map them to reference space
  for ( int i=0 ; i < nnodes ; ++i )
  {

    if ( fe->getCellType()==MINT_PYRAMID && i==4 )
    {
      // skip inverse map at the apex of the pyramid, system is singular
      continue;
    }

    const double * expected_xr = reference_nodes.getColumn( i );
    const double * xp          = physical_nodes.getColumn( i );
    int rc = fe->computeReferenceCoords( xp, xr, TOL );
    EXPECT_TRUE( rc==mint::INSIDE_ELEMENT );

    // check mapping
    for ( int j=0 ; j < ndims ; ++j )
    {
      EXPECT_NEAR( expected_xr[ j ], xr[ j ], TOL );
    }

  } // END for all physical nodes

  // STEP 4: calculate centroid in physical space by averaging coordinates
  double * centroid = new double[ ndims ];
  compute_centroid( fe, centroid );

  // STEP 4: map centroid
  int rc = fe->computeReferenceCoords( centroid, xr, TOL );
  EXPECT_TRUE( rc==mint::INSIDE_ELEMENT );

  // STEP 5: check centroid mapping
  const double * refcenter = fe->getReferenceCenter();
  for ( int i=0 ; i < ndims ; ++i )
  {
    EXPECT_NEAR( refcenter[ i ], xr[ i ], TOL );
  }

  // STEP 6: clean up
  delete [] xr;
  delete [] centroid;
}

/*!
 * \brief Performs basic checks on a FiniteElement instance bound to a
 *  corresponding basis function.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see check_reference_element()
 */
template < int BasisType, int CellType >
void check_shape( )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  typedef typename mint::FEBasis< BasisType, CellType > FEMType;
  typedef typename FEMType::ShapeFunctionType ShapeFunctionType;
  typedef typename mint::UnstructuredMesh< CellType > MeshType;

  // STEP 0: construct finite element mesh
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE(  m != AXOM_NULLPTR );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 1: test FE instance
  check_reference_element< ShapeFunctionType >( fe );

  // STEP 2: clean up
  delete fe;
  delete m;
}

/*!
 * \brief Ensures the jacobian is positive within the element by evaluating
 *  the determinant of the jacobian at the element nodes, center and computed
 *  interior points.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 */
template < int BasisType, int CellType >
void check_jacobian( double TOL=1.e-9 )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  const double LTOL = 0.0-TOL;
  double det = 0.0;

  typedef typename mint::UnstructuredMesh< CellType > MeshType;

  // STEP 0: construct a mesh with a single element
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE(  m != AXOM_NULLPTR );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 1: construct a Matrix object to store the jacobian
  const int ndims = fe->getPhysicalDimension();
  numerics::Matrix< double > J( ndims, ndims );

  // STEP 2: test jacobian at the reference nodes
  const int ndofs = fe->getNumDofs();
  const int rdim  = fe->getReferenceDimension();
  numerics::Matrix< double > n(rdim,ndofs, fe->getReferenceNodes(), true );
  for ( int i=0 ; i < ndofs ; ++i )
  {
    const double * xi = n.getColumn( i );
    fe->jacobian( xi, J );

    det = numerics::determinant( J );
    EXPECT_GT( det, LTOL );
  }

  // STEP 3: test jacobian at the reference center
  const double * xi_c = fe->getReferenceCenter();
  fe->jacobian( xi_c, J );
  det = numerics::determinant( J );
  EXPECT_GT( det, LTOL );

  // STEP 4: test jacobian at interior points; the interior points are
  // computed by taking the midpoint of a reference node and the centroid
  double * rp = new double[ ndims ];
  for ( int i=0 ; i < ndofs ; ++i )
  {

    const double * xi = n.getColumn( i );
    for ( int j=0 ; j < ndims ; ++j )
    {
      rp[ j ] = 0.5 * ( xi[ j ] +  xi_c[ j ] );
    }

    fe->jacobian( rp, J );
    det = numerics::determinant( J );
    EXPECT_GT( det, LTOL );
  }

  // STEP 5: clean up
  delete []  rp;
  delete fe;
  delete m;
}

/*!
 * \brief Checks the forward map of a FiniteElement object.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see test_forward_map()
 */
template < int BasisType, int CellType >
void check_forward_map( double TOL=1.e-9 )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  typedef typename mint::UnstructuredMesh< CellType > MeshType;

  // STEP 0: construct a mesh with a single element
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE(  m != AXOM_NULLPTR );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 1: check forward mapping
  test_forward_map( fe, TOL );

  // STEP 2: clean up
  delete fe;
  delete m;
}

/*!
 * \brief Checks the inverse map of a FiniteElement object.
 *
 * \param [in] TOL optional user-supplied tolerange. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 *
 * \see test_inverse_map()
 */
template < int BasisType, int CellType >
void check_inverse_map( double TOL=1.e-9 )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  typedef typename mint::UnstructuredMesh< CellType > MeshType;

  // STEP 0: construct a mesh with a single element
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE(  m != AXOM_NULLPTR );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 1: check inverse map
  test_inverse_map( fe, TOL );

  // STEP 2: clean up
  delete fe;
  delete m;
}

/*!
 * \brief Checks correctness of the inverse map for point-in-cell queries.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 */
template < int BasisType, int CellType >
void point_in_cell( double TOL=1.e-9 )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  typedef typename mint::UnstructuredMesh< CellType > MeshType;

  // STEP 0: construct a mesh with a single element
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE(  m != AXOM_NULLPTR );
  EXPECT_TRUE(  fe != AXOM_NULLPTR );
  EXPECT_EQ(  mint::cell::num_nodes[ CellType ],  m->getNumberOfNodes() );
  EXPECT_EQ(  1,                                  m->getNumberOfCells() );

  // STEP 0: test variables
  const int nnodes = fe->getNumNodes();
  const int ndims  = fe->getPhysicalDimension();
  double * xi       = new double[ ndims ];
  double * xc       = new double[ ndims ];
  double * rp       = new double[ ndims ];
  int status       = 0;

  numerics::Matrix< double > nodes( ndims,nnodes,fe->getPhysicalNodes(),true );

  // STEP 1: Ensure element nodes are inside
  for ( int i=0 ; i < nnodes ; ++i )
  {

    if ( fe->getCellType()==MINT_PYRAMID && i==4 )
    {
      // skip inverse map at the apex of the pyramid, system is singular
      continue;
    }

    const double * xp = nodes.getColumn( i );

    status = fe->computeReferenceCoords( xp, xi, TOL );
    EXPECT_TRUE( status==mint::INSIDE_ELEMENT );
  }

  // STEP 2: Ensure center is inside
  fe->computePhysicalCoords( fe->getReferenceCenter(), xc );

  status = fe->computeReferenceCoords( xc, xi, TOL );
  EXPECT_TRUE( status==mint::INSIDE_ELEMENT );

  // STEP 3: Ensure other interior nodes are inside
  numerics::Matrix< double > dofs( ndims, nnodes, fe->getReferenceNodes(),
                                   true );

  for ( int i=0 ; i < nnodes ; ++i )
  {

    const double * dof = dofs.getColumn( i );

    for ( int j=0 ; j < ndims ; ++j )
    {
      rp[ j ] = 0.5 * ( dof[ j ] +  fe->getReferenceCenter()[ j ] );
    }

    fe->computePhysicalCoords( rp, xc );

    status = fe->computeReferenceCoords( xc, xi, TOL );
    EXPECT_TRUE( status==mint::INSIDE_ELEMENT );

  } // END for

  // STEP 4: test outside points by shifting reference element
  const double SHIFT = 10;
  for ( int i=0 ; i < nnodes ; ++i )
  {

    const double * dof = dofs.getColumn( i );

    // Test +SHIFT
    for ( int j=0 ; j < ndims ; ++j )
    {
      rp[ j ] = dof[ j ] + SHIFT;
    }

    fe->computePhysicalCoords( rp, xc );

    status = fe->computeReferenceCoords( xc, xi, TOL );
    EXPECT_TRUE( status==mint::OUTSIDE_ELEMENT );

    // Test -SHIFT
    for ( int j=0 ; j < ndims ; ++j )
    {
      rp[ j ] = dof[ j ] - SHIFT;
    }

    fe->computePhysicalCoords( rp, xc );

    status = fe->computeReferenceCoords( xc, xi, TOL );
    EXPECT_TRUE( status==mint::OUTSIDE_ELEMENT );
  }

  // STEP 5: clean up
  delete [] xi;
  delete [] xc;
  delete [] rp;
  delete fe;
  delete m;
}

/*!
 * \brief Computes a fake node-centered field, using the nodal coordinates, to
 *  test interpolation within the element.
 *
 * \param [in] x pointer to the coordinates of the node.
 * \param [in] N the number of dimensions
 *
 * \return f the computed field.
 *
 * \see check_interp()
 */
double analytic_function( const double * x, int N )
{
  double f = 0.0;
  for ( int i=0 ; i < N ; ++i )
  {
    f += x[ i ];
  }
  return ( f );
}

/*!
 * \brief Given field values and associated weights at the nodes of an element,
 *  this method interpolates the field.
 *
 * \param [in] f user-supplied buffer where the nodal field values are stored
 * \param [in] wgts user-supplied buffer of corresponding interpolation weigths
 * \param [in] N the number of nodes, i.e., the length of f and wgts
 *
 * \pre f != AXOM_NULLPTR
 * \pre wgts != AXOM_NULLPTR
 *
 * \return finterp the interpolated field value
 *
 * \see check_interp()
 */
double interp( const double * f, const double * wgts, int N )
{
  EXPECT_TRUE(  f != AXOM_NULLPTR );
  EXPECT_TRUE(  wgts != AXOM_NULLPTR );

  double finterp = 0.0;
  for ( int i=0 ; i < N ; ++i )
  {
    finterp += wgts[ i ]*f[ i ];
  }
  return ( finterp );
}

/*!
 * \brief Checks the interpolation within the element against a field that
 *  is defined analytically.
 *
 * \param [in] TOL optional user-supplied tolerance. Default is 1.e-9.
 *
 * \tparam BasisType basis bound to the FiniteElemen, e.g., MINT_LAGRANGE_BASIS
 * \tparam CellType the corresponding cell type, e.g., MINT_QUAD
 */
template < int BasisType, int CellType >
void check_interp( double TOL=1.e-9 )
{
  EXPECT_TRUE(  (CellType >= 0) && (CellType < MINT_NUM_CELL_TYPES) );
  EXPECT_TRUE(  (BasisType >= 0) && (BasisType < MINT_NUM_BASIS_TYPES) );

  SLIC_INFO( "checking " << mint::basis_name[ BasisType ] << " / "
                         << mint::cell::name[ CellType ] );

  // STEP 0: construct a mesh with a single element
  mint::UnstructuredMesh< CellType > * m = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;

  get_fe_mesh< BasisType, CellType >( m, fe );
  EXPECT_TRUE( m != AXOM_NULLPTR );
  EXPECT_TRUE( fe != AXOM_NULLPTR );
  EXPECT_EQ( mint::cell::num_nodes[ CellType ], m->getNumberOfNodes() );
  EXPECT_EQ( 1, m->getNumberOfCells() );

  const int ndims  = fe->getPhysicalDimension();
  const int nnodes = fe->getNumNodes();
  double * wgts = new double[ nnodes ];

  // STEP 1: setup a nodal field to interpolate
  // Note: I don't know why we need to cast as m as a Mesh to add a field
  // it works just fine other places. -BC
  m->getNodeFieldData().addField(
       new mint::FieldVariable< double >("foo", nnodes ) );
  double* f = m->getNodeFieldData().getField( "foo" )->getDoublePtr();

//  mint::Field * F = mesh->addNodeField< double >( "foo", 1 );
//  double * f = F->getDoublePtr();

  numerics::Matrix< double > nodes( ndims,nnodes,fe->getPhysicalNodes(),true );

  for ( int i=0 ; i < nnodes ; ++i )
  {
    const double * x = nodes.getColumn( i );
    f[ i ] = analytic_function( x, ndims );
  }

  // STEP 2: check interpolation at the nodes
  const int ndofs = fe->getNumDofs();
  numerics::Matrix< double > dofs( ndims,ndofs,fe->getReferenceNodes(),true );
  for ( int i=0 ; i < ndofs ; ++i )
  {
    const double * xi = dofs.getColumn( i );
    fe->evaluateShapeFunctions( xi, wgts );
    double finterp = interp( f, wgts, ndofs );
    EXPECT_NEAR( f[i], finterp, TOL );
  }

  // STEP 3: compute centroid
  double * xc = new double[ ndims ];
  fe->computePhysicalCoords( fe->getReferenceCenter(), xc );

  // STEP 4: interpolate
  fe->evaluateShapeFunctions( fe->getReferenceCenter(), wgts );
  double finterp = interp( f, wgts, nnodes );

  // STEP 5: check interpolation
  const double fexpected = analytic_function(xc,ndims);
  EXPECT_NEAR( fexpected, finterp, TOL );

  // STEP 6: clean up
  delete m;
  delete fe;
  delete [] xc;
  delete [] wgts;
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( mint_single_fe, check_override_max_newton )
{
  typedef mint::UnstructuredMesh< MINT_QUAD > MeshType;

  const int MAX_NEWTON = 42; // test value to override max newton

  // STEP 0: construct a mesh with a single element
  MeshType * m             = AXOM_NULLPTR;
  mint::FiniteElement * fe = AXOM_NULLPTR;
  get_fe_mesh< MINT_LAGRANGE_BASIS, MINT_QUAD >( m, fe );

  EXPECT_FALSE( MAX_NEWTON==fe->getMaxSolverIterations() );

  // STEP 1: override max newton iterations
  fe->setMaxSolverIterations( MAX_NEWTON );
  EXPECT_EQ( MAX_NEWTON, fe->getMaxSolverIterations() );

  // clean up
  delete fe;
  delete m;
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, matrix_constructor_deepcopy )
{
  // STEP 0: constants used in the test
  const int NROWS = 2;
  const int NCOLS = 3;
  const int NSIZE = NROWS * NCOLS;

  // STEP 1: setup element coordinates matrix
  numerics::Matrix< double > M( NROWS, NCOLS );
  M.fillColumn( 0,  1.0 );
  M.fillColumn( 1,  2.0 );
  M.fillColumn( 2,  3.0 );

  // STEP 2: construct FE object by making a deep-copy
  mint::FiniteElement fe( M, MINT_TRIANGLE );
  EXPECT_FALSE( fe.usesExternalBuffer() );

  // STEP 3: ensure FE and Matrix objects are pointing to the same buffer
  double * physical_nodes = fe.getPhysicalNodes();
  EXPECT_FALSE( physical_nodes == M.data() );

  // STEP 4: ensure the contents are the same
  const double * expected = M.data();
  for ( int i=0 ; i < NSIZE ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected[ i ], physical_nodes[ i ] );
  }

  // STEP 5: change the matrix data
  M.swapColumns(  0,  2 );
  M.swapColumns(  1,  2 );

  // STEP 6: ensure the contents *are not* the same
  for ( int i=0 ; i < NSIZE ; ++i )
  {
    EXPECT_FALSE( utilities::isNearlyEqual( M.data()[i], physical_nodes[i] ) );
  }

}

//------------------------------------------------------------------------------
TEST( mint_single_fe, matrix_constructor_shallowcopy)
{
  // STEP 0: constants used in the test
  const int NROWS = 2;
  const int NCOLS = 3;
  const int NSIZE = NROWS * NCOLS;

  // STEP 1: setup element coordinates matrix
  numerics::Matrix< double > M( 2, 3 );
  M.fillColumn( 0,  1.0 );
  M.fillColumn( 1,  2.0 );
  M.fillColumn( 2,  3.0 );

  // STEP 2: construct FE object by making a shallow-copy
  mint::FiniteElement * fe = new mint::FiniteElement( M, MINT_TRIANGLE, true );
  EXPECT_TRUE( fe->usesExternalBuffer() );

  // STEP 3: ensure FE and Matrix objects are pointing to the same buffer
  double * physical_nodes = fe->getPhysicalNodes();
  EXPECT_TRUE( physical_nodes == M.data() );

  // STEP 4: ensure the contents are the same
  const double * expected = M.data();
  for ( int i=0 ; i < NSIZE ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected[ i ], physical_nodes[ i ] );
  }

  // STEP 5: change the matrix data
  M.swapColumns(  0,  2 );
  M.swapColumns(  1,  2 );

  // STEP 6: ensure the contents *still* the same
  for ( int i=0 ; i < NSIZE ; ++i )
  {
    EXPECT_DOUBLE_EQ( expected[ i ], physical_nodes[ i ] );
  }

  // STEP 7: delete the FE object, ensure Matrix buffer is not corrupted
  delete fe;
  fe = AXOM_NULLPTR;
  EXPECT_FALSE( M.data()==AXOM_NULLPTR );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_shape_function )
{
  check_shape< MINT_LAGRANGE_BASIS, MINT_QUAD >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_TET >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_HEX >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_PRISM >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_PYRAMID >( );

  check_shape< MINT_LAGRANGE_BASIS, MINT_QUAD9 >( );
  check_shape< MINT_LAGRANGE_BASIS, MINT_HEX27 > ( );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_jacobian )
{
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_QUAD >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_TET >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_HEX >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_PRISM >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_PYRAMID >( );

  check_jacobian< MINT_LAGRANGE_BASIS, MINT_QUAD9 >( );
  check_jacobian< MINT_LAGRANGE_BASIS, MINT_HEX27 >( );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_forward_map )
{
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_QUAD >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_TET >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_HEX >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_PRISM >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_PYRAMID >( 1.e-5 );

  check_forward_map< MINT_LAGRANGE_BASIS, MINT_QUAD9 >( );
  check_forward_map< MINT_LAGRANGE_BASIS, MINT_HEX27 >( );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_inverse_map )
{
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_QUAD >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_TET >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_HEX >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_PRISM >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_PYRAMID >( 1.e-5 );

  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_QUAD9 >( );
  check_inverse_map< MINT_LAGRANGE_BASIS, MINT_HEX27 >( );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_point_in_cell )
{
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_QUAD >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_TET >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_HEX >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_PRISM >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_PYRAMID >( );

  point_in_cell< MINT_LAGRANGE_BASIS, MINT_QUAD9 >( );
  point_in_cell< MINT_LAGRANGE_BASIS, MINT_HEX27 >( 1.e-7 );
}

//------------------------------------------------------------------------------
TEST( mint_single_fe, check_fe_interp )
{
  check_interp< MINT_LAGRANGE_BASIS, MINT_QUAD >(     1.e-24 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_TRIANGLE >( 1.e-12 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_TET >(      1.e-12 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_HEX >(      1.e-24 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_PRISM >(    1.e-12 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_PYRAMID >(  1.e-12 );

  check_interp< MINT_LAGRANGE_BASIS, MINT_QUAD9 >(    1.e-24 );
  check_interp< MINT_LAGRANGE_BASIS, MINT_HEX27 >(    1.e-24 );
}

//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
