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

#include "FiniteElement.hpp"

// Axom includes
#include "axom/Macros.hpp" // for AXOM macros
#include "axom/Types.hpp"  // for AXOM_NULLPTR

// Axom Utils includes
#include "axom_utils/Utilities.hpp"    // for abs()
#include "axom_utils/linear_solve.hpp" // for linear_solve()

// Mint includes
#include "mint/CellType.hpp"           // for cell type definitions
#include "mint/Lagrange.hpp"           // For Lagrange ShapeFunctions
#include "mint/Mesh.hpp"               // for mesh data structure
#include "mint/ShapeFunction.hpp"      // For ShapeFunction definition

// Slic includes
#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace detail
{

//------------------------------------------------------------------------------
bool diverged( const double* xi, int N )
{
  const double DIVERGED = 1.e6;

  bool divergence_detected = false;
  for ( int i=0 ; !divergence_detected && (i < N) ; ++i )
  {
    divergence_detected = divergence_detected || (xi[ i ] > DIVERGED);
  }

  return divergence_detected;
}

}

//------------------------------------------------------------------------------
// FINITE ELEMENT CLASS IMPLEMENTATION
//------------------------------------------------------------------------------
FiniteElement::FiniteElement( const Mesh* mesh, int cellIdx ) :
  m_dim( mesh->getDimension() ),
  m_ctype( mesh->getMeshCellType( cellIdx ) ),
  m_shape_func_type( MINT_UNDEFINED_BASIS ),
  m_maxNewtonIterations( -1 ),
  m_numnodes( mesh->getMeshNumberOfCellNodes( cellIdx ) ),
  m_jac( AXOM_NULLPTR ),
  m_xyz( AXOM_NULLPTR ),
  m_phi( AXOM_NULLPTR ),
  m_phidot( AXOM_NULLPTR ),
  m_usingExternal( false ),
  m_shapeFunction( AXOM_NULLPTR ),
  m_shapeFunctionDerivatives( AXOM_NULLPTR ),
  m_reference_min( -1 ),
  m_reference_max( -1 ),
  m_reference_dim( -1 ),
  m_numdofs( -1 ),
  m_reference_coords( AXOM_NULLPTR ),
  m_reference_center( AXOM_NULLPTR )
{
  this->setUp( );
  this->getCellCoords( mesh, cellIdx );
}

//------------------------------------------------------------------------------
FiniteElement::FiniteElement( numerics::Matrix< double >& M,
                              int cellType,
                              bool useExternal ) :
  m_dim( M.getNumRows() ),
  m_ctype( cellType ),
  m_shape_func_type( MINT_UNDEFINED_BASIS ),
  m_maxNewtonIterations( -1 ),
  m_numnodes( M.getNumColumns() ),
  m_jac( AXOM_NULLPTR ),
  m_xyz( AXOM_NULLPTR ),
  m_phi( AXOM_NULLPTR ),
  m_phidot( AXOM_NULLPTR ),
  m_usingExternal( useExternal ),
  m_shapeFunction( AXOM_NULLPTR ),
  m_shapeFunctionDerivatives( AXOM_NULLPTR ),
  m_reference_min( -1 ),
  m_reference_max( -1 ),
  m_reference_dim( -1 ),
  m_numdofs( -1 ),
  m_reference_coords( AXOM_NULLPTR ),
  m_reference_center( AXOM_NULLPTR )
{
  this->setUp( );

  if ( m_usingExternal )
  {

    // shallow copy the data
    m_xyz = M.data( );

  }
  else
  {

    // make a deep copy of the data
    const int N = m_dim*m_numnodes;
    for ( int i=0 ; i < N ; ++i )
    {
      m_xyz[ i ] = M.data()[ i ];
    }

  }

}

//------------------------------------------------------------------------------
FiniteElement::~FiniteElement()
{
  this->tearDown( );
}

//------------------------------------------------------------------------------
int FiniteElement::computeReferenceCoords( const double* xp,
                                           double* xr,
                                           double TOL )
{
  SLIC_ASSERT(  xp != AXOM_NULLPTR );
  SLIC_ASSERT(  xr != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getBasisType() != MINT_UNDEFINED_BASIS );

  if ( this->getBasisType() == MINT_UNDEFINED_BASIS )
  {
    SLIC_WARNING( "No associated FiniteElement basis!" );
    return INVERSE_MAP_FAILED;
  }

  const int MAX_DIM = 3;
  double l1norm     = 0.0;
  bool converged    = false;

  double x[ MAX_DIM ];     // residual vector
  double psi[ MAX_DIM ];   // rhs

  // STEP 0: matrix instance to store the jacobian, J
  numerics::Matrix< double > J( m_dim, m_dim, m_jac, true );
  numerics::Matrix< double > coord_matrix( m_dim, m_numnodes, m_xyz,  true );

  // STEP 1: set initial guess for Newton-Raphson at the parametric center
  const int ref_dim = this->getReferenceDimension();
  for ( int i=0 ; i < m_dim ; ++i )
  {
    xr[ i ] = ( i < ref_dim ) ? this->getReferenceCenter()[ i ] : 0.0;
  }

  // STEP 2: start Newton-Raphson iteration
  for ( int iter=0 ; !converged && (iter < m_maxNewtonIterations) ; ++iter )
  {

    // evaluate the shape functions and jacobian @ \xi
    this->evaluateShapeFunctions( xr, m_phi );
    this->jacobian( xr, J );

    // compute the right-hand side term, -\psi
    for ( int i=0 ; i < m_dim ; ++i )
    {

      psi[ i ] = (-1.0)*xp[ i ];
      for ( int j=0 ; j < m_numnodes ; ++j )
      {
        psi[ i ] += m_phi[ j ] * coord_matrix(i,j);
      }

      psi[ i ] = (-1.0) * psi[ i ];
    }

    // calculate residual
    int rc = numerics::linear_solve( J, psi, x );

    // compute l1norm and generate improvements
    l1norm = 0.0;
    for ( int i=0 ; i < m_dim ; ++i )
    {
      l1norm  += utilities::abs( x[ i ] );
      xr[ i ] += x[ i ];
    }

    if ( rc != 0 )
    {
      double det = numerics::determinant( J );
      SLIC_WARNING( "Newton-Raphson failed, system appears singular!" );
      SLIC_INFO(  "singular system => det(J)=" << det );
      SLIC_INFO(  "l1norm=" << l1norm );
      SLIC_INFO(  "Newton iteration=" << iter );

      double xr1 = xr[ 0 ];
      double xr2 = xr[ 1 ];
      double xr3 = ( m_dim==3 ) ? xr[ 2 ] : 0.0;
      SLIC_INFO( "xr=[" << xr1 << " " << xr2 << " " << xr3 << "]" );
      break;
    }

    // convergence check
    converged = ( l1norm < TOL ) ? true : false;

    // divergence check
    if ( !converged && detail::diverged( xr, m_dim ) )
    {
      SLIC_WARNING( "Newton-Raphson divergence detected!" );
      SLIC_INFO( "l1norm=" << l1norm << " iter= " << iter );
      break;
    }

  } // END newton iterations

  if ( !converged )
  {

    SLIC_WARNING( "Newton-Raphson did not converge!" );
    SLIC_INFO( "l1norm=" << l1norm );
    return INVERSE_MAP_FAILED;
  }

  int rc =
    ( this->inReferenceElement(xr,TOL) ) ? INSIDE_ELEMENT : OUTSIDE_ELEMENT;

  return rc;
}

//------------------------------------------------------------------------------
void FiniteElement::computePhysicalCoords( const double* xr, double* xp )
{
  SLIC_ASSERT(  xr != AXOM_NULLPTR );
  SLIC_ASSERT(  xp != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getBasisType() != MINT_UNDEFINED_BASIS );

  if ( this->getBasisType() == MINT_UNDEFINED_BASIS )
  {
    SLIC_WARNING( "No associated FiniteElement basis!" );
    return;
  }

  // STEP 0: evaluate the shape functions, update m_phi
  this->evaluateShapeFunctions( xr, m_phi );

  // STEP 1: get a coordinates matrix instance
  numerics::Matrix< double > coords_matrix( m_dim, m_numnodes, m_xyz, true );

  // STEP 2: compute the global (physical) coordinates
  numerics::vector_multiply( coords_matrix, m_phi, xp );
}

//------------------------------------------------------------------------------
void FiniteElement::jacobian( const double* lc,
                              numerics::Matrix< double >& J  )
{
  SLIC_ASSERT(  lc != AXOM_NULLPTR );
  SLIC_ASSERT(  J.getNumColumns() == this->getReferenceDimension() );
  SLIC_ASSERT(  J.getNumRows() == this->getPhysicalDimension() );
  SLIC_ASSERT(  this->getBasisType() != MINT_UNDEFINED_BASIS );

  if ( this->getBasisType() == MINT_UNDEFINED_BASIS )
  {
    SLIC_WARNING( "No associated FiniteElement basis!" );
    return;
  }

  // STEP 0: evaluate the derivatives
  this->evaluateDerivatives( lc, m_phidot );
  numerics::Matrix< double > derivs_matrix(
    m_numdofs, this->getReferenceDimension(), m_phidot, true );

  // STEP 1: get the coordinates
  numerics::Matrix< double > coords_matrix( m_dim, m_numnodes, m_xyz, true );

  // STEP 2: compute the jacobian
  numerics::matrix_multiply( coords_matrix, derivs_matrix, J );
}

//------------------------------------------------------------------------------
void FiniteElement::evaluateShapeFunctions( const double* xr, double* phi )
{
  SLIC_ASSERT(  xr != AXOM_NULLPTR );
  SLIC_ASSERT(  phi != AXOM_NULLPTR );
  SLIC_ASSERT(  m_shapeFunction != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getBasisType() != MINT_UNDEFINED_BASIS );

  if ( this->getBasisType() == MINT_UNDEFINED_BASIS )
  {
    SLIC_WARNING( "No associated FiniteElement basis!" );
    return;
  }

  m_shapeFunction( xr, phi );
}

//------------------------------------------------------------------------------
void FiniteElement::evaluateDerivatives( const double* xr, double* phidot )
{
  SLIC_ASSERT(  xr != AXOM_NULLPTR );
  SLIC_ASSERT(  phidot != AXOM_NULLPTR );
  SLIC_ASSERT(  m_shapeFunctionDerivatives != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getBasisType() != MINT_UNDEFINED_BASIS );

  if ( this->getBasisType() == MINT_UNDEFINED_BASIS )
  {
    SLIC_WARNING( "No associated FiniteElement basis!" );
    return;
  }

  m_shapeFunctionDerivatives( xr, phidot );
}

//------------------------------------------------------------------------------
//  PRIVATE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void FiniteElement::setUp()
{
  m_jac              = new double[ m_dim*m_dim ];

  if ( !m_usingExternal )
  {
    m_xyz = new double[ m_numnodes*m_dim ];
  }

  m_phi              = new double[ m_numnodes ];
  m_phidot           = new double[ m_numnodes*m_dim ];
  m_reference_coords = new double[ m_numnodes*m_dim ];
  m_reference_center = new double[ m_dim ];
}

//------------------------------------------------------------------------------
void FiniteElement::tearDown()
{
  delete [] m_jac;

  if ( !m_usingExternal )
  {
    delete [] m_xyz;
  }

  delete [] m_phi;
  delete [] m_phidot;
  delete [] m_reference_coords;
  delete [] m_reference_center;
}

//------------------------------------------------------------------------------
void FiniteElement::getCellCoords( const Mesh* m, int cellIdx )
{
  SLIC_ASSERT(  m != AXOM_NULLPTR );
  SLIC_ASSERT(  (cellIdx >= 0) && (cellIdx < m->getMeshNumberOfCells() ) );

  // TODO: iron out this code
  int cell[ MINT_MAX_NUM_NODES ];
  m->getMeshCell( cellIdx, cell );

  for ( int i=0 ; i < m_numnodes ; ++i )
  {
    m->getMeshNode( cell[ i ], &m_xyz[ i*m_dim ]  );
  }

}

//------------------------------------------------------------------------------
bool FiniteElement::inReferenceElement( const double* xi, double TOL )
{
  SLIC_ASSERT( xi != AXOM_NULLPTR );

  const double LTOL = m_reference_min - TOL;
  const double HTOL = m_reference_max + TOL;

  bool is_inside = true;

  switch ( m_ctype )
  {
  case MINT_TRIANGLE:
  case MINT_TET:
  case MINT_PRISM:
  case MINT_PYRAMID:
    this->evaluateShapeFunctions( xi, m_phi );
    for ( int i=0 ; is_inside && (i < m_numdofs) ; ++i )
    {
      is_inside = is_inside && (m_phi[i] > LTOL) && (m_phi[i] < HTOL);
    }
    break;
  default:
    for ( int i=0 ; is_inside && i < m_dim ; ++i )
    {
      is_inside = is_inside && ( xi[i] > LTOL ) && ( xi[i] < HTOL );
    }
  } // END switch

  return is_inside;
}

} /* namespace mint */
} /* namespace axom */
