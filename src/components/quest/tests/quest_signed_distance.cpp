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

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

#include "mint/CellType.hpp"
#include "mint/MeshType.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/UnstructuredMesh.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Sphere.hpp"
#include "primal/Point.hpp"

#include "quest/SignedDistance.hpp"

// Google Test includes
#include "gtest/gtest.h"

// C/C++ includes
#include <cmath>

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef axom::mint::UniformMesh UniformMesh;

using axom::primal::BoundingBox;
using axom::primal::Sphere;
using axom::primal::Point;

// pi / 180
#define DEG_TO_RAD 0.01745329251

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace detail
{

/*!
 * \brief Gets a surface mesh instance for the sphere.
 * \param [in] mesh pointer to the mesh instance.
 * \pre mesh != AXOM_NULLPTR
 */
void getMesh( TriangleMesh* mesh )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  const int THETA_RES             = 25;
  const int PHI_RES               = 25;
  const double RADIUS             = 0.5;
  const double SPHERE_CENTER[ 3 ] = { 0.0, 0.0, 0.0 };
  const double theta_start        = 0;
  const double theta_end          = 360 * DEG_TO_RAD;
  const double phi_start          = 0;
  const double phi_end            = 180 * DEG_TO_RAD;

  double x[3];
  double n[3];
  int c[3];

  // North pole point
  x[0] = SPHERE_CENTER[0];
  x[1] = SPHERE_CENTER[1];
  x[2] = SPHERE_CENTER[2] + RADIUS;
  mesh->insertNode( x[0], x[1], x[2] );

  // South pole point
  x[2] = SPHERE_CENTER[2] - RADIUS;
  mesh->insertNode( x[0], x[1], x[2] );

  // Calculate spacing
  const double dphi   = ( phi_end-phi_start ) /
                        ( static_cast< double>(PHI_RES-1) );
  const double dtheta = ( theta_end-theta_start ) /
                        ( static_cast<double>(THETA_RES-1) );

  // Generate points
  for ( int i=0 ; i < THETA_RES ; ++i )
  {

    const double theta = theta_start + i*dtheta;

    for ( int j=0 ; j < PHI_RES-2 ; ++j )
    {

      const double phi = phi_start + j*dphi;
      const double radius = RADIUS * sin( phi );

      n[0] = radius * cos( theta );
      n[1] = radius * sin( theta );
      n[2] = RADIUS * cos( phi );
      x[0] = n[0] + SPHERE_CENTER[0];
      x[1] = n[1] + SPHERE_CENTER[1];
      x[2] = n[2] + SPHERE_CENTER[2];

      mesh->insertNode( x[0], x[1], x[2] );

    }  // END for all j

  } // END for all i

  const int phiResolution = PHI_RES-2; // taking in to account two pole points.
  int stride = phiResolution * THETA_RES;

  // Generate mesh connectivity around north pole
  for ( int i=0 ; i < THETA_RES ; ++i )
  {
    c[2] = phiResolution*i + /* number of poles */ 2;
    c[1] = ( phiResolution*(i+1) % stride ) + /* number of poles */ 2;
    c[0] = 0;
    mesh->insertCell( c, MINT_TRIANGLE, 3 );
  } // END for

  // Generate mesh connectivity around south pole
  int offset = PHI_RES - 1;
  for ( int i=0 ; i < THETA_RES ; ++i )
  {
    c[2] = phiResolution*i + offset;
    c[1] = ( phiResolution*(i+1) % stride ) + offset;
    c[0] = 1;
    mesh->insertCell( c, MINT_TRIANGLE, 3 );
  }

  // Generate mesh connectivity in between poles
  for ( int i=0 ; i < THETA_RES ; ++i )
  {

    for ( int j=0 ; j < PHI_RES-3 ; ++j )
    {

      c[ 0 ] = phiResolution*i + j + 2;
      c[ 1 ] = c[0] + 1;
      c[ 2 ] = ( ( phiResolution*(i+1)+j) % stride ) + 3;

      mesh->insertCell( c, MINT_TRIANGLE, 3 );

      c[ 1 ] = c[ 2 ];
      c[ 2 ] = c[ 1 ] - 1;
      mesh->insertCell( c, MINT_TRIANGLE, 3 );
    }  // END for all j

  } // END for all i

}

/*!
 * \brief Returns the bounding box of the mesh.
 * \param [in] mesh pointer to the mesh instance.
 * \return bb bounding box of the mesh
 */
BoundingBox< double,3 > getBounds( const axom::mint::Mesh* mesh )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  BoundingBox< double,3 > bb;
  Point< double,3 > pt;

  const int nnodes = mesh->getMeshNumberOfNodes();
  for ( int inode=0 ; inode < nnodes ; ++inode )
  {
    mesh->getMeshNode( inode, pt.data() );
    bb.addPoint( pt );
  }

  return( bb );
}

/*!
 * \brief Generates a uniform mesh surrounding the given triangle mesh.
 * \param [in] mesh pointer to the input mesh.
 * \param [in] umesh pointer to the uniform mesh;
 */
void getUniformMesh( const TriangleMesh* mesh, UniformMesh*& umesh )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  SLIC_ASSERT( umesh == AXOM_NULLPTR );

  const int N = 16; // number of points along each dimension

  BoundingBox< double,3 > bb = getBounds( mesh );
  bb.expand( 2.0 );

  double h[3];
  h[0] = (bb.getMax()[0]-bb.getMin()[0]) / N;
  h[1] = (bb.getMax()[1]-bb.getMin()[1]) / N;
  h[2] = (bb.getMax()[2]-bb.getMin()[2]) / N;

  int ext[6];
  ext[0] = ext[2] = ext[4] = 0;
  ext[1] = ext[3] = ext[5] = N-1;

  umesh = new UniformMesh(3,bb.getMin().data(),h,ext);
}


} /* end detail namespace */

//------------------------------------------------------------------------------
TEST( quest_signed_distance, sphere_test )
{
  const double sphere_radius   = 0.5;
  const double l1norm_expected = 6.08188;
  const double l2norm_expected = 0.123521;
  const double linf_expected   = 0.00532092;
  const double TOL             = 1.e-3;

  SLIC_INFO( "Constructing sphere mesh..." );
  TriangleMesh* surface_mesh = new TriangleMesh( 3 );
  detail::getMesh( surface_mesh );

  SLIC_INFO( "Generating uniform mesh..." );
  UniformMesh* umesh = AXOM_NULLPTR;
  detail::getUniformMesh( surface_mesh, umesh );

  const int nnodes = umesh->getNumberOfNodes();

  SLIC_INFO( "Generate BVHTree..." );
  axom::quest::SignedDistance< 3 > signed_distance( surface_mesh, 25, 10 );

  SLIC_INFO( "Compute signed distance..." );
  Sphere< double,3 > analytic_sphere( sphere_radius );
  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf   = std::numeric_limits< double >::min( );

  for ( int inode=0 ; inode < nnodes ; ++inode )
  {

    Point< double,3 > pt;
    umesh->getNode( inode, pt.data() );

    double computed = signed_distance.computeDistance( pt );
    double exact    = analytic_sphere.computeSignedDistance( pt.data() );

    // compute error
    double dx    = computed - exact;

    double absdx = std::fabs( dx );
    EXPECT_NEAR( exact, computed, 1.e-2 );

    l1norm += absdx;
    l2norm += dx*dx;

    if ( absdx > linf )
    {
      linf = absdx;
    }

  } // END for all nodes

  l2norm = std::sqrt( l2norm );

  SLIC_INFO( "l1 =" << l1norm );
  SLIC_INFO( "l2 =" << l2norm );
  SLIC_INFO( "linf = " << linf );

  EXPECT_NEAR( l1norm_expected, l1norm, TOL );
  EXPECT_NEAR( l2norm_expected, l2norm, TOL );
  EXPECT_NEAR( linf_expected, linf, TOL );

  delete surface_mesh;
  delete umesh;

  SLIC_INFO( "Done." );
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
