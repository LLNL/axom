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

#ifndef QUEST_TEST_UTILITIES_HPP_
#define QUEST_TEST_UTILITIES_HPP_

#include "axom/core/utilities/Utilities.hpp" // for random_real

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/operators/orientation.hpp"

#include "axom/mint/config.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/mesh/Field.hpp"
#include "axom/mint/mesh/FieldData.hpp"

using axom::primal::Point;
using axom::primal::Triangle;

/*!
 * \file
 *
 * This file contains several utility functions to facilitate in the development
 * of unit tests in quest, e.g., generating a random points and simple
 * unstructured mesh types, such as, the surface mesh of a sphere with a
 * given center and radius and an octahedron.
 */

namespace axom
{
namespace quest
{
namespace utilities
{

// pi / 180
constexpr double DEG_TO_RAD = 0.01745329251;

/*!
 * \brief Gets a surface mesh instance for the sphere.
 *
 * \param [in,out] mesh pointer to the mesh instance, where to store the mesh
 * \param [in] SPHERE_CENTER the center of the sphere
 * \param [in] RADIUS the radius of the sphere
 * \param [in] THETA_RES the theta resolution for generating the sphere
 * \param [in] PHI_RES the phi resolution for generating the sphere
 *
 * \pre mesh != AXOM_NULLPTR
 * \pre mesh->getDimension() == 3
 * \pre mesh->getMeshType() == mint::UNSTRUCTURED_MESH
 * \pre mesh->hasMixedCellTypes() == false
 * \pre mesh->getCellType() == mint::TRIANGLE
 */
void getSphereSurfaceMesh( mint::UnstructuredMesh< mint::SINGLE_SHAPE >* mesh,
                           const double* SPHERE_CENTER,
                           double RADIUS,
                           int THETA_RES,
                           int PHI_RES )
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  SLIC_ASSERT( mesh->getDimension()==3 );
  SLIC_ASSERT( mesh->getMeshType()==mint::UNSTRUCTURED_MESH );
  SLIC_ASSERT( mesh->hasMixedCellTypes()==false );
  SLIC_ASSERT( mesh->getCellType()==mint::TRIANGLE );

  const double totalNodes = THETA_RES * PHI_RES;
  const double totalCells = (2 * THETA_RES) + (THETA_RES * PHI_RES);
  mesh->reserve( totalNodes, totalCells );

  const double theta_start = 0;
  const double theta_end   = 360 * DEG_TO_RAD;
  const double phi_start   = 0;
  const double phi_end     = 180 * DEG_TO_RAD;

  double x[3];
  double n[3];
  mint::IndexType c[3];

  // North pole point
  x[0] = SPHERE_CENTER[0];
  x[1] = SPHERE_CENTER[1];
  x[2] = SPHERE_CENTER[2] + RADIUS;
  mesh->appendNode( x[0], x[1], x[2] );

  // South pole point
  x[2] = SPHERE_CENTER[2] - RADIUS;
  mesh->appendNode( x[0], x[1], x[2] );

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

      mesh->appendNode( x[0], x[1], x[2] );

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
    mesh->appendCell( c );
  } // END for

  // Generate mesh connectivity around south pole
  int offset = PHI_RES - 1;
  for ( int i=0 ; i < THETA_RES ; ++i )
  {
    c[2] = phiResolution*i + offset;
    c[1] = ( phiResolution*(i+1) % stride ) + offset;
    c[0] = 1;
    mesh->appendCell( c );
  }

  // Generate mesh connectivity in between poles
  for ( int i=0 ; i < THETA_RES ; ++i )
  {

    for ( int j=0 ; j < PHI_RES-3 ; ++j )
    {

      c[ 0 ] = phiResolution*i + j + 2;
      c[ 1 ] = c[0] + 1;
      c[ 2 ] = ( ( phiResolution*(i+1)+j) % stride ) + 3;

      mesh->appendCell( c );

      c[ 1 ] = c[ 2 ];
      c[ 2 ] = c[ 1 ] - 1;
      mesh->appendCell( c );
    }  // END for all j
  } // END for all i
}

/*!
 * \brief Simple utility to generate a Point whose entries
 * are random values in the range [beg, end]
 */
template<int DIM>
Point<double,DIM> randomSpacePt(double beg, double end)
{
  Point<double,DIM> pt;
  for(int i=0 ; i< DIM ; ++i)
    pt[i] = axom::utilities::random_real( beg, end );

  return pt;
}


/*!
 * \brief Simple utility to find the centroid of two points
 */
template<int DIM>
Point<double,DIM> getCentroid( const Point<double,DIM>& pt0,
                               const Point<double,DIM>& pt1)
{
  return (pt0.array() + pt1.array()) /2.;
}

/*!
 * \brief Simple utility to find the centroid of three points
 */
template<int DIM>
Point<double,DIM> getCentroid( const Point<double,DIM>& pt0,
                               const Point<double,DIM>& pt1,
                               const Point<double,DIM>& pt2  )
{
  return (pt0.array() + pt1.array() + pt2.array()) /3.;
}

/*!
 * \brief Utility function to generate a triangle mesh of an octahedron
 * Vertices of the octahedron are at +-i, +-j and +-k.
 * \note The caller must delete the mesh
 */
axom::mint::Mesh* make_octahedron_mesh()
{
  typedef axom::mint::IndexType VertexIndex;
  typedef Point<double, 3> SpacePt;
  typedef Triangle<double, 3> SpaceTriangle;

  enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

  // The six vertices of the octahedron
  const int NUM_VERTS = 6;
  SpacePt verts[NUM_VERTS]
    =  { SpacePt::make_point( 1., 0., 0.),
         SpacePt::make_point(-1., 0., 0.),
         SpacePt::make_point( 0,  1., 0.),
         SpacePt::make_point( 0, -1., 0.),
         SpacePt::make_point( 0,  0,  1.),
         SpacePt::make_point( 0,  0, -1.)  };

  // The eight triangles of the octahedron
  // Explicit representation of triangle-vertex incidence relation
  // Note: We are orienting the triangles with normals pointing outside
  const int NUM_TRIS = 8;
  const int VERTS_PER_TRI = 3;
  VertexIndex tvRelation[NUM_TRIS*VERTS_PER_TRI]
    = { POS_Z, POS_X, POS_Y,
        POS_Z, POS_Y, NEG_X,
        POS_Z, NEG_X, NEG_Y,
        POS_Z, NEG_Y, POS_X,
        NEG_Z, POS_Y, POS_X,
        NEG_Z, NEG_X, POS_Y,
        NEG_Z, NEG_Y, NEG_X,
        NEG_Z, POS_X, NEG_Y };

  // Note (KW 3/2016) -- We are not currently using this
  // Explicit representation of edge-vertex incidence relation
  // Note: we don't care about the orientation here
  //const int NUM_EDGES = 12;
  //const int VERTS_PER_EDGE = 3;
  //VertexIndex evRelation[NUM_EDGES*VERTS_PER_EDGE]
  //              = { POS_Z, POS_X  // Four edges incident in +Z
  //                , POS_Z, POS_Y
  //                , POS_Z, NEG_X
  //                , POS_Z, NEG_Y
  //                , NEG_Z, POS_X  // Four edges incident in -Z
  //                , NEG_Z, POS_Y
  //                , NEG_Z, NEG_X
  //                , NEG_Z, NEG_Y
  //                , POS_Y, POS_X  // Four edges not incident in Z
  //                , NEG_Y, POS_Y
  //                , POS_Y, NEG_X
  //                , NEG_Y, NEG_Y };

  // First, confirm that all triangle normals point away from the origin
  for(int i =0 ; i < NUM_TRIS ; ++i)
  {
    int baseIndex = i*VERTS_PER_TRI;
    SpaceTriangle tri( verts[ tvRelation[ baseIndex + 0]]
                       , verts[ tvRelation[ baseIndex + 1]]
                       , verts[ tvRelation[ baseIndex + 2]] );

    SLIC_ASSERT( axom::primal::ON_NEGATIVE_SIDE ==
                 axom::primal::orientation( SpacePt(), tri) );
  }

  // Now create an unstructured triangle mesh from the two arrays
  typedef axom::mint::UnstructuredMesh< mint::SINGLE_SHAPE > UMesh;
  UMesh* triMesh = new UMesh(3, mint::TRIANGLE);

  // insert verts
  for(int i=0 ; i< NUM_VERTS ; ++i)
  {
    triMesh->appendNode(verts[i][0], verts[i][1], verts[i][2]);
  }

  // insert triangles
  for(int i=0 ; i< NUM_TRIS ; ++i)
  {
    triMesh->appendCell(&tvRelation[i*VERTS_PER_TRI]);
  }

  SLIC_ASSERT( NUM_VERTS == triMesh->getNumberOfNodes() );
  SLIC_ASSERT( NUM_TRIS == triMesh->getNumberOfCells() );

  return triMesh;
}

/*!
 * \brief Utility function to generate the triangle mesh of a tetrahedron
 * Vertices are close to, but not exactly: (0, 0, 20),
 * (-18.21, 4.88, -6.66), (4.88, -18.21, -6.66), (13.33, 13.33, -6.66)
 * \note The caller must delete the mesh
 */
axom::mint::Mesh* make_tetrahedron_mesh()
{
  typedef axom::mint::UnstructuredMesh< mint::SINGLE_SHAPE > UMesh;

  UMesh* surface_mesh = new UMesh(3, mint::TRIANGLE);
  surface_mesh->appendNode( -0.000003, -0.000003, 19.999999);
  surface_mesh->appendNode(-18.213671,  4.880339, -6.666668);
  surface_mesh->appendNode(  4.880339,-18.213671, -6.666668);
  surface_mesh->appendNode( 13.333334, 13.333334, -6.666663);
  axom::mint::IndexType cell[3];
  cell[0] = 0;    cell[1] = 1;    cell[2] = 2;
  surface_mesh->appendCell(cell);
  cell[0] = 0;    cell[1] = 3;    cell[2] = 1;
  surface_mesh->appendCell(cell);
  cell[0] = 0;    cell[1] = 2;    cell[2] = 3;
  surface_mesh->appendCell(cell);
  cell[0] = 1;    cell[1] = 3;    cell[2] = 2;
  surface_mesh->appendCell(cell);

  return surface_mesh;
}

} // end namespace utilities
} // end namespace quest
} // end namespace axom


#endif // QUEST_TEST_UTILITIES_HPP_
