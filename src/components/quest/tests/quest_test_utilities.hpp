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


#include "slic/slic.hpp"

#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/orientation.hpp"

#include "mint/DataTypes.hpp"
#include "mint/Mesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/Field.hpp"
#include "mint/FieldData.hpp"

using axom::primal::Point;
using axom::primal::Triangle;

/**
 * \file
 *
 * This file contains several utility functions for testing the quest component,
 * e.g. generating random doubles and Points, creating a simple mesh of an octahedron...
 * We may later decide to move some of these into the actual component if they are deemed useful.
 */

namespace axom
{
namespace quest
{
namespace utilities
{


/**
 * \brief Simple utility function to generate a random double in the range of beg to end
 * \param beg The lower value of the range (default 0.)
 * \param end The upper value of the range (default 1.)
 * \pre If the range is zero, we set it to 1
 * \return A double in the range [beg,end]
 */
double randomDouble(double beg = 0., double end = 1.)
{
  double range = end-beg;

  if(range == 0)
    range = 1.;

  return beg + (rand() / ( RAND_MAX / range ) );
}

/**
 * \brief Simple utility to generate a Point whose entries
 * are random values in the range [beg, end]
 */
template<int DIM>
Point<double,DIM> randomSpacePt(double beg, double end)
{
  Point<double,DIM> pt;
  for(int i=0 ; i< DIM ; ++i)
    pt[i] = randomDouble(beg,end);

  return pt;
}


/**
 * \brief Simple utility to find the centroid of two points
 */
template<int DIM>
Point<double,DIM> getCentroid( const Point<double,DIM>& pt0,
                               const Point<double,DIM>& pt1)
{
  return (pt0.array() + pt1.array()) /2.;
}

/**
 * \brief Simple utility to find the centroid of three points
 */
template<int DIM>
Point<double,DIM> getCentroid( const Point<double,DIM>& pt0,
                               const Point<double,DIM>& pt1,
                               const Point<double,DIM>& pt2  )
{
  return (pt0.array() + pt1.array() + pt2.array()) /3.;
}

/**
 * \brief Utility function to generate a triangle mesh of an octahedron
 * Vertices of the octahedron are at +-i, +-j and +-k.
 * \note The caller must delete the mesh
 */
axom::mint::Mesh * make_octahedron_mesh()
{
  typedef axom::mint::localIndex VertexIndex;
  typedef Point<double, 3> SpacePt;
  typedef Triangle<double, 3> SpaceTriangle;

  enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

  // The six vertices of the octahedron
  const int NUM_VERTS = 6;
  SpacePt verts[NUM_VERTS]
    =  { SpacePt::make_point( 1., 0., 0.)
         , SpacePt::make_point(-1., 0., 0.)
         , SpacePt::make_point( 0,  1., 0.)
         , SpacePt::make_point( 0, -1., 0.)
         , SpacePt::make_point( 0,  0,  1.)
         , SpacePt::make_point( 0,  0, -1.)  };

  // The eight triangles of the octahedron
  // Explicit representation of triangle-vertex incidence relation
  // Note: We are orienting the triangles with normals pointing outside
  const int NUM_TRIS = 8;
  const int VERTS_PER_TRI = 3;
  VertexIndex tvRelation[NUM_TRIS*VERTS_PER_TRI]
    = { POS_Z, POS_X, POS_Y
        , POS_Z, POS_Y, NEG_X
        , POS_Z, NEG_X, NEG_Y
        , POS_Z, NEG_Y, POS_X
        , NEG_Z, POS_Y, POS_X
        , NEG_Z, NEG_X, POS_Y
        , NEG_Z, NEG_Y, NEG_X
        , NEG_Z, POS_X, NEG_Y };

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

  // First, confirm that all triangles have normals that point away from the origin
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
  typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
  TriangleMesh * triMesh = new TriangleMesh(3);

  // insert verts
  for(int i=0 ; i< NUM_VERTS ; ++i)
  {
    triMesh->addNode(verts[i][0], verts[i][1], verts[i][2]);
  }

  // insert triangles
  for(int i=0 ; i< NUM_TRIS ; ++i)
  {
    triMesh->addCell(&tvRelation[i*VERTS_PER_TRI], MINT_TRIANGLE);
  }

  SLIC_ASSERT( NUM_VERTS == triMesh->getMeshNumberOfNodes() );
  SLIC_ASSERT( NUM_TRIS == triMesh->getMeshNumberOfCells() );

  return triMesh;
}

/**
 * \brief Utility function to generate the triangle mesh of a tetrahedron
 * Vertices are close to, but not exactly: (0, 0, 20),
 * (-18.21, 4.88, -6.66), (4.88, -18.21, -6.66), (13.33, 13.33, -6.66)
 * \note The caller must delete the mesh
 */
axom::mint::Mesh* make_tetrahedron_mesh()
{
  typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

  TriangleMesh* surface_mesh = new TriangleMesh(3);
  surface_mesh->insertNode( -0.000003, -0.000003, 19.999999);
  surface_mesh->insertNode(-18.213671,  4.880339, -6.666668);
  surface_mesh->insertNode(  4.880339,-18.213671, -6.666668);
  surface_mesh->insertNode( 13.333334, 13.333334, -6.666663);
  int cell[3];
  cell[0] = 0;    cell[1] = 1;    cell[2] = 2;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 0;    cell[1] = 3;    cell[2] = 1;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 0;    cell[1] = 2;    cell[2] = 3;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 1;    cell[1] = 3;    cell[2] = 2;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);

  return surface_mesh;
}

} // end namespace utilities
} // end namespace quest
} // end namespace axom


#endif // QUEST_TEST_UTILITIES_HPP_
