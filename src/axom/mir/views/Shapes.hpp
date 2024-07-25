// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_SHAPES_HPP_
#define AXOM_MIR_VIEWS_SHAPES_HPP_

#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

#include <iostream>

namespace axom
{
namespace mir
{
namespace views
{

// Shape ids. These are used to identify shapes. These are used as indices
// in bit fields in some algorithms.
enum
{
   Point_ShapeID = 0,
   Line_ShapeID = 1,
   Tri_ShapeID = 2,
   Quad_ShapeID = 3,
   Polygon_ShapeID = 4,
   Tet_ShapeID = 5,
   Pyramid_ShapeID = 6,
   Wedge_ShapeID = 7,
   Hex_ShapeID = 8,
   Polyhedron_ShapeID = 9,
   Mixed_ShapeID = 10,

   Invalid_ShapeID = 20
};

// TODO: PointTraits


/*

  0*-----------* 1

 */
struct LineTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Line_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 1; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static IndexType faces[][2] = {{0, 1}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int /*edgeIndex*/)
  {
    return axom::StackArray<IndexType, 2>{0,1};
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "line"; }
};

/*
  2*
   |\
   | \
   |  \
   |   \
   |    \
  0*-----* 1

 */
struct TriTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Tri_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static IndexType faces[][3] = {{0, 1, 2}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tri"; }
};

/*
  3*-----------* 2
   |           |
   |           |
   |           |
   |           |
   |           |
  0*-----------* 1

 */
struct QuadTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Quad_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static IndexType faces[][4] = {{0, 1, 2, 3}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,3}, {3,0}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "quad"; }
};

/*
      3
      *
     /|\         face 0: 0,2,1
    / | \        face 1: 0,1,3
   /  |  \       face 2: 1,2,3
  /   |   \      face 3: 2,0,3
0*----|----* 2
  \   |   /      edge 0: 0,1
   \  |  /       edge 1: 1,2
    \ | /        edge 2: 2,0
     \|/         edge 3: 0,3
      *          edge 4: 1,3
       1         edge 5: 2,3

 */
struct TetTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Tet_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static IndexType faces[][3] = {{0,1,3}, {1,2,3}, {2,0,3}, {0,2,1}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 0},{0, 3}, {1, 3}, {2, 3}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "tet"; }
};

/*

3*-----------* 2  face 0: 3,2,1,0
 |\         /|    face 1: 0,1,4
 | \       / |    face 2: 1,2,4
 |  \     /  |    face 3: 2,3,4
 |   \   /   |    face 4: 3,0,4
 |    \ /    |
 |     * 4   |    edge 0: 0,1
 |    / \    |    edge 1: 1,2
 |   /   \   |    edge 2: 2,3
 |  /     \  |    edge 3: 3,0
 | /       \ |    edge 4: 0,4
 |/         \|    edge 5: 1,4
0*-----------* 1  edge 6: 2,4
                  edge 7: 3,4

 */
struct PyramidTraits
{
  AXOM_HOST_DEVICE constexpr static int id() { return Pyramid_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex) { return faceIndex == 0 ? 4 : 3; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static int faces[][4] = {{3,2,1,0}, {0,1,4,-1}, {1,2,4,-1}, {2,3,4,-1}, {3,0,4,-1}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "pyramid"; }
};

/*

3*---------* 5  face 0: 0,2,1
 |\       /|    face 1: 3,4,5
 | \     / |    face 2: 0,1,4,3
 |  \   /  |    face 3: 1,2,5,4
 |   \ /   |    face 4: 2,0,3,5
 |    *4   |
 |    |    |    edge 0: 0,1
0*----|----* 2  edge 1: 1,2
  \   |   /     edge 2: 2,0
   \  |  /      edge 3: 3,4
    \ | /       edge 4: 4,5
     \|/        edge 5: 5,3
      * 1       edge 6: 0,3
                edge 7: 1,4
                edge 8: 2,3

 */
struct WedgeTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Wedge_ShapeID; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int faceIndex) { return (faceIndex == 0 || faceIndex == 1) ? 3 : 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 5; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 9; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static int faces[][4] = {{0,2,1,-1}, {3,4,5,-1}, {0,1,4,3}, {1,2,5,4}, {2,0,3,5}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0,1}, {1,2}, {2,0}, {3,4}, {4,5}, {5,3}, {0,3}, {1,4}, {2,3}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "wedge"; }
};

/*
      4*------------* 7
      /|           /|     
     / |          / |
    /  |         /  |
  5*------------*6  |
   |   |        |   |
   |   |        |   |
   |  0*--------|---* 3
   |  /         |  /
   | /          | /
   |/           |/
   *------------*
   1            2

 */
struct HexTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Hex_ShapeID; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return false; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 3; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodes() { return 8; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfNodesInFace(int /*faceIndex*/) { return 4; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 4; }

  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 6; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfEdges() { return 12; }
  AXOM_HOST_DEVICE constexpr static IndexType zoneOffset(int zoneIndex) { return numberOfNodes() * zoneIndex; }

  constexpr static IndexType faces[][4] = {{3, 0, 4, 7}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

  AXOM_HOST_DEVICE constexpr static axom::StackArray<IndexType, 2> getEdge(int edgeIndex)
  {
    const axom::StackArray<IndexType, 2> edges[] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
    return edges[edgeIndex];
  }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "hex"; }
};

/*
  
n-1 *-... * 2
    |     |
    |     |
  0 *-----* 1

 */
struct PolygonTraits
{
  AXOM_HOST_DEVICE constexpr static int  id() { return Polygon_ShapeID; }
  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE constexpr static IndexType dimension() { return 2; }
  AXOM_HOST_DEVICE constexpr static IndexType numberOfFaces() { return 1; }
  AXOM_HOST_DEVICE constexpr static IndexType maxNodesInFace() { return 20; }
  AXOM_HOST_DEVICE constexpr static const char *name() { return "polygon"; }

};

template <typename ConnType>
struct PolygonShape : public PolygonTraits
{
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE PolygonShape(const ConnectivityView &ids) : m_idsView(ids)
  {
  }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE const ConnectivityView &getIds() const { return m_idsView; }

  /**
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   *
   * \return An array view (wrapping m_faceIds) that contains the ids for the face.
   */
  AXOM_HOST_DEVICE ConnectivityView getFace(int /*faceIndex*/) const
  {
    return m_idsView;
  }

  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex) const
  {
    const auto p0 = edgeIndex % m_idsView.size();
    const auto p1 = (edgeIndex + 1) % m_idsView.size();
    return axom::StackArray<IndexType, 2>{p0, p1};
  }

private:
  ConnectivityView m_idsView;
};

/**
 * \brief This class extends the ShapeTraits with object state so it can represent a zone.
 */
template <typename ShapeTraits, typename ConnType>
struct Shape : public ShapeTraits
{
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief Construct a shape.
   */
  AXOM_HOST_DEVICE Shape(const ConnectivityView &ids) : m_idsView(ids), m_faceIds()
  {
    SLIC_ASSERT(m_idsView.size() == ShapeTraits::numberOfNodes());
  }

  /**
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityType getId(size_t index) const { return m_idsView[index]; }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE const ConnectivityView &getIds() const { return m_idsView; }

  /**
   * \brief Get the unique ids that make up this shape. For basic shapes, assume they are unique.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityView getUniqueIds() const { return m_idsView; }

  /**
   * \brief Get the ids for the requested face.
   *
   * \param faceIndex The index of the desired face.
   *
   * \return An array view (wrapping m_faceIds) that contains the ids for the face.
   */
  AXOM_HOST_DEVICE
  ConnectivityView
  getFace(int faceIndex) const
  {
    if constexpr(ShapeTraits::dimension() == 2)
      return m_idsView;
    else
    {  
      const auto nnodes = ShapeTraits::numberOfNodesInFace(faceIndex);
      for(IndexType i = 0; i < nnodes; i++)
        m_faceIds[i] = m_idsView[ShapeTraits::faces[faceIndex][i]];
      return ConnectivityView(m_faceIds.m_data, nnodes);
    }
  }
private:
  ConnectivityView m_idsView;
  mutable axom::StackArray<ConnectivityType, ShapeTraits::maxNodesInFace()> m_faceIds;
};

// Make some concrete shape classes based on the shape traits.
template <typename ConnectivityType>
using LineShape = Shape<LineTraits, ConnectivityType>;

template <typename ConnectivityType>
using TriShape = Shape<TriTraits, ConnectivityType>;

template <typename ConnectivityType>
using QuadShape = Shape<QuadTraits, ConnectivityType>;

template <typename ConnectivityType>
using TetShape = Shape<TetTraits, ConnectivityType>;

template <typename ConnectivityType>
using PyramidShape = Shape<PyramidTraits, ConnectivityType>;

template <typename ConnectivityType>
using WedgeShape = Shape<WedgeTraits, ConnectivityType>;

template <typename ConnectivityType>
using HexShape = Shape<HexTraits, ConnectivityType>;


/**
 * \brief This is a shape that can act as any of the other shapes.
 *
 * \note This is a substitute for polymorphism so we can run on device.
 */
template <typename ConnType>
struct VariableShape
{
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief Constructor
   *
   * \param shapeId The shape id that describes the points.
   * \param ids The ids that describe the shape.
   */
  AXOM_HOST_DEVICE
  VariableShape(int shapeId, const ConnectivityView &ids) : m_shapeId(shapeId), m_idsView(ids)
  {
    //SLIC_ASSERT(shapeId >= Point_ShapeID && shapeID <= Hex_ShapeID);
  }

  /**
   * \brief Returns the shape id of the actual shape represented by the variable shape.
   * \return The actual shape represented.
   */
  AXOM_HOST_DEVICE int id() const { return m_shapeId; }

  AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return false; }
  AXOM_HOST_DEVICE constexpr static bool is_variable_size() { return true; }

  AXOM_HOST_DEVICE IndexType dimension() const
  {
    IndexType dim = 2;
    switch(m_shapeId)
    {
    case Line_ShapeID:    dim = LineTraits::dimension(); break;
    case Tri_ShapeID:     dim = TriTraits::dimension(); break;
    case Quad_ShapeID:    dim = QuadTraits::dimension(); break;
    case Polygon_ShapeID: dim = 2; break;
    case Tet_ShapeID:     dim = TetTraits::dimension(); break;
    case Pyramid_ShapeID: dim = PyramidTraits::dimension(); break;
    case Wedge_ShapeID:   dim = WedgeTraits::dimension(); break;
    case Hex_ShapeID:     dim = HexTraits::dimension(); break;
    }
    return dim;
  }

  AXOM_HOST_DEVICE IndexType numberOfNodes() const
  {
    return m_idsView.size();
  }

  AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case Line_ShapeID:    nnodes = LineTraits::numberOfNodesInFace(faceIndex); break;
    case Tri_ShapeID:     nnodes = TriTraits::numberOfNodesInFace(faceIndex); break;
    case Quad_ShapeID:    nnodes = QuadTraits::numberOfNodesInFace(faceIndex); break;
    case Polygon_ShapeID: nnodes = (faceIndex == 0) ? m_idsView.size() : 0; break;
    case Tet_ShapeID:     nnodes = TetTraits::numberOfNodesInFace(faceIndex); break;
    case Pyramid_ShapeID: nnodes = PyramidTraits::numberOfNodesInFace(faceIndex); break;
    case Wedge_ShapeID:   nnodes = WedgeTraits::numberOfNodesInFace(faceIndex); break;
    case Hex_ShapeID:     nnodes = HexTraits::numberOfNodesInFace(faceIndex); break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType maxNodesInFace() const
  {
    IndexType nnodes = 0;
    switch(m_shapeId)
    {
    case Line_ShapeID:    nnodes = LineTraits::maxNodesInFace(); break;
    case Tri_ShapeID:     nnodes = TriTraits::maxNodesInFace(); break;
    case Quad_ShapeID:    nnodes = QuadTraits::maxNodesInFace(); break;
    case Polygon_ShapeID: nnodes = PolygonShape<int>::maxNodesInFace(); break;
    case Tet_ShapeID:     nnodes = TetTraits::maxNodesInFace(); break;
    case Pyramid_ShapeID: nnodes = PyramidTraits::maxNodesInFace(); break;
    case Wedge_ShapeID:   nnodes = WedgeTraits::maxNodesInFace(); break;
    case Hex_ShapeID:     nnodes = HexTraits::maxNodesInFace(); break;
    }
    return nnodes;
  }

  AXOM_HOST_DEVICE IndexType numberOfFaces() const
  {
    IndexType nfaces = 0;
    switch(m_shapeId)
    {
    case Line_ShapeID:    nfaces = LineTraits::numberOfFaces(); break;
    case Tri_ShapeID:     nfaces = TriTraits::numberOfFaces(); break;
    case Quad_ShapeID:    nfaces = QuadTraits::numberOfFaces(); break;
    case Polygon_ShapeID: nfaces = 1; break;
    case Tet_ShapeID:     nfaces = TetTraits::numberOfFaces(); break;
    case Pyramid_ShapeID: nfaces = PyramidTraits::numberOfFaces(); break;
    case Wedge_ShapeID:   nfaces = WedgeTraits::numberOfFaces(); break;
    case Hex_ShapeID:     nfaces = HexTraits::numberOfFaces(); break;
    }
    return nfaces;
  }

  AXOM_HOST_DEVICE IndexType numberOfEdges() const
  {
    IndexType nedges = 0;
    switch(m_shapeId)
    {
    case Line_ShapeID:    nedges = LineTraits::numberOfEdges(); break;
    case Tri_ShapeID:     nedges = TriTraits::numberOfEdges(); break;
    case Quad_ShapeID:    nedges = QuadTraits::numberOfEdges(); break;
    case Polygon_ShapeID: nedges = m_idsView.size(); break;
    case Tet_ShapeID:     nedges = TetTraits::numberOfEdges(); break;
    case Pyramid_ShapeID: nedges = PyramidTraits::numberOfEdges(); break;
    case Wedge_ShapeID:   nedges = WedgeTraits::numberOfEdges(); break;
    case Hex_ShapeID:     nedges = HexTraits::numberOfEdges(); break;
    }
    return nedges;
  }

  AXOM_HOST_DEVICE axom::StackArray<IndexType, 2> getEdge(int edgeIndex) const
  {
    axom::StackArray<IndexType, 2> edge;
    switch(m_shapeId)
    {
    case Line_ShapeID:    edge = LineTraits::getEdge(edgeIndex); break;
    case Tri_ShapeID:     edge = TriTraits::getEdge(edgeIndex); break;
    case Quad_ShapeID:    edge = QuadTraits::getEdge(edgeIndex); break;
    case Polygon_ShapeID:
    {
      const auto n = m_idsView.size();
      edge[0] = edgeIndex % n;
      edge[1] = (edgeIndex + 1) % n;
      break;
    }
    case Tet_ShapeID:     edge = TetTraits::getEdge(edgeIndex); break;
    case Pyramid_ShapeID: edge = PyramidTraits::getEdge(edgeIndex); break;
    case Wedge_ShapeID:   edge = WedgeTraits::getEdge(edgeIndex); break;
    case Hex_ShapeID:     edge = HexTraits::getEdge(edgeIndex); break;
    }
    return edge;
  }

  /**
   * \brief Get a specific id that makes up this shape.
   *
   * \return The i'th id that makes up this shape.
   */
  AXOM_HOST_DEVICE ConnectivityType getId(IndexType index) const { return m_idsView[index]; }

  /**
   * \brief Get the ids that make up this shape.
   *
   * \return A view containing the ids that make up this shape.
   */
  AXOM_HOST_DEVICE const ConnectivityView &getIds() const { return m_idsView; }

  AXOM_HOST_DEVICE constexpr static const char *name() { return "mixed"; }
private:
  int m_shapeId;
  ConnectivityView m_idsView;
};

/**
 * \brief Given a shape name (matches Blueprint shape name), return the Shape id() value.
 *
 * \param name The shape name.
 *
 * \return The shape id that matches the name, or 0 if there is no match.
 */
inline int shapeNameToID(const std::string &name)
{
  int id = 0;
  if(name == LineTraits::name())
    id = Line_ShapeID;
  else if(name == TriTraits::name())
    id = Tri_ShapeID;
  else if(name == QuadTraits::name())
    id = Quad_ShapeID;
  else if(name == PolygonTraits::name())
    id = Polygon_ShapeID;
  else if(name == TetTraits::name())
    id = Tet_ShapeID;
  else if(name == PyramidTraits::name())
    id = Pyramid_ShapeID;
  else if(name == WedgeTraits::name())
    id = Wedge_ShapeID;
  else if(name == HexTraits::name())
    id = Hex_ShapeID;
  return id;
}


} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
