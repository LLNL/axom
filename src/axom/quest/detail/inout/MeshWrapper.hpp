// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctree_detail.hpp
 *
 * \brief Defines helper classes for the InOutOctree.
 */

#ifndef INOUT_OCTREE_MESH_WRAPPER__HXX_
#define INOUT_OCTREE_MESH_WRAPPER__HXX_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "fmt/fmt.hpp"

namespace axom
{
namespace quest
{
/**
   * \brief A utility class that wraps the access to the mesh data
   *
   * This class helps separate the specifics of accessing the underlying mesh
   * for an InOutOctree. It is customized for unstructured Triangle meshes,
   * but we will later want to apply the InOutOctree to other mesh types,
   * e.g. Segments in 2D, Bilinear quads in 3D.
   *
   * We can later specialize this class for triangle meshes and implement other
   * customized mesh wrappers.
   */
class MeshWrapper
{
public:
  static constexpr int DIM = 3;

  using VertexIndex = axom::IndexType;
  using TriangleIndex = axom::IndexType;
  using SurfaceMesh = mint::Mesh;

  using SpacePt = primal::Point<double, DIM>;
  using GeometricBoundingBox = primal::BoundingBox<double, DIM>;

  /** \brief A vertex index to indicate that there is no associated vertex */
  static const VertexIndex NO_VERTEX = -1;

  /** \brief A vertex index to indicate that there is no associated vertex */
  static const TriangleIndex NO_TRIANGLE = -1;

  /** \brief A constant for the number of boundary vertices in a triangle */
  static const int NUM_TRI_VERTS = 3;

  using SpaceTriangle = primal::Triangle<double, DIM>;

  using MeshVertexSet = slam::PositionSet<>;
  using MeshElementSet = slam::PositionSet<>;

  using VertexIndexMap = slam::Map<slam::Set<VertexIndex>, VertexIndex>;
  using VertexPositionMap = slam::Map<slam::Set<VertexIndex>, SpacePt>;

  using STLIndirection =
    slam::policies::STLVectorIndirection<VertexIndex, VertexIndex>;
  using TVStride = slam::policies::CompileTimeStride<VertexIndex, NUM_TRI_VERTS>;
  using ConstantCardinality =
    slam::policies::ConstantCardinality<VertexIndex, TVStride>;
  using TriangleVertexRelation = slam::StaticRelation<VertexIndex,
                                                      VertexIndex,
                                                      ConstantCardinality,
                                                      STLIndirection,
                                                      MeshElementSet,
                                                      MeshVertexSet>;
  using TriVertIndices = typename TriangleVertexRelation::RelationSubset;

public:
  /** \brief Constructor for a mesh wrapper */
  MeshWrapper(SurfaceMesh*& meshPtr)
    : m_surfaceMesh(meshPtr)
    , m_vertexSet(0)
    , m_elementSet(0)
    , m_vertexPositions(&m_vertexSet)
    , m_triangleToVertexRelation()
    , m_meshWasReindexed(false)
  { }

  /** Const accessor to the vertex set of the wrapped surface mesh */
  const MeshVertexSet& vertexSet() const { return m_vertexSet; }

  /** Const accessor to the element set of the wrapped surface mesh */
  const MeshElementSet& elementSet() const { return m_elementSet; }

  /** Accessor to the vertex set of the wrapped surface mesh */
  MeshVertexSet& vertexSet() { return m_vertexSet; }

  /** Accessor to the element set of the wrapped surface mesh */
  MeshElementSet& elementSet() { return m_elementSet; }

  /** Accessor for the number of vertices in the wrapped surface mesh */
  int numMeshVertices() const
  {
    if(m_meshWasReindexed)
      return m_vertexSet.size();
    else
      return m_surfaceMesh->getNumberOfNodes();
  }

  /** Accessor for the number of elements in the wrapped surface mesh */
  int numMeshElements() const
  {
    if(m_meshWasReindexed)
      return m_elementSet.size();
    else
      return m_surfaceMesh->getNumberOfCells();
  }

  /** Predicate to determine if the wrapped surface mesh has been reindexed */
  bool meshWasReindexed() const { return m_meshWasReindexed; }

  /**
     * \brief Helper function to retrieve the position of the vertex from the
     * mesh
     * \param idx The index of the vertex within the surface mesh
     */
  SpacePt getMeshVertexPosition(VertexIndex idx) const
  {
    if(m_meshWasReindexed)
    {
      return m_vertexPositions[idx];
    }
    else
    {
      SLIC_ASSERT(m_surfaceMesh->getDimension() == SpacePt::dimension());

      SpacePt pt;
      m_surfaceMesh->getNode(idx, pt.data());

      return pt;
    }
  }

  /**
     * \brief Returns the spatial position of the vertex
     *  with index idx of the wrapped surface mesh
     */
  const SpacePt& vertexPosition(VertexIndex idx) const
  {
    return m_vertexPositions[idx];
  }

  /**
     * \brief Returns the indices of the boundary vertices of the element
     * of the wrapped surface mesh with the given index
     *
     * \param idx The index of an element within the surface mesh
     */
  TriVertIndices triangleVertexIndices(TriangleIndex idx) const
  {
    return m_triangleToVertexRelation[idx];
  }

  /**
     * \brief Helper function to compute the bounding box of a triangle
     *
     * \param idx The triangle's index within the surface mesh
     */
  GeometricBoundingBox triangleBoundingBox(TriangleIndex idx) const
  {
    // Get the ids of the verts bounding this triangle
    TriVertIndices vertIds = triangleVertexIndices(idx);

    GeometricBoundingBox bb(vertexPosition(vertIds[0]));
    bb.addPoint(vertexPosition(vertIds[1]));
    bb.addPoint(vertexPosition(vertIds[2]));

    return bb;
  }

  /**
     * \brief Utility function to retrieve the positions of the triangle's
     * vertices
     *
     * \return A triangle instance whose vertices are positioned in space
     */
  SpaceTriangle trianglePositions(TriangleIndex idx) const
  {
    TriVertIndices verts = triangleVertexIndices(idx);
    return SpaceTriangle(vertexPosition(verts[0]),
                         vertexPosition(verts[1]),
                         vertexPosition(verts[2]));
  }

  /**
     * \brief Checks whether the indexed triangle contains a reference to the
     * given vertex
     */
  bool incidentInVertex(const TriVertIndices& triVerts, VertexIndex vIdx) const
  {
    return (triVerts[0] == vIdx) || (triVerts[1] == vIdx) ||
      (triVerts[2] == vIdx);
  }

  /**
     * \brief Finds the index of a vertex in triangle t1 that is not in triangle
     * t0
     *
     * \param t0 The index of the first triangle
     * \param t1 The index of the second triangle
     * \pre t0 and t1 must be distinct triangles
     * \return The index of a vertex in t1 that is not in t0,
     *         (NO_VERTEX if one does not exist)
     */
  VertexIndex distinctVertex(TriangleIndex t0, TriangleIndex t1) const
  {
    SLIC_ASSERT_MSG(
      t0 != t1,
      "Expected two different triangle indices in "
        << "quest::MeshWrapper::findDistinctVertex, got triangle index " << t0
        << " twice.");

    // Find a vertex from the local surface that is not incident in the first
    // triangle
    TriVertIndices tvRel0 = triangleVertexIndices(t0);
    TriVertIndices tvRel1 = triangleVertexIndices(t1);

    for(int i = 0; i < NUM_TRI_VERTS; ++i)
      if(!incidentInVertex(tvRel0, tvRel1[i])) return tvRel1[i];

    SLIC_ASSERT_MSG(
      false,
      fmt::format("There should be a vertex in triangle {} that was not in {}",
                  t1,
                  t0));
    return NO_VERTEX;
  }

  /**
     * \brief Determine if the two given triangles have a vertex in common
     *
     * \param t0 The index of the first triangle
     * \param t1 The index of the second triangle
     * \param [out] sharedVert The index of the shared vertex, if it exists
     * \return true if the two triangles have a vertex
     * in common (returned in sharedVert), false otherwise
     */
  bool haveSharedVertex(TriangleIndex t0,
                        TriangleIndex t1,
                        VertexIndex& sharedVert) const
  {
    // There are two triangles -- check that they have at least one common
    // vertex
    TriVertIndices tvRel0 = triangleVertexIndices(t0);
    TriVertIndices tvRel1 = triangleVertexIndices(t1);

    for(int i = 0; i < NUM_TRI_VERTS; ++i)
    {
      if(incidentInVertex(tvRel0, tvRel1[i]))
      {
        sharedVert = tvRel1[i];
        return true;
      }
    }
    return false;
  }

  /**
     * \brief Determine if the three given triangles have a vertex in common
     *
     * \param t0 The index of the first triangle
     * \param t1 The index of the second triangle
     * \param t2 The index of the second triangle
     * \param [out] sharedVert The index of the shared vertex, if it exists
     * \return true if the three triangles have a vertex in common
     *  (returned in sharedVert), false otherwise
     */
  bool haveSharedVertex(TriangleIndex t0,
                        TriangleIndex t1,
                        TriangleIndex t2,
                        VertexIndex& sharedVert) const
  {
    TriVertIndices t0Verts = triangleVertexIndices(t0);
    TriVertIndices t1Verts = triangleVertexIndices(t1);
    TriVertIndices t2Verts = triangleVertexIndices(t2);

    for(int i = 0; i < NUM_TRI_VERTS; ++i)
    {
      // check if a vertex from the third triangle is in the first and second
      if(incidentInVertex(t0Verts, t2Verts[i]) &&
         incidentInVertex(t1Verts, t2Verts[i]))
      {
        sharedVert = t2Verts[i];
        return true;
      }
    }

    return false;
  }

  /**
     * \brief Reindexes the mesh vertices and triangle indices using the given
     * map
     *
     * \param numVertices The number of vertices in the new mesh
     * \param vertexIndexMap A mapping from the old vertex indices to the new
     * ones
     * \note This step clears out the original mesh,
     * which can be reconstructed using the regenerateSurfaceMesh() function
     */
  void reindexMesh(int numVertices, const VertexIndexMap& vertexIndexMap)
  {
    // Create a vertex set on the new vertices and grab coordinates from the
    // old ones
    m_vertexSet = MeshVertexSet(numVertices);
    m_vertexPositions = VertexPositionMap(&m_vertexSet);

    int numOrigVertices = numMeshVertices();
    for(int i = 0; i < numOrigVertices; ++i)
    {
      const VertexIndex& vInd = vertexIndexMap[i];
      m_vertexPositions[vInd] = getMeshVertexPosition(i);
    }

    // Update the vertex IDs of the triangles to the new vertices
    // and create a SLAM relation on these
    int numOrigTris = numMeshElements();

    m_tv_data.clear();
    m_tv_data.reserve(NUM_TRI_VERTS * numOrigTris);
    for(axom::IndexType i = 0; i < numOrigTris; ++i)
    {
      // Grab relation from mesh
      using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
      axom::IndexType* vertIds =
        static_cast<UMesh*>(m_surfaceMesh)->getCellNodeIDs(i);

      // Remap the vertex IDs
      for(int j = 0; j < NUM_TRI_VERTS; ++j)
        vertIds[j] = vertexIndexMap[vertIds[j]];

      // Add to relation if not degenerate triangles
      // (namely, we need 3 unique vertex IDs)
      if((vertIds[0] != vertIds[1]) && (vertIds[1] != vertIds[2]) &&
         (vertIds[2] != vertIds[0]))
      {
        m_tv_data.push_back(vertIds[0]);
        m_tv_data.push_back(vertIds[1]);
        m_tv_data.push_back(vertIds[2]);
      }
    }

    m_elementSet =
      MeshElementSet(static_cast<int>(m_tv_data.size()) / NUM_TRI_VERTS);
    m_triangleToVertexRelation =
      TriangleVertexRelation(&m_elementSet, &m_vertexSet);
    m_triangleToVertexRelation.bindIndices(static_cast<int>(m_tv_data.size()),
                                           &m_tv_data);

    // Delete old mesh, and NULL its pointer
    delete m_surfaceMesh;
    m_surfaceMesh = nullptr;

    m_meshWasReindexed = true;
  }

  /** \brief Add the vertex positions and triangle boundary
     *   relations to the surface mesh */
  void regenerateSurfaceMesh()
  {
    if(m_surfaceMesh != nullptr)
    {
      delete m_surfaceMesh;
      m_surfaceMesh = nullptr;
    }

    using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
    UMesh* triMesh =
      new UMesh(3, mint::TRIANGLE, m_vertexSet.size(), m_elementSet.size());

    // Add vertices to the mesh (i.e. vertex positions)
    for(int i = 0; i < m_vertexSet.size(); ++i)
    {
      const SpacePt& pt = vertexPosition(i);
      triMesh->appendNode(pt[0], pt[1], pt[2]);
    }

    // Add triangles to the mesh (i.e. boundary vertices)
    for(int i = 0; i < m_elementSet.size(); ++i)
    {
      const TriangleIndex* tv = &triangleVertexIndices(i)[0];
      triMesh->appendCell(tv);
    }

    m_surfaceMesh = triMesh;
  }

private:
  SurfaceMesh*& m_surfaceMesh;  // ref to pointer to allow changing the mesh

  MeshVertexSet m_vertexSet;
  MeshElementSet m_elementSet;
  VertexPositionMap m_vertexPositions;

  std::vector<VertexIndex> m_tv_data;
  TriangleVertexRelation m_triangleToVertexRelation;

  bool m_meshWasReindexed;
};

}  // namespace quest
}  // namespace axom
#endif  // INOUT_OCTREE_MESH_WRAPPER__HXX_
