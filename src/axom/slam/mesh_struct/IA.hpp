// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IA.hpp
 *
 * \brief Contains the header information of IA class
 */

#ifndef SLAM_IA_H_
#define SLAM_IA_H_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/StaticRelation.hpp"
#include "axom/slam/DynamicConstantRelation.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/DynamicMap.hpp"
#include "axom/slam/FieldRegistry.hpp"

namespace axom
{
namespace slam
{
/**
 * \class IA
 *
 * \brief Implements the Indexed mesh data structure with Adjacencies (IA)
 *
 * \tparam TDIM The topological dimension (2=triangle, 3=tetrahedra).
 * \tparam SDIM The dimension of the embedding space.
 * \tparam PointType A point type for the mesh position data
 *
 * \details The IAMesh class is an adjacency-based topological mesh data
 * structure for simplicial complexes. It encodes the boundary relation
 * from elements to their vertices, the partial coboundary relation from
 * vertices to one incident element and the adjacency relation between
 * elements along their facets (faces of dimension TDIM-1).
 *
 * PointType is required to have the following interface:
 * TODO: FILL THIS IN
 */

template <unsigned int TDIM = 2,
          unsigned int SDIM = 3,
          typename PointType = axom::slam::util::Point3<double>>
class IAMesh
{
public:
  static constexpr int COORDS_PER_VERT = SDIM;     //2 or 3 dimensional space
  static constexpr int VERTS_PER_ELEM = TDIM + 1;  //3 = triangle, 4 = tetrahedron

  using IndexType = axom::IndexType;
  using DataType = double;

  using IndexArray = std::vector<IndexType>;

  using Point = PointType;
  using PositionType = IndexType;
  using ElementType = IndexType;

  /// types for Sets
  using SetBase = slam::Set<PositionType, ElementType>;
  using VertexSet = slam::DynamicSet<PositionType, ElementType>;
  using ElementSet = slam::DynamicSet<PositionType, ElementType>;

  static constexpr int INVALID_VERTEX_INDEX = VertexSet::INVALID_ENTRY;
  static constexpr int INVALID_ELEMENT_INDEX = ElementSet::INVALID_ENTRY;

  /// types for Relations
  using STLIndirection =
    slam::policies::STLVectorIndirection<PositionType, ElementType>;
  using VariableCardinality =
    slam::policies::VariableCardinality<PositionType, STLIndirection>;

  using EVStride =
    slam::policies::CompileTimeStride<PositionType, VERTS_PER_ELEM>;
  using EEStride = slam::policies::CompileTimeStride<PositionType, 1>;
  using ConstantCardinalityZ =
    slam::policies::ConstantCardinality<PositionType, EVStride>;
  using ConstantCardinality1 =
    slam::policies::ConstantCardinality<PositionType, EEStride>;

  using ElementBoundaryRelation =
    slam::DynamicConstantRelation<PositionType, ElementType, ConstantCardinalityZ>;
  using VertexCoboundaryRelation =
    slam::DynamicConstantRelation<PositionType, ElementType, ConstantCardinality1>;
  using ElementAdjacencyRelation =
    slam::DynamicConstantRelation<PositionType, ElementType, ConstantCardinalityZ>;

  /// types for Maps
  using PositionMap = DynamicMap<VertexSet, Point>;

  //using VertexField = slam::Map< SetBase, DataType >;
  //using ElementField = slam::Map< SetBase, DataType >;

  using IndexBuf = slam::FieldRegistry<SetBase, IndexType>;
  IndexBuf index_buffer;

  VertexSet vertex_set;    //Set of vertices
  ElementSet element_set;  //Set of elements

  ElementBoundaryRelation ev_rel;   //Element to vertex relation.
  VertexCoboundaryRelation ve_rel;  //Vertex to one element partial relation.
  ElementAdjacencyRelation ee_rel;  //Element to neighboring element relation
  PositionMap vcoord_map;           //map of coordinates per vertex.

public:
  /**
   * \brief Default Constructor for an empty mesh
   */
  IAMesh();

  /**
   * \brief Construct an IA mesh with the given point coordinate and vertex
   * indices for each elements
   */
  IAMesh(std::vector<double>& points, std::vector<IndexType>& ev_vec);

  /**
   * \brief Copy constructor
   */
  IAMesh(const IAMesh&);

  IAMesh& operator=(const IAMesh&);

  /**
   * \brief check that the mesh data stored is valid.
   *
   * \note Valid meshes are not necessarily manifold.
   */
  bool isValid(bool verboseOutput = false) const;

  /**
   * \brief return true if the encoded mesh is a pure pseudo-manifold
   * simplicial complexes (with boundary)
   *
   * The D-dimensional mesh is *pure* when all top elements (those not
   * on the boundary of other elements) are D-dimensional.
   * A D-dimensional mesh is pseudo-manifold when
   * (a) each facet (i.e. face of dimension D-1) is incident in one or
   * two D-dimensional elements, and (b) the mesh is D-connected,
   * i.e. that you can traverse from any element of the mesh to any other
   * element through the element facets.
   * It is manifold when the link of every vertex is a (combinatorial)
   * sphere or a disk.
   *
   * \warning This function is only partly implemented
   */
  bool isManifold(bool verboseOutput = false) const;

  /**
   * \brief Given an element index, return a list incident vertices.
   *
   * \note If the index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexArray getVerticesInElement(IndexType element_idx) const;

  /**
   * \brief Given a vertex index, return a list of incident elements
   *
   * \note If the index is invalid or out of bounds,
   * an empty list is returned. This function may return incorrect/partial
   * results if the mesh is not a manifold mesh.
   */
  IndexArray getElementsWithVertex(IndexType vertex_idx) const;

  /**
   * \brief Given an element index, and a face index i,
   * return a list of the vertices on the i^th face of that element
   *
   * \details If the element index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexArray getElementFace(IndexType element_idx, IndexType face_idx) const;

  /**
   * \brief Given an element index, return a list of adjacent elements
   *
   * \details If the index is invalid or out of bounds,
   * an empty list is returned.
   *
   * TODO: Add what happens when neighbor elements are invalid
   */
  IndexArray getElementNeighbors(IndexType element_idx) const;

  /**
   * \brief Given a vertex index, return its Point coordinate.
   */
  const Point& getVertexPoint(IndexType vertex_idx) const;

  /**
   * \brief Return true if the mesh has no vertex or element
   */
  bool isEmpty() const;

  /**
   * \brief Returns the number of elements in the mesh.
   */
  IndexType getNumberOfElements() const
  {
    return element_set.numberOfValidEntries();
  }

  /**
   * \brief Returns the number of vertices in the mesh.
   */
  IndexType getNumberOfVertices() const
  {
    return vertex_set.numberOfValidEntries();
  }

  /**
   * \brief Returns true if the element indexed is a valid element.
   *
   * An element index is valid when the element is not deleted.
   */
  bool isValidElementEntry(IndexType element_idx) const
  {
    return element_set.isValidEntry(element_idx);
  }

  /**
   * \brief Returns true if the vertex indexed is a valid vertex.
   *
   * An vertex index is valid when the vertex is not deleted.
   */
  bool isValidVertexEntry(IndexType vertex_idx) const
  {
    return vertex_set.isValidEntry(vertex_idx);
  }

  /**
   * \brief Fix the element adjacency relation in the neighborhood of a vertex
   *
   * \details Given a vertex index and a list of all the elements incident in
   * that vertex, fix the element->element relation data.
   * Sometimes when modifying the mesh, the mesh becomes non-manifold.
   * Adding elements may result in incorrect element->element data.
   */
  void fixVertexNeighborhood(IndexType vertex_idx,
                             const std::vector<IndexType>& new_elements);

  /**
   * \brief Return a valid element index
   */
  IndexType getValidElementIndex() const
  {
    for(int i = element_set.size() - 1; true; i--)
    {
      if(isValidElementEntry(i))
      {
        return i;
      }
    }
  }

  /**
   * \brief Add a vertex to the mesh,
   *
   * \param p Coordinates for the new vertex
   */
  IndexType addVertex(const Point& p);

  /**
   * \brief Add an element to the mesh.
   *
   * This is a convenience function that only uses the first VERTS_PER_ELEM
   * vertex identifiers.
   *
   * \param n0 - first vertex id of the element
   * \param n1 - second vertex id of the element
   * \param n2 - third vertex id of the element
   * \param n3 - fourth vertex id of the element
   */
  IndexType addElement(IndexType n0,
                       IndexType n1,
                       IndexType n2 = INVALID_VERTEX_INDEX,
                       IndexType n3 = INVALID_VERTEX_INDEX);

  /**
   * \brief Add an element to the mesh.
   *
   * \param nlist - A pointer to the vertex indices of the new element.
   * The array size should be at least VERTS_PER_ELEM
   */
  IndexType addElement(const IndexType* nlist);

  /**
   * \brief Removes an element from the mesh
   *
   * \details If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   *
   * \param element_idx The index of the element to remove.
   * \warning Removing an element could make one of its vertices non-manifold
   * (its link will have more than one connected components) since the
   * IA only retains 1 element reference per vertex.
   */
  void removeElement(IndexType element_idx);

  /**
   * \brief Removes a vertex and all incident elements from the mesh
   *
   * \details If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   * \param vertex_idx The index of the vertex to remove.
   */
  void removeVertex(IndexType vertex_idx);

  /**
   * \brief Removes all the invalid entries in the mesh and reduce memory used
   * \details This function may invalidates all indices in user code.
   */
  void compact();

  /**
   * \brief Prints the IA mesh structure, for debug purpose.
   */
  void print_all() const;

private:
  //helper functions to find element to element relations
  using ElementIndexType = PositionType;
  using FaceIndexType = PositionType;
  using ElementAndFaceIdxType = std::pair<ElementIndexType, FaceIndexType>;
  using V2EMapType = std::map<IndexArray, ElementAndFaceIdxType>;

  /**
   * \brief Helper function to find adjacent elements when adding a new element
   *
   * \param vertpair_to_elem_map A map from facets encoded using a sorted tuple of
   * vertex IDs to the face encoded as a element index and a face index
   * \param element_i The index of the element whose neighbor we are seeking
   * \param side_i The local index of the desired face w.r.t. \a element_i
   */
  ElementAndFaceIdxType ElemNbrFinder(V2EMapType& vertpair_to_elem_map,
                                      IndexType element_i,
                                      IndexType side_i);

};  //end class IAMesh

}  // end namespace slam
}  // end namespace axom

#include "axom/slam/mesh_struct/IA_impl.hpp"

#endif  //  SLAM_IA_H_
