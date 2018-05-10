/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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

/**
 * \file IA.hpp
 *
 * \brief Contains the header information of IA class
 */

#ifndef SLAM_IA_H_
#define SLAM_IA_H_

//#include <cstddef>

#include "axom/Macros.hpp"
#include "axom_utils/Utilities.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "primal/Point.hpp"
#include "primal/Triangle.hpp"

#include "slam/IndirectionPolicies.hpp"
#include "slam/CardinalityPolicies.hpp"
#include "slam/Set.hpp"
#include "slam/RangeSet.hpp"
#include "slam/DynamicSet.hpp"
#include "slam/StaticRelation.hpp"
#include "slam/DynamicConstantRelation.hpp"


#include "slam/Utilities.hpp"
#include "slam/Map.hpp"
#include "slam/DynamicMap.hpp"
#include "slam/FieldRegistry.hpp"


namespace axom
{
namespace slam
{

/**
 * \class IA
 *
 * \brief Implements the Indexed Data structure with Adjacencies (IA)
 *
 * \detail The IAMesh class is an adjacency-based topological mesh data
 * structure for simplicial complexes.
 * \tparam TOPOLOGICAL_DIMENSION The dimension of the
 * simplices (2=triangle, 3=tetrahedra).
 * \tparam SPATIAL_DIMENSION The dimension of the embedding space.
 */

template<
  unsigned int TOPOLOGICAL_DIMENSION = 2,
  unsigned int SPATIAL_DIMENSION = 2,
  typename PointType = axom::slam::util::Point3<double>
  >
class IAMesh
{
public:

  enum
  {
    COORDS_PER_VERT = SPATIAL_DIMENSION,        //2 or 3 dimensional space
    VERTS_PER_ELEM = TOPOLOGICAL_DIMENSION + 1  //3 = triangle, 4 = tetrahedron
  };

  typedef int IndexType;
  typedef double DataType;

  typedef std::vector<IndexType>    IndexListType;

  typedef PointType Point;
  typedef IndexType PositionType;

  /// types for set
  typedef axom::slam::DynamicSet<>
    VertexSet;
  typedef axom::slam::DynamicSet<>
    ElementSet;

  /// types for relations
  typedef axom::slam::policies::
    STLVectorIndirection<PositionType,PositionType>  STLIndirection;
  typedef axom::slam::policies::
    VariableCardinality<PositionType,STLIndirection> VariableCardinality;

  typedef axom::slam::policies::
    CompileTimeStride<PositionType,VERTS_PER_ELEM>  EVStride;
  typedef axom::slam::policies::
    CompileTimeStride<PositionType,1 >              EEStride;
  typedef axom::slam::policies::
    ConstantCardinality<PositionType,EVStride>      ConstantCardinalityZ;
  typedef axom::slam::policies::
    ConstantCardinality<PositionType,EEStride>      ConstantCardinality1;
  typedef axom::slam::
    DynamicConstantRelation<ConstantCardinalityZ>   ElementToVertexRelation;
  typedef axom::slam::
    DynamicConstantRelation<ConstantCardinality1>   VertexToOneElementRelation;
  typedef axom::slam::
    DynamicConstantRelation<ConstantCardinalityZ>   ElementToElementRelation;

  /// types for maps
  typedef axom::slam::DynamicMap< Point >           PositionMap;

  //typedef axom::slam::Map< DataType >             VertexField;
  //typedef axom::slam::Map< DataType >             ElementField;

  typedef axom::slam::FieldRegistry<int>            IndexBuf;
  IndexBuf index_buffer;

  VertexSet vertex_set;                //Set of vertices
  ElementSet element_set;              //Set of elements

  ElementToVertexRelation ev_rel;      //Element to vertex relation.
  VertexToOneElementRelation ve_rel;   //Vertex to one of the element relation.
  ElementToElementRelation ee_rel;     //Element to neighboring element relation
  PositionMap vcoord_map;              //map of coordinates per vertex.

public:
  /**
   * \brief Default Constructor for an empty mesh
   */
  IAMesh();

  /**
   * \brief Construct an IA mesh with the given point coordinate and vertex
   * indices for each elements
   */
  IAMesh(std::vector<double>& points, std::vector<int>& ev_vec);

  /**
   * \brief Copy constructor
   */
  IAMesh(const IAMesh& );

  IAMesh&    operator= (const IAMesh & );

  /**
   * \brief check that the mesh data stored is valid.
   *
   * \note Valid meshes are not necesarily manifold.
   */
  bool       isValid(bool verboseOutput = false) const;

  /**
   * \brief return true if the mesh data is a pure pseudo-manifold simplicial
   * complexes with boundary.
   *
   * A pure pseudo-manifold simplicial complexes with boundary means:
   * Pseudo-manifold indicates that each facet (of dimension D-1) is
   * incident in one or two D-dimensional elements,
   * and that the mesh is D-connected, I.e. that you can traverse
   * from any element of the mesh to any other element
   * through the element facets.
   * Manifold would indicate that the link of every vertex is a (combinatorial)
   * sphere or a disk.
   * Pure indicates that the top elements of the mesh (those not on the boundary
   * of any other element) all have the same dimension.
   */
  bool       isManifold(bool verboseOutput = false) const;

  /**
   * \brief Given an element index, return a list incident vertices.
   *
   * \note If the index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexListType  getVerticesInElement(IndexType element_idx) const;

  /**
   * \brief Given a vertex index, return a list of inident elements
   *
   * \note If the index is invalid or out of bounds,
   * an empty list is returned. This function may return incorrect/partial
   * results if the mesh is not a manifold mesh.
   */
  IndexListType  getElementsWithVertex(IndexType vertex_idx) const;

  /**
   * \brief Given an element index, and a face index i,
   * return a list of the vertices on the ith face of that element
   *
   * \detail If the element index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexListType  getElementFace(IndexType element_idx,
                                IndexType face_idx) const;

  /**
   * \brief Given an element index, return a list of adjacent elements
   *
   * \detail If the index is invalid or out of bounds,
   * an empty list is returned.
   */
  IndexListType  getElementNeighbors(IndexType element_idx) const;

  /**
   * \brief Given a vertex index, return its Point coordinate.
   */
  const Point &  getVertexPoint(IndexType vertex_idx) const;

  /**
   * \brief Return true if the mesh has no vertex or element
   */
  bool           isEmpty() const;

  /**
   * \brief Returns the number of elements in the mesh.
   */
  IndexType  getNumberOfElements() const
  {
    return element_set.numberOfValidEntries();
  }

  /**
   * \brief Returns the number of vertices in the mesh.
   */
  IndexType  getNumberOfVertices() const
  {
    return vertex_set.numberOfValidEntries();
  }

  /**
   * \brief Returns true if the element indexed is a valid element.
   *
   * An element index is valid when the element is not deleted.
   */
  bool       isValidElementEntry(IndexType element_idx) const
  {
    return element_set.isValidEntry(element_idx);
  }

  /**
   * \brief Returns true if the vertex indexed is a valid vertex.
   *
   * An vertex index is valid when the vertex is not deleted.
   */
  bool       isValidVertexEntry(IndexType vertex_idx) const
  {
    return vertex_set.isValidEntry(vertex_idx);
  }

  /**
   * \brief Fix the element neighbor relation,
   *
   * \detail Given a vertex index and a list of all the elements that connects
   * to that vertex, fix the element->element relation data.
   * Sometimes when modifying the mesh, when the mesh becomes non-manifold.
   * Adding elements may result in incorrect element->element data.
   */
  void       fixVertexNeighborhood(IndexType vertex_idx,
                                   const std::vector<IndexType>& new_elements);

  /**
   * \brief Return a valid element index
   */
  IndexType  getValidElementIndex() const
  {
    for( int i = element_set.size() - 1 ; true ; i--)
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
   * \param p The coordinates of the point for a new vertex
   */
  IndexType  addVertex(const Point& p);

  /**
   * \brief Add a (triangular) element to the mesh.
   * @param n0 - first vertex id of the element
   * @param n1 - second vertex id of the element
   * @param n2 - third vertex id of the element
   */
  IndexType  addElement(IndexType n0, IndexType n1, IndexType n2);

  /**
   * \brief Add a (tetrahedral) element to the mesh.
   * @param n0 - first vertex id of the element
   * @param n1 - second vertex id of the element
   * @param n2 - third vertex id of the element
   * @param n3 - fourth vertex id of the element
   */
  IndexType  addElement(IndexType n0, IndexType n1, IndexType n2, IndexType n3);

  /**
   * \brief Add an element to the mesh.
   *
   * \param nlist - A pointer to the vertex indices of the new element.
   * The array size should be at least VERTS_PER_ELEM
   */
  IndexType  addElement(const IndexType* nlist);

  /**
   * \brief Removes an element
   *
   * \detail If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   *
   * \param element_idx The index of the element to remove.
   * \warning Removing an element could make one of its vertices non-manifold
   * (its link will have more than one connected components,
   * IA only retain 1 element reference per vertex).
   */
  void       removeElement(IndexType element_idx);

  /**
   * \brief Removes a vertex
   *
   * \details If the index is invalid or out of bounds,
   * no changes are made to the mesh.
   * \param vertex_idx The index of the vertex to remove.
   */
  void       removeVertex(IndexType vertex_idx);

  /**
   * \brief Removes all the invalid entries in the mesh and reduce memory used
   * \details This function may invalidates all indices in user code.
   */
  void       compact();

  /**
   * \brief Prints the IA mesh structure, for debug purpose.
   */
  void       print_all() const;

private:
  //helper functions to find element to element relations
  typedef PositionType ElementIndexType;
  typedef PositionType FaceIndexType;
  typedef std::pair<ElementIndexType,FaceIndexType> ElementAndFaceIdxType;
  typedef std::map<IndexListType, ElementAndFaceIdxType>       V2EMapType;

  /**
   * \brief Helper function to find adjacent elements when adding a new element
   */
  ElementAndFaceIdxType   ElemNbrFinder( V2EMapType &, IndexType element_i,
                                         IndexType side_i );


}; //end class IAMesh

} // end namespace slam
} // end namespace axom

#include "slam/IA_impl.hpp"

#endif //  SLAM_IA_H_
