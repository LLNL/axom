// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutOctree.hpp
 *
 * \brief Defines an InOutOctree for containment queries on a surface.
 */

#ifndef INOUT_OCTREE__HXX_
#define INOUT_OCTREE__HXX_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "fmt/fmt.hpp"

#include <vector>    // For InOutLeafData triangle lists
#include <iterator>  // For back_inserter
#include <limits>    // numeric_limits traits
#include <sstream>
#include <unordered_map>

#define DEBUG_VERT_IDX -2  // 1160
#define DEBUG_TRI_IDX -2   // 1654

#define DEBUG_BLOCK_2 BlockIndex::invalid_index()
//                     BlockIndex(GridPt::make_point(32, 61, 20), 6)
#define DEBUG_BLOCK_1 BlockIndex::invalid_index()
//                     BlockIndex(GridPt::make_point(32, 60, 20), 6)

#ifndef DUMP_VTK_MESH
//    #define DUMP_VTK_MESH
#endif

#ifndef DUMP_OCTREE_INFO
//    #define DUMP_OCTREE_INFO 1
#endif

#ifndef DEBUG_OCTREE_ACTIVE
//    #define DEBUG_OCTREE_ACTIVE
#endif

#if defined(DEBUG_OCTREE_ACTIVE) && defined(AXOM_DEBUG)
  #define QUEST_OCTREE_DEBUG_LOG_IF(_cond, _msg) \
    if(_cond) SLIC_DEBUG(_msg)
#else
  #define QUEST_OCTREE_DEBUG_LOG_IF(_cond, _msg) ((void)0)
#endif

namespace axom
{
namespace quest
{
namespace detail
{
// Predeclare utility classes for validating and logging stats of InOutOctrees
template <int DIM>
class InOutOctreeStats;
template <int DIM>
class InOutOctreeValidator;
template <int DIM>
class InOutOctreeMeshDumper;

}  // end namespace detail

/**
 * \brief Compact BlockDataType for an InOutOctree
 *
 * Storage requirement is one integer per block to hold the color of a block
 * and for gray block, the index of the associated triangles
 */
class InOutBlockData
{
  // Some internal constants for keeping tracking of the associated block
  // A block is a leaf block when its m_idx is not INTERNAL_BLOCK
  // Leaf blocks can be uncolored or colored (without additional data)
  //      or m_idx be the index of the data associated with a gray block
  enum
  {
    LEAF_BLOCK_UNCOLORED = -1,
    LEAF_BLOCK_WHITE = -2,
    LEAF_BLOCK_BLACK = -3,
    INTERNAL_BLOCK = -4,
    NON_BLOCK = -5
  };

public:
  enum LeafColor
  {
    Undetermined = -2,
    White = -1,
    Gray = 0,
    Black = 1
  };

public:
  /**
   * \brief Default constructor for an InOutBlockData
   *
   * \note Default constructed instances are assumed to be leaf blocks
   */
  InOutBlockData() : m_idx(LEAF_BLOCK_UNCOLORED) { }

  /** \brief Constructor from a given index */
  explicit InOutBlockData(int dataIdx) : m_idx(dataIdx) { }

  /** \brief Copy constructor for an InOutBlockData instance */
  InOutBlockData(const InOutBlockData& other) : m_idx(other.m_idx) { }

  /** \brief Assignment operator for an InOutBlockData instance */
  InOutBlockData& operator=(const InOutBlockData& other)
  {
    this->m_idx = other.m_idx;
    return *this;
  }

public:  // API for a BlockData
  /**
   * \brief Predicate to determine if the associated block is a leaf
   *
   * \return True, if the block is a leaf, False otherwise
   */
  bool isLeaf() const { return m_idx > INTERNAL_BLOCK; }

  /** \brief Marks the associated block as internal */
  void setInternal() { m_idx = INTERNAL_BLOCK; }

  /** \brief Marks the associated block as a non-block (i.e. not in the tree) */
  void setNonBlock() { m_idx = NON_BLOCK; }

  /**
   * \brief Predicate to determine if the associated block is in the tree
   *
   * \return True, if the block is in the tree (internal or leaf), False
   * otherwise
   */
  bool isBlock() const { return m_idx != NON_BLOCK; }

public:  // Other functions
  /**
   * Clears the data associated with the block
   * \note This function is currently a no-op
   * */
  void clear()
  {
    // No-op for now -- eventually, will need to do something about the index
  }

  /**
   * Predicate to determine if the associated block has data (i.e. it is a gray
   * block)
   * \return True, if the block has data, False otherwise
   * */
  bool hasData() const { return m_idx >= 0; }

  /**
   * Returns the index of the data associated with the block
   */
  const int& dataIndex() const
  {
    //SLIC_ASSERT(hasData());
    return m_idx;
  }

  /**
   * \brief Sets the block as gray, and provides index of its associated data
   *
   * \param idx The index of the data associated with the gray leaf block
   * \pre The block must be a leaf block
   * \pre The passed in index, idx, must be a non-negative integer
   */
  void setGray(int idx)
  {
    SLIC_ASSERT(isLeaf());
    SLIC_ASSERT(idx >= 0);
    m_idx = idx;
  }

  /** Marks the block as Black (the entire domain is inside the surface) */
  void setBlack()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_BLACK;
  }

  /** Marks the block as Black (the entire domain is outside the surface) */
  void setWhite()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_WHITE;
  }

  /** Sets the data associated with the block to the given index idx */
  void setData(int idx) { m_idx = idx; }

  /** Marks the block as uncolored */
  void setUncoloredLeaf()
  {
    SLIC_ASSERT(isLeaf());
    m_idx = LEAF_BLOCK_UNCOLORED;
  }

  /**
   * \brief Find the 'color' of this LeafBlock
   *
   * 'Black' indicates that the entire block is within the surface
   * 'White' indicates that the entire block is outside the surface
   * 'Gray' indicates that the block intersects the surface geometry
   * Leaves that haven't been colored yet are 'Undetermined'
   */
  LeafColor color() const
  {
    if(hasData()) return Gray;

    switch(m_idx)
    {
    case LEAF_BLOCK_BLACK:
      return Black;
    case LEAF_BLOCK_WHITE:
      return White;
    case LEAF_BLOCK_UNCOLORED:
      return Undetermined;
    }

    SLIC_ASSERT_MSG(false, "Invalid state in InOuLeafData::color()");
    return Undetermined;
  }

  /** Predicate to determine if the associated block has a color
   * \return True if the block has a color, false otherwise
   * \sa color()
   */
  bool isColored() const { return color() != Undetermined; }

  /** Friend function to compare equality of two InOutBlockData instances  */
  friend bool operator==(const InOutBlockData& lhs, const InOutBlockData& rhs)
  {
    return lhs.m_idx == rhs.m_idx;
  }

private:
  int m_idx;
};

/**
 * Free function to print an InOutBlockData to an output stream
 * \param os The output stream to write to
 * \param iob The InOUtBlockData instance that we are writing
 */
inline std::ostream& operator<<(std::ostream& os, const InOutBlockData& iob)
{
  os << "InOutBlockData{"
     << "isLeaf: " << (iob.isLeaf() ? "yes" : "no");

  bool showData = true;

  if(iob.isLeaf())
  {
    os << ", color: ";
    switch(iob.color())
    {
    case InOutBlockData::Gray:
      os << "Gray";
      break;
    case InOutBlockData::White:
      os << "White";
      showData = false;
      break;
    case InOutBlockData::Black:
      os << "Black";
      showData = false;
      break;
    default:
      os << "Undetermined";
      break;
    }
  }

  if(showData)
  {
    os << ", dataIndex: ";
    if(!iob.hasData())
      os << "<no data>";
    else
      os << iob.dataIndex();
  }

  os << "}";

  return os;
}

/**
 * \brief Verbose BlockDataType for an InOutOctree
 *
 * \note Used when generating the octree.
 */
class DynamicGrayBlockData
{
public:
  enum
  {
    NO_VERTEX = -1
  };

  using VertexIndex = axom::IndexType;
  using TriangleIndex = axom::IndexType;

  using TriangleList = std::vector<TriangleIndex>;

public:
  /**
   * \brief Default constructor for an InOutLeafData
   */
  DynamicGrayBlockData() : m_vertIndex(NO_VERTEX), m_isLeaf(true) { }

  /**
   * \brief Constructor for an InOutLeafData
   *
   * \param vInd The index of a vertex
   * (optional; default is to not set a vertex)
   */
  DynamicGrayBlockData(VertexIndex vInd, bool isLeaf)
    : m_vertIndex(vInd)
    , m_isLeaf(isLeaf)
  { }

  /**
   * \brief Copy constructor for an DynamicGrayBlockData instance
   */
  DynamicGrayBlockData(const DynamicGrayBlockData& other)
    : m_vertIndex(other.m_vertIndex)
    , m_tris(other.m_tris)
    , m_isLeaf(other.m_isLeaf)
  { }

  /**
   * \brief Assignment operator for an InOutLeafData instance
   */
  DynamicGrayBlockData& operator=(const DynamicGrayBlockData& other)
  {
    this->m_vertIndex = other.m_vertIndex;

    this->m_tris.reserve(other.m_tris.size());
    std::copy(other.m_tris.begin(),
              other.m_tris.end(),
              std::back_inserter(this->m_tris));

    this->m_isLeaf = other.m_isLeaf;

    return *this;
  }

  //        /**
  //         * \brief Removes all indexed data from this leaf
  //         */
  //        void clear()
  //        {
  //            m_isLeaf = false;
  //            m_vertIndex = NO_VERTEX;
  //            m_tris.clear();
  //            m_tris = TriangleList(0);    // reconstruct to deallocate memory
  //        }

  /**
   * \brief Equality operator to determine if two
   * DynamicGrayBlockData instances are equivalent
   */
  friend bool operator==(const DynamicGrayBlockData& lhs,
                         const DynamicGrayBlockData& rhs)
  {
    return
      //(static_cast<const BlockData&>(lhs) == static_cast<const
      // BlockData&>(rhs))
      //&&
      (lhs.m_vertIndex == rhs.m_vertIndex) &&
      (lhs.m_tris.size() == rhs.m_tris.size())  // Note: We are not
                                                // checking the contents
      // && (lhs.m_tris == rhs.m_tris)                //       of the triangle
      // array, only the size
      && lhs.m_isLeaf == rhs.m_isLeaf;
  }

public:  // Functions related to whether this is a leaf
  /** Predicate to determine if the associated block is a leaf in the octree */
  bool isLeaf() const { return m_isLeaf; }

  /** Sets a flag to determine whether the associated block is a leaf or
     internal */
  void setLeafFlag(bool isLeaf) { m_isLeaf = isLeaf; }

public:  // Functions related to the associated vertex
  /**
   * \brief Checks whether there is a vertex associated with this leaf
   */
  bool hasVertex() const { return m_vertIndex >= 0; }

  /** Sets the vertex associated with this leaf */
  void setVertex(VertexIndex vInd) { m_vertIndex = vInd; }

  /** Clears the associated vertex index */
  void clearVertex() { m_vertIndex = NO_VERTEX; }

  /** Accessor for the index of the vertex associated with this leaf */
  VertexIndex& vertexIndex() { return m_vertIndex; }

  /** Constant accessor for the index of the vertex associated with this leaf */
  const VertexIndex& vertexIndex() const { return m_vertIndex; }

public:  // Functions related to the associated triangles
  /** Check whether this Leaf has any associated triangles */
  bool hasTriangles() const { return !m_tris.empty(); }

  /**
   * Reserves space for a given number of triangles
   * \param count The number of triangles for which to reserve space
   */
  void reserveTriangles(int count) { m_tris.reserve(count); }

  /** Find the number of triangles associated with this leaf */
  int numTriangles() const { return static_cast<int>(m_tris.size()); }

  /** Associates the surface triangle with the given index with this block */
  void addTriangle(TriangleIndex tInd) { m_tris.push_back(tInd); }

  /** Returns a const reference to the list of triangle indexes associated with
     the block */
  const TriangleList& triangles() const { return m_tris; }

  /** Returns a reference to the list of triangle indexes associated with the
     block */
  TriangleList& triangles() { return m_tris; }

private:
  VertexIndex m_vertIndex;
  TriangleList m_tris;
  bool m_isLeaf;
};

/**
 * Free function to print a DynamicGrayBlockData instance to an output stream
 */
inline std::ostream& operator<<(std::ostream& os,
                                const DynamicGrayBlockData& bData)
{
  os << "DynamicGrayBlockData{";

  os << "isLeaf: " << (bData.isLeaf() ? "yes" : "no");

  os << ", vertex: ";
  if(bData.hasVertex())
    os << bData.vertexIndex();
  else
    os << "<none>";

  os << ", triangles: ";
  if(bData.hasTriangles())
  {
    int numTri = bData.numTriangles();
    os << "(" << numTri << ") {";
    for(int i = 0; i < numTri; ++i)
      os << bData.triangles()[i] << ((i == numTri - 1) ? "} " : ",");
  }

  os << "}";

  return os;
}

/**
 * \class
 * \brief Handles generation of a point containment spatial index over a surface
 * mesh
 *
 * The point containment queries determine whether a given arbitrary point in
 * space lies inside or outside of the surface.  This class depends on a
 * watertight surface mesh.  In order to repair common mesh defects, this
 * class modifies the Mesh passed to it.  Please discard all other copies
 * of the Mesh pointer.
 */
template <int DIM>
class InOutOctree : public spin::SpatialOctree<DIM, InOutBlockData>
{
private:
  friend class detail::InOutOctreeStats<DIM>;
  friend class detail::InOutOctreeValidator<DIM>;
  friend class detail::InOutOctreeMeshDumper<DIM>;

public:
  using OctreeBaseType = spin::OctreeBase<DIM, InOutBlockData>;
  using SpatialOctreeType = spin::SpatialOctree<DIM, InOutBlockData>;

  using GeometricBoundingBox = typename SpatialOctreeType::GeometricBoundingBox;
  using SpacePt = typename SpatialOctreeType::SpacePt;
  using SpaceVector = typename SpatialOctreeType::SpaceVector;
  using BlockIndex = typename SpatialOctreeType::BlockIndex;
  using GridPt = typename OctreeBaseType::GridPt;
  using SpaceRay = primal::Ray<double, DIM>;

private:
  enum GenerationState
  {
    INOUTOCTREE_UNINITIALIZED,
    INOUTOCTREE_VERTICES_INSERTED,
    INOUTOCTREE_MESH_REORDERED,
    INOUTOCTREE_ELEMENTS_INSERTED,
    INOUTOCTREE_LEAVES_COLORED
  };

  static double DEFAULT_VERTEX_WELD_THRESHOLD;
  static double DEFAULT_BOUNDING_BOX_SCALE_FACTOR;

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
    using VertexIndex = axom::IndexType;
    using TriangleIndex = axom::IndexType;
    using SurfaceMesh = mint::Mesh;

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
    using TVStride =
      slam::policies::CompileTimeStride<VertexIndex, NUM_TRI_VERTS>;
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
        fmt::format(
          "There should be a vertex in triangle {} that was not in {}",
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

public:
  using SurfaceMesh = typename MeshWrapper::SurfaceMesh;

  using VertexIndex = axom::IndexType;
  using TriangleIndex = axom::IndexType;
  using IndexRegistry = slam::FieldRegistry<slam::Set<VertexIndex>, VertexIndex>;

  using SpaceTriangle = typename MeshWrapper::SpaceTriangle;

  using MeshVertexSet = typename MeshWrapper::MeshVertexSet;
  using MeshElementSet = typename MeshWrapper::MeshElementSet;
  using TriVertIndices = typename MeshWrapper::TriVertIndices;
  using VertexIndexMap = typename MeshWrapper::VertexIndexMap;

  // Type aliases for the Relations from Gray leaf blocks to mesh entities
  static const int MAX_VERTS_PER_BLOCK = 1;
  using VertexBlockMap = slam::Map<slam::Set<>, BlockIndex>;
  using STLIndirection =
    slam::policies::STLVectorIndirection<VertexIndex, VertexIndex>;

  using GrayLeafSet = slam::PositionSet<>;
  using BVStride =
    slam::policies::CompileTimeStride<VertexIndex, MAX_VERTS_PER_BLOCK>;
  using ConstantCardinality =
    slam::policies::ConstantCardinality<VertexIndex, BVStride>;
  using GrayLeafVertexRelation = slam::StaticRelation<VertexIndex,
                                                      VertexIndex,
                                                      ConstantCardinality,
                                                      STLIndirection,
                                                      GrayLeafSet,
                                                      MeshVertexSet>;

  using VariableCardinality =
    slam::policies::VariableCardinality<VertexIndex, STLIndirection>;
  using GrayLeafElementRelation = slam::StaticRelation<VertexIndex,
                                                       VertexIndex,
                                                       VariableCardinality,
                                                       STLIndirection,
                                                       GrayLeafSet,
                                                       MeshElementSet>;
  using TriangleIndexSet = typename GrayLeafElementRelation::RelationSubset;

  using GrayLeafsLevelMap = slam::Map<slam::Set<>, GrayLeafSet>;
  using GrayLeafVertexRelationLevelMap =
    slam::Map<slam::Set<>, GrayLeafVertexRelation>;
  using GrayLeafElementRelationLevelMap =
    slam::Map<slam::Set<>, GrayLeafElementRelation>;

public:
  /**
   * \brief Construct an InOutOctree to handle containment queries on a surface
   * mesh
   *
   * \param [in] bb The spatial extent covered by the octree
   * \note We slightly scale the bounding box so all mesh elements are
   * guaranteed to be enclosed by the octree and to remedy problems we've
   * encountered related to meshes that are aligned with the octree grid
   * \note The InOutOctree modifies its mesh in an effort to repair common
   * problems.  Please make sure to discard all old copies of the meshPtr.
   */
  InOutOctree(const GeometricBoundingBox& bb, SurfaceMesh*& meshPtr)
    : SpatialOctreeType(
        GeometricBoundingBox(bb).scale(DEFAULT_BOUNDING_BOX_SCALE_FACTOR))
    , m_meshWrapper(meshPtr)
    , m_vertexToBlockMap(&m_meshWrapper.vertexSet())
    //
    , m_grayLeafsMap(&this->m_levels)
    , m_grayLeafToVertexRelationLevelMap(&this->m_levels)
    , m_grayLeafToElementRelationLevelMap(&this->m_levels)
    //
    , m_generationState(INOUTOCTREE_UNINITIALIZED)
  {
    setVertexWeldThreshold(DEFAULT_VERTEX_WELD_THRESHOLD);
  }

  /**
   * \brief Generate the spatial index over the triangle mesh
   */
  void generateIndex();

  /**
   * \brief The point containment query.
   *
   * \param pt The point at which we are checking for containment
   * \return True if the point is within (or on) the surface, false otherwise
   * \note Points outside the octree bounding box are considered outside
   */
  bool within(const SpacePt& pt) const;

  /**
   * \brief Sets the threshold for welding vertices during octree construction
   *
   * \param [in] thresh The cutoff distance at which we consider two vertices
   * to be identical during the octree construction
   *
   * \pre thresh >= 0
   * \pre This function cannot be called after the octree has been constructed
   *
   * \note The InOutOctree requires the input surface to be watertight so this
   * parameter should be set with care. A welding threshold that is too
   * high could unnecessarily merge vertices and create topological defects,
   * while a value that is too low risks leaving gaps in meshes with tolerances
   * between vertices. The default value tends to work well in practice.
   *
   * \note The code actually uses the square of the threshold for comparisons
   */
  void setVertexWeldThreshold(double thresh)
  {
    SLIC_WARNING_IF(thresh < 0.,
                    "Distance threshold for vertices cannot be negative.");

    SLIC_WARNING_IF(m_generationState > INOUTOCTREE_UNINITIALIZED,
                    "Can only set the vertex welding threshold "
                      << "before initializing the InOutOctree");

    m_vertexWeldThresholdSquared = thresh * thresh;
  }

private:
  /**
   * \brief Helper function to insert a vertex into the octree
   *
   * \param idx The index of the vertex that we are inserting
   * \param startingLevel (optional, default = 0) The octree level at which
   * to begin the search for the leaf node covering this vertex
   */
  void insertVertex(VertexIndex idx, int startingLevel = 0);

  /**
   * \brief Insert all triangles of the mesh into the octree, generating a PM
   * octree
   */
  void insertMeshTriangles();

  /**
   * \brief Set a color for each leaf block of the octree.
   *
   * Black blocks are entirely within the surface, white blocks are entirely
   * outside the surface
   * and Gray blocks intersect the surface.
   */
  void colorOctreeLeaves();

  /**
   * Use octree index over mesh vertices to convert the 'triangle soup'
   * from the stl file into an indexed triangle mesh representation.
   * In particular, all vertices in the mesh that are nearly coincident will be
   * merged,
   * and degenerate triangles (where the three vertices do not have unique
   * indices)
   * will be removed.
   */
  void updateSurfaceMeshVertices();

private:
  /**
   * \brief Checks if all indexed triangles in the block share a common vertex
   *
   * \param leafBlock [in] The current octree block
   * \param leafData [inout] The data associated with this block
   * \note A side effect of this function is that we set the leafData's vertex
   * to the common
   * vertex if one is found
   * \return True, if all triangles indexed by this leaf share a common vertex,
   * false otherwise.
   */
  bool allTrianglesIncidentInCommonVertex(const BlockIndex& leafBlock,
                                          DynamicGrayBlockData& leafData) const;

  /**
   * \brief Finds a color for the given block blk and propagates to neighbors
   *
   * \note Propagates color to same-level and coarser level neighbors
   * \param blk The block to color
   * \param blkData The data associated with this block
   * \return True if we were able to find a color for blk, false otherwise
   */
  bool colorLeafAndNeighbors(const BlockIndex& blk, InOutBlockData& blkData);

  /**
   * \brief Predicate to determine if the vertex is indexed by the blk
   *
   * \pre This function assumes the vertices have been inserted and the mesh has
   * been reordered
   * \param vIdx The index of the vertex to check
   * \param blk The block that we are checking against
   * \return true if vIdx is indexed by blk, false otherwise
   */
  bool blockIndexesVertex(VertexIndex vIdx, const BlockIndex& blk) const
  {
    SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);

    // Needs to account for non-leaf ancestors of the block
    return vIdx >= 0 && m_vertexToBlockMap[vIdx].isDescendantOf(blk);
  }

  /**
   * \brief Predicate to determine if any of the elements vertices are indexed
   * by the given BlockIndex
   *
   * \pre This function assumes the vertices have been inserted and the mesh has
   * been reordered
   * \param tIdx The index of the triangle to check
   * \param blk The block that we are checking against
   * \return true if one of the triangle's vertices are indexed by blk, false
   * otherwise
   */
  bool blockIndexesElementVertex(TriangleIndex tIdx, const BlockIndex& blk) const
  {
    SLIC_ASSERT(m_generationState >= INOUTOCTREE_MESH_REORDERED);

    TriVertIndices tVerts = m_meshWrapper.triangleVertexIndices(tIdx);
    for(int i = 0; i < tVerts.size(); ++i)
    {
      // Using the vertex-to-block cache to avoid numerical degeneracies
      if(blockIndexesVertex(tVerts[i], blk)) return true;
    }
    return false;
  }

  /**
   * \brief Determines whether the specified point is within the gray leaf
   *
   * \param queryPt The point we are querying
   * \param leafBlk The block of the gray leaf
   * \param data The data associated with the leaf block
   * \return True, if the point is inside the local surface associated with this
   * block, false otherwise
   */
  bool withinGrayBlock(const SpacePt& queryPt,
                       const BlockIndex& leafBlk,
                       const InOutBlockData& data) const;

  /**
   * \brief Returns the index of the mesh vertex associated with the given leaf
   * block
   *
   * \pre leafBlk is a leaf block of the octree
   * \param leafBlk The BlockIndex of a leaf block in the octree
   * \param leafData The data associated with this leaf block
   * \return The index of the mesh vertex associated with this leaf block
   */
  VertexIndex leafVertex(const BlockIndex& leafBlk,
                         const InOutBlockData& leafData) const;

  /**
   * \brief Returns the set of mesh triangle indices associated with the given
   * leaf block
   *
   * \pre leafBlk is a leaf block of the octree
   * \param leafBlk The BlockIndex of a leaf block in the octree
   * \param leafData The data associated with this leaf block
   * \return The set of mesh triangle indices associated with this leaf block
   */
  TriangleIndexSet leafTriangles(const BlockIndex& leafBlk,
                                 const InOutBlockData& leafData) const;

private:
  DISABLE_COPY_AND_ASSIGNMENT(InOutOctree);
  DISABLE_MOVE_AND_ASSIGNMENT(InOutOctree);

  /** \brief Checks internal consistency of the octree representation */
  void checkValid() const;

  /** \brief Helper function to verify that all leaves at the given level have a
     color */
  void checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const;

  void dumpOctreeMeshVTK(const std::string& name) const;
  void dumpTriMeshVTK(const std::string& name) const;

  /**
   * \brief Utility function to dump any Inside blocks whose neighbors are
   * outside (and vice-versa)
   *
   * \note There should not be any such blocks in a valid InOutOctree
   */
  void dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const;

  /**
   * Utility function to print some statistics about the InOutOctree instance
   */
  void printOctreeStats() const;

  /**
   * \brief Utility function to compute the angle-weighted pseudonormal for a
   * vertex in the mesh
   *
   * \note Not optimized
   * \note The returned normal is not normalized.
   */
  SpaceVector vertexNormal(VertexIndex vIdx) const
  {
    SpaceVector vec;

    BlockIndex vertexBlock = m_vertexToBlockMap[vIdx];
    TriangleIndexSet triSet = leafTriangles(vertexBlock, (*this)[vertexBlock]);
    for(int i = 0; i < triSet.size(); ++i)
    {
      TriangleIndex tIdx = triSet[i];
      TriVertIndices tv = m_meshWrapper.triangleVertexIndices(tIdx);
      if(m_meshWrapper.incidentInVertex(tv, vIdx))
      {
        int idx = (vIdx == tv[0]) ? 0 : (vIdx == tv[1] ? 1 : 2);

        SpaceTriangle tr = m_meshWrapper.trianglePositions(tIdx);
        vec += tr.angle(idx) * tr.normal();
      }
    }

    return vec;
  }

  /**
   * \brief Utility function to compute the normal for an edge of the mesh.
   *
   * The computed edge normal is the average of its face normals.
   * There should be two of these in a closed manifold surface mesh.
   * \note Not optimized
   * \note The returned normal is not normalized.
   */
  SpaceVector edgeNormal(VertexIndex vIdx1, VertexIndex vIdx2) const
  {
    SpaceVector vec;

    BlockIndex vertexBlock = m_vertexToBlockMap[vIdx1];
    TriangleIndexSet triSet = leafTriangles(vertexBlock, (*this)[vertexBlock]);
    for(int i = 0; i < triSet.size(); ++i)
    {
      TriangleIndex tIdx = triSet[i];
      TriVertIndices tv = m_meshWrapper.triangleVertexIndices(tIdx);
      if(m_meshWrapper.incidentInVertex(tv, vIdx1) &&
         m_meshWrapper.incidentInVertex(tv, vIdx2))
      {
        vec += m_meshWrapper.trianglePositions(tIdx).normal();
      }
    }

    return vec;
  }

protected:
  MeshWrapper m_meshWrapper;

  VertexBlockMap m_vertexToBlockMap;

  GrayLeafsLevelMap m_grayLeafsMap;
  GrayLeafVertexRelationLevelMap m_grayLeafToVertexRelationLevelMap;
  GrayLeafElementRelationLevelMap m_grayLeafToElementRelationLevelMap;

  GenerationState m_generationState;

  IndexRegistry m_indexRegistry;

  double m_vertexWeldThresholdSquared;

  /// Bounding box scaling factor for dealing with grazing triangles
  double m_boundingBoxScaleFactor {DEFAULT_BOUNDING_BOX_SCALE_FACTOR};
};

template <int DIM>
double InOutOctree<DIM>::DEFAULT_VERTEX_WELD_THRESHOLD = 1E-9;

template <int DIM>
double InOutOctree<DIM>::DEFAULT_BOUNDING_BOX_SCALE_FACTOR = 1.000123;

namespace
{
#ifdef AXOM_DEBUG
/**
 * \brief Utility function to print the vertex indices of a cell
 */
inline std::ostream& operator<<(std::ostream& os,
                                const InOutOctree<3>::TriVertIndices& tvInd)
{
  os << "[";
  for(int i = 0; i < tvInd.size(); ++i)
    os << tvInd[i] << ((i == tvInd.size() - 1) ? "]" : ",");

  return os;
}

inline std::ostream& operator<<(std::ostream& os,
                                const InOutOctree<3>::TriangleIndexSet& tSet)
{
  os << "[";
  for(int i = 0; i < tSet.size(); ++i)
    os << tSet[i] << ((i == tSet.size() - 1) ? "]" : ",");

  return os;
}

#endif
}  // namespace

template <int DIM>
void InOutOctree<DIM>::generateIndex()
{
  using Timer = axom::utilities::Timer;

  // Loop through mesh vertices
  SLIC_INFO(
    fmt::format("  Generating InOutOctree over surface mesh with {} vertices "
                "and {} elements.",
                m_meshWrapper.numMeshVertices(),
                m_meshWrapper.numMeshElements()));

  Timer timer;

  // STEP 1 -- Add mesh vertices to octree
  timer.start();
  int numMeshVerts = m_meshWrapper.numMeshVertices();
  for(int idx = 0; idx < numMeshVerts; ++idx)
  {
    insertVertex(idx);
  }
  timer.stop();
  m_generationState = INOUTOCTREE_VERTICES_INSERTED;
  SLIC_INFO(
    fmt::format("\t--Inserting vertices took {} seconds.", timer.elapsed()));

  // STEP 1(b) -- Update the mesh vertices and cells with after vertex welding
  // from octree
  timer.start();
  updateSurfaceMeshVertices();
  timer.stop();
  m_generationState = INOUTOCTREE_MESH_REORDERED;

  SLIC_INFO("\t--Updating mesh took " << timer.elapsed() << " seconds.");
  SLIC_INFO(
    fmt::format("  After inserting vertices, reindexed mesh has {} vertices "
                "and {} triangles.",
                m_meshWrapper.numMeshVertices(),
                m_meshWrapper.numMeshElements()));

#ifdef DUMP_OCTREE_INFO
  // -- Print some stats about the octree
  SLIC_INFO("** Octree stats after inserting vertices");
  dumpTriMeshVTK("surfaceMesh");
  dumpOctreeMeshVTK("prOctree");
  printOctreeStats();
#endif
  checkValid();

  // STEP 2 -- Add mesh triangles to octree
  timer.start();
  insertMeshTriangles();
  timer.stop();
  m_generationState = INOUTOCTREE_ELEMENTS_INSERTED;
  SLIC_INFO("\t--Inserting triangles took " << timer.elapsed() << " seconds.");

  // STEP 3 -- Color the blocks of the octree
  // -- Black (in), White(out), Gray(Intersects surface)
  timer.start();
  colorOctreeLeaves();

  timer.stop();
  m_generationState = INOUTOCTREE_LEAVES_COLORED;
  SLIC_INFO("\t--Coloring octree leaves took " << timer.elapsed() << " seconds.");

// -- Print some stats about the octree
#ifdef DUMP_OCTREE_INFO
  SLIC_INFO("** Octree stats after inserting triangles");
  dumpOctreeMeshVTK("pmOctree");
  dumpDifferentColoredNeighborsMeshVTK("differentNeighbors");
  printOctreeStats();
#endif
  checkValid();

  // CLEANUP -- Finally, fix up the surface mesh after octree operations
  timer.start();
  m_meshWrapper.regenerateSurfaceMesh();
  timer.stop();
  SLIC_INFO("\t--Regenerating the mesh took " << timer.elapsed() << " seconds.");

  SLIC_INFO("  Finished generating the InOutOctree.");
}

template <int DIM>
void InOutOctree<DIM>::insertVertex(VertexIndex idx, int startingLevel)
{
  const SpacePt pt = m_meshWrapper.getMeshVertexPosition(idx);

  BlockIndex block = this->findLeafBlock(pt, startingLevel);
  InOutBlockData& blkData = (*this)[block];

  QUEST_OCTREE_DEBUG_LOG_IF(
    idx == DEBUG_VERT_IDX,
    fmt::format("\t -- inserting pt {} with index {}. "
                "Looking at block {} w/ blockBB {} indexing leaf vertex {}",
                pt,
                idx,
                block,
                this->blockBoundingBox(block),
                blkData.dataIndex()));

  if(!blkData.hasData())
  {
    blkData.setData(idx);

    // Update the vertex-to-block map for this vertex
    if(m_generationState >= INOUTOCTREE_MESH_REORDERED)
      m_vertexToBlockMap[idx] = block;
  }
  else
  {
    // check if we should merge the vertices
    VertexIndex origVertInd = blkData.dataIndex();
    if(squared_distance(pt, m_meshWrapper.getMeshVertexPosition(origVertInd)) >=
       m_vertexWeldThresholdSquared)
    {
      blkData.clear();
      this->refineLeaf(block);

      insertVertex(origVertInd, block.childLevel());
      insertVertex(idx, block.childLevel());
    }
  }

  QUEST_OCTREE_DEBUG_LOG_IF(
    blkData.dataIndex() == DEBUG_VERT_IDX,
    fmt::format("-- vertex {} is indexed in block {}. Leaf vertex is {}",
                idx,
                block,
                blkData.dataIndex()));
}

template <int DIM>
void InOutOctree<DIM>::insertMeshTriangles()
{
  using Timer = axom::utilities::Timer;
  using LeavesLevelMap = typename OctreeBaseType::OctreeLevelType;

  SLIC_ASSERT(m_meshWrapper.meshWasReindexed());

  // Temporary arrays of DyamicGrayBlockData for current and next level
  using DynamicLevelData = std::vector<DynamicGrayBlockData>;
  const int NUM_INIT_DATA_ENTRIES = 1 << 10;
  DynamicLevelData currentLevelData;
  DynamicLevelData nextLevelData;
  currentLevelData.reserve(NUM_INIT_DATA_ENTRIES);
  nextLevelData.reserve(NUM_INIT_DATA_ENTRIES);

  /// --- Initialize root level data
  BlockIndex rootBlock = this->root();
  InOutBlockData& rootData = (*this)[rootBlock];

  currentLevelData.push_back(DynamicGrayBlockData());
  DynamicGrayBlockData& dynamicRootData = currentLevelData[0];
  if(rootData.hasData()) dynamicRootData.setVertex(rootData.dataIndex());
  dynamicRootData.setLeafFlag(rootData.isLeaf());
  rootData.setData(0);

  // Add all triangle references to the root
  int const numTris = m_meshWrapper.numMeshElements();
  dynamicRootData.triangles().reserve(numTris);
  for(int idx = 0; idx < numTris; ++idx)
  {
    dynamicRootData.addTriangle(idx);
  }

  // Iterate through octree levels
  // and insert triangles into the blocks that they intersect
  for(int lev = 0; lev < this->m_levels.size(); ++lev)
  {
    Timer levelTimer(true);

    auto& gvRelData = m_indexRegistry.addNamelessBuffer();
    auto& geIndRelData = m_indexRegistry.addNamelessBuffer();
    auto& geSizeRelData = m_indexRegistry.addNamelessBuffer();
    geSizeRelData.push_back(0);

    int nextLevelDataBlockCounter = 0;

    auto& levelLeafMap = this->getOctreeLevel(lev);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      InOutBlockData& blkData = *it;

      if(!blkData.hasData()) continue;

      BlockIndex blk(it.pt(), lev);
      DynamicGrayBlockData& dynamicLeafData =
        currentLevelData[blkData.dataIndex()];

      bool isInternal = !dynamicLeafData.isLeaf();
      bool isLeafThatMustRefine = !isInternal &&
        !allTrianglesIncidentInCommonVertex(blk, dynamicLeafData);

      QUEST_OCTREE_DEBUG_LOG_IF(
        DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk,
        fmt::format("Attempting to insert triangles from block {}."
                    "\n\tDynamic data: {}"
                    "\n\tBlock data: {}"
                    "\n\tAbout to finalize? {}",
                    blk,
                    dynamicLeafData,
                    blkData,
                    (!isInternal && !isLeafThatMustRefine ? " yes" : "no")));

      // Leaf blocks that don't refine are 'finalized'
      // -- add  them to the current level's relations
      if(!isInternal && !isLeafThatMustRefine)
      {
        if(dynamicLeafData.hasTriangles())
        {
          // Set the leaf data in the octree
          blkData.setData(static_cast<int>(gvRelData.size()));

          // Add the vertex index to the gray blocks vertex relation
          gvRelData.push_back(dynamicLeafData.vertexIndex());

          // Add the triangles to the gray block's element relations
          std::copy(dynamicLeafData.triangles().begin(),
                    dynamicLeafData.triangles().end(),
                    std::back_inserter(geIndRelData));
          geSizeRelData.push_back(static_cast<int>(geIndRelData.size()));

          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == blk || DEBUG_BLOCK_2 == blk,
            fmt::format("[Added block {} into tree as a gray leaf]."
                        "\n\tDynamic data: {}"
                        "\n\tBlock data: {}",
                        blk,
                        dynamicLeafData,
                        blkData));
        }
      }
      else
      {
        /// Otherwise, we must distribute the block data among the children

        // Refine the leaf if necessary
        if(isLeafThatMustRefine)
        {
          const VertexIndex vIdx = dynamicLeafData.vertexIndex();

          this->refineLeaf(blk);
          dynamicLeafData.setLeafFlag(false);

          // Reinsert the vertex into the tree, if vIdx was indexed by blk
          if(blockIndexesVertex(vIdx, blk))
            insertVertex(vIdx, blk.childLevel());
        }
        else if(isInternal)
        {
          // Need to mark the leaf as internal since we were using its data
          // as an index into the DynamicGrayBlockData array
          blkData.setInternal();
        }

        SLIC_ASSERT_MSG(
          this->isInternal(blk),
          fmt::format(
            "Block {} was refined, so it should be marked as internal.",
            blk));

        /// Setup caches for data associated with children
        BlockIndex childBlk[BlockIndex::NUM_CHILDREN];
        GeometricBoundingBox childBB[BlockIndex::NUM_CHILDREN];
        DynamicGrayBlockData childData[BlockIndex::NUM_CHILDREN];
        DynamicGrayBlockData* childDataPtr[BlockIndex::NUM_CHILDREN];

        const typename LeavesLevelMap::BroodData& broodData =
          this->getOctreeLevel(lev + 1).getBroodData(blk.pt());

        for(int j = 0; j < BlockIndex::NUM_CHILDREN; ++j)
        {
          childBlk[j] = blk.child(j);
          childBB[j] = this->blockBoundingBox(childBlk[j]);

          // expand bounding box slightly to deal with grazing triangles
          childBB[j].scale(m_boundingBoxScaleFactor);

          const InOutBlockData& childBlockData = broodData[j];
          if(!childBlockData.hasData())
          {
            childData[j] = DynamicGrayBlockData();
            childData[j].setLeafFlag(childBlockData.isLeaf());
          }
          else
          {
            childData[j] = DynamicGrayBlockData(childBlockData.dataIndex(),
                                                childBlockData.isLeaf());
          }

          childDataPtr[j] = &childData[j];
        }

        // Check that the vector has enough capacity for all eight children
        // This ensures that our child data pointers will not be invalidated
        if(nextLevelData.capacity() <
           (nextLevelData.size() + BlockIndex::NUM_CHILDREN))
          nextLevelData.reserve(nextLevelData.size() * 4);

        // Add all triangles to intersecting children blocks
        DynamicGrayBlockData::TriangleList& parentTriangles =
          dynamicLeafData.triangles();
        int numTriangles = static_cast<int>(parentTriangles.size());
        for(int i = 0; i < numTriangles; ++i)
        {
          TriangleIndex tIdx = parentTriangles[i];
          SpaceTriangle spaceTri = m_meshWrapper.trianglePositions(tIdx);
          GeometricBoundingBox tBB = m_meshWrapper.triangleBoundingBox(tIdx);

          for(int j = 0; j < BlockIndex::numChildren(); ++j)
          {
            bool shouldAddTriangle =
              blockIndexesElementVertex(tIdx, childBlk[j]) ||
              (childDataPtr[j]->isLeaf() ? intersect(spaceTri, childBB[j])
                                         : intersect(tBB, childBB[j]));

            QUEST_OCTREE_DEBUG_LOG_IF(
              DEBUG_BLOCK_1 == childBlk[j] || DEBUG_BLOCK_2 == childBlk[j],
              //&& tIdx == DEBUG_TRI_IDX
              fmt::format("Attempting to insert triangle {} @ {} w/ BB {}"
                          "\n\t into block {} w/ BB {} and data {} "
                          "\n\tShould add? {}",
                          tIdx,
                          spaceTri,
                          tBB,
                          childBlk[j],
                          childBB[j],
                          *childDataPtr[j],
                          (shouldAddTriangle ? " yes" : "no")));

            if(shouldAddTriangle)
            {
              // Place the DynamicGrayBlockData in the array before adding its data
              if(!childDataPtr[j]->hasTriangles())
              {
                // Copy the DynamicGrayBlockData into the array
                nextLevelData.push_back(childData[j]);

                // Update the child data pointer
                childDataPtr[j] = &nextLevelData[nextLevelDataBlockCounter];

                // Set the data in the octree to this index and update the index
                (*this)[childBlk[j]].setData(nextLevelDataBlockCounter++);
              }

              childDataPtr[j]->addTriangle(tIdx);

              QUEST_OCTREE_DEBUG_LOG_IF(
                DEBUG_BLOCK_1 == childBlk[j] || DEBUG_BLOCK_2 == childBlk[j],
                //&& tIdx == DEBUG_TRI_IDX
                fmt::format("Added triangle {} @ {} with verts {}"
                            "\n\tinto block {} with data {}.",
                            tIdx,
                            spaceTri,
                            m_meshWrapper.triangleVertexIndices(tIdx),
                            childBlk[j],
                            *(childDataPtr[j])));
            }
          }
        }
      }
    }

    if(!levelLeafMap.empty())
    {
      // Create the relations from gray leaves to mesh vertices and elements
      m_grayLeafsMap[lev] = GrayLeafSet(static_cast<int>(gvRelData.size()));

      m_grayLeafToVertexRelationLevelMap[lev] =
        GrayLeafVertexRelation(&m_grayLeafsMap[lev], &m_meshWrapper.vertexSet());
      m_grayLeafToVertexRelationLevelMap[lev].bindIndices(
        static_cast<int>(gvRelData.size()),
        &gvRelData);

      m_grayLeafToElementRelationLevelMap[lev] =
        GrayLeafElementRelation(&m_grayLeafsMap[lev],
                                &m_meshWrapper.elementSet());
      m_grayLeafToElementRelationLevelMap[lev].bindBeginOffsets(
        m_grayLeafsMap[lev].size(),
        &geSizeRelData);
      m_grayLeafToElementRelationLevelMap[lev].bindIndices(
        static_cast<int>(geIndRelData.size()),
        &geIndRelData);
    }

    currentLevelData.clear();
    nextLevelData.swap(currentLevelData);

    if(!levelLeafMap.empty())
      SLIC_DEBUG("\tInserting triangles into level "
                 << lev << " took " << levelTimer.elapsed() << " seconds.");
  }
}

template <int DIM>
void InOutOctree<DIM>::colorOctreeLeaves()
{
  // Note (KW): Existence of leaf implies that either
  // * it is gray
  // * one of its siblings is gray
  // * one of its siblings has a gray descendant

  using Timer = axom::utilities::Timer;
  using GridPtVec = std::vector<GridPt>;
  GridPtVec uncoloredBlocks;

  // Bottom-up traversal of octree
  for(int lev = this->maxLeafLevel() - 1; lev >= 0; --lev)
  {
    uncoloredBlocks.clear();
    Timer levelTimer(true);

    auto& levelLeafMap = this->getOctreeLevel(lev);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      if(!it->isLeaf()) continue;

      BlockIndex leafBlk(it.pt(), lev);
      InOutBlockData& blockData = *it;
      if(!colorLeafAndNeighbors(leafBlk, blockData))
        uncoloredBlocks.push_back(leafBlk.pt());
    }

    // Iterate through the uncolored blocks until all have a color
    // This terminates since we know that one of its siblings
    // (or their descendants) is gray
    while(!uncoloredBlocks.empty())
    {
      int prevCount = static_cast<int>(uncoloredBlocks.size());
      AXOM_DEBUG_VAR(prevCount);

      GridPtVec prevVec;
      prevVec.swap(uncoloredBlocks);
      auto end = prevVec.end();
      for(auto it = prevVec.begin(); it < end; ++it)
      {
        BlockIndex leafBlk(*it, lev);
        if(!colorLeafAndNeighbors(leafBlk, (*this)[leafBlk]))
          uncoloredBlocks.push_back(*it);
      }

      SLIC_ASSERT_MSG(
        static_cast<int>(uncoloredBlocks.size()) < prevCount,
        fmt::format("Problem coloring leaf blocks at level {}. "
                    "There are {} blocks that are still not colored. "
                    "First problem block is: {}",
                    lev,
                    uncoloredBlocks.size(),
                    BlockIndex(uncoloredBlocks[0], lev)));
    }

    if(!levelLeafMap.empty())
    {
      checkAllLeavesColoredAtLevel(lev);
      SLIC_DEBUG(fmt::format("\tColoring level {} took {} seconds.",
                             lev,
                             levelTimer.elapsed()));
    }
  }
}

template <int DIM>
bool InOutOctree<DIM>::colorLeafAndNeighbors(const BlockIndex& leafBlk,
                                             InOutBlockData& leafData)
{
  bool isColored = leafData.isColored();

  QUEST_OCTREE_DEBUG_LOG_IF(
    leafBlk == DEBUG_BLOCK_1 || leafBlk == DEBUG_BLOCK_2,
    fmt::format("Trying to color {} with data: {}", leafBlk, leafData));

  if(!isColored)
  {
    // Leaf does not yet have a color... try to find its color from same-level
    // face neighbors
    for(int i = 0; !isColored && i < leafBlk.numFaceNeighbors(); ++i)
    {
      BlockIndex neighborBlk = leafBlk.faceNeighbor(i);
      if(this->isLeaf(neighborBlk))
      {
        const InOutBlockData& neighborData = (*this)[neighborBlk];

        QUEST_OCTREE_DEBUG_LOG_IF(
          DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk ||
            DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
          fmt::format("Spreading color to block {} with data {}, "
                      "bounding box {} w/ midpoint {}"
                      "\n\t\t from block {} with data {}, "
                      "bounding box {} w/ midpoint {}.",
                      leafBlk,
                      leafData,
                      this->blockBoundingBox(leafBlk),
                      this->blockBoundingBox(leafBlk).getCentroid(),
                      neighborBlk,
                      neighborData,
                      this->blockBoundingBox(neighborBlk),
                      this->blockBoundingBox(neighborBlk).getCentroid()));

        switch(neighborData.color())
        {
        case InOutBlockData::Black:
          leafData.setBlack();
          break;
        case InOutBlockData::White:
          leafData.setWhite();
          break;
        case InOutBlockData::Gray:
        {
          SpacePt faceCenter =
            SpacePt::midpoint(this->blockBoundingBox(leafBlk).getCentroid(),
                              this->blockBoundingBox(neighborBlk).getCentroid());
          if(withinGrayBlock(faceCenter, neighborBlk, neighborData))
            leafData.setBlack();
          else
            leafData.setWhite();
        }
        break;
        case InOutBlockData::Undetermined:
          break;
        }

        isColored = leafData.isColored();

        QUEST_OCTREE_DEBUG_LOG_IF(
          isColored &&
            (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk),
          fmt::format("Leaf block was colored -- {} now has data {}",
                      leafBlk,
                      leafData));
      }
    }
  }

  // If the block has a color, try to color its face neighbors at the same or
  // coarser resolution
  if(isColored)
  {
    for(int i = 0; i < leafBlk.numFaceNeighbors(); ++i)
    {
      BlockIndex neighborBlk = this->coveringLeafBlock(leafBlk.faceNeighbor(i));
      if(neighborBlk != BlockIndex::invalid_index())
      {
        InOutBlockData& neighborData = (*this)[neighborBlk];
        if(!neighborData.isColored())
        {
          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk ||
              DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
            fmt::format("Spreading color from block {} with data {}, "
                        "bounding box {} w/ midpoint {}"
                        "\n\t\t to block {} with data {}, "
                        "bounding box {} w/ midpoint {}.",
                        leafBlk,
                        leafData,
                        this->blockBoundingBox(leafBlk),
                        this->blockBoundingBox(leafBlk).getCentroid(),
                        neighborBlk,
                        neighborData,
                        this->blockBoundingBox(neighborBlk),
                        this->blockBoundingBox(neighborBlk).getCentroid()));

          switch(leafData.color())
          {
          case InOutBlockData::Black:
            neighborData.setBlack();
            break;
          case InOutBlockData::White:
            neighborData.setWhite();
            break;
          case InOutBlockData::Gray:
          {
            // Check the center of the shared face between the same-level neighbor of the gray block
            SpacePt faceCenter = SpacePt::midpoint(
              this->blockBoundingBox(leafBlk).getCentroid(),
              this->blockBoundingBox(leafBlk.faceNeighbor(i)).getCentroid());

            if(withinGrayBlock(faceCenter, leafBlk, leafData))
              neighborData.setBlack();
            else
              neighborData.setWhite();
          }
          break;
          case InOutBlockData::Undetermined:
            break;
          }

          QUEST_OCTREE_DEBUG_LOG_IF(
            neighborData.isColored() &&
              (DEBUG_BLOCK_1 == neighborBlk || DEBUG_BLOCK_2 == neighborBlk),
            fmt::format("Neighbor block was colored -- {} now has data {}",
                        neighborBlk,
                        neighborData));
        }
      }
    }
  }

  return isColored;
}

template <int DIM>
typename InOutOctree<DIM>::VertexIndex InOutOctree<DIM>::leafVertex(
  const BlockIndex& leafBlk,
  const InOutBlockData& leafData) const
{
  if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
  {
    SLIC_ASSERT(leafData.hasData());
    return m_grayLeafToVertexRelationLevelMap[leafBlk.level()]
                                             [leafData.dataIndex()][0];
  }
  else
  {
    return leafData.dataIndex();
  }
}

template <int DIM>
typename InOutOctree<DIM>::TriangleIndexSet InOutOctree<DIM>::leafTriangles(
  const BlockIndex& leafBlk,
  const InOutBlockData& leafData) const
{
  SLIC_ASSERT(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED &&
              leafData.hasData());

  return m_grayLeafToElementRelationLevelMap[leafBlk.level()][leafData.dataIndex()];
}

template <int DIM>
bool InOutOctree<DIM>::withinGrayBlock(const SpacePt& queryPt,
                                       const BlockIndex& leafBlk,
                                       const InOutBlockData& leafData) const
{
  /// Finds a ray from queryPt to a point of a triangle within leafBlk.
  /// Then find the first triangle along this ray. The orientation of the ray
  /// against this triangle's normal indicates queryPt's containment.
  /// It is inside when the dot product is positive.

  SLIC_ASSERT(leafData.color() == InOutBlockData::Gray);
  SLIC_ASSERT(leafData.hasData());

  GeometricBoundingBox blockBB = this->blockBoundingBox(leafBlk);

  SpacePt triPt;

  TriangleIndexSet triSet = leafTriangles(leafBlk, leafData);
  const int numTris = triSet.size();
  for(int i = 0; i < numTris; ++i)
  {
    /// Get the triangle
    TriangleIndex idx = triSet[i];
    SpaceTriangle tri = m_meshWrapper.trianglePositions(idx);

    /// Find a point from this triangle within the bounding box of the mesh
    primal::Polygon<double, DIM> poly = primal::clip(tri, blockBB);
    if(poly.numVertices() == 0)
    {
      // Account for cases where the triangle only grazes the bounding box.
      // Here, intersect(tri,blockBB) is true, but the clipping algorithm
      // produces an empty polygon.  To resolve this, clip against a
      // slightly expanded bounding box
      GeometricBoundingBox expandedBB = blockBB;
      expandedBB.scale(m_boundingBoxScaleFactor);

      poly = primal::clip(tri, expandedBB);

      // If that still doesn't work, move on to the next triangle
      if(poly.numVertices() == 0)
      {
        continue;
      }
    }

    triPt = poly.centroid();

    /// Use a ray from the query point to the triangle point to find an
    /// intersection. Note: We have to check all triangles to ensure that
    /// there is not a closer triangle than tri along this direction.
    TriangleIndex tIdx = MeshWrapper::NO_TRIANGLE;
    double minRayParam = std::numeric_limits<double>::infinity();
    SpaceRay ray(queryPt, SpaceVector(queryPt, triPt));

    QUEST_OCTREE_DEBUG_LOG_IF(
      DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
      fmt::format("Checking if pt {} is within block {} with data {}, "
                  "ray is {}, triangle point is {} on triangle with index {}.",
                  queryPt,
                  leafBlk,
                  leafData,
                  ray,
                  triPt,
                  idx));

    double rayParam = 0;
    if(primal::intersect(tri, ray, rayParam))
    {
      minRayParam = rayParam;
      tIdx = idx;

      QUEST_OCTREE_DEBUG_LOG_IF(
        DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
        fmt::format("... intersection for triangle w/ index {} at ray "
                    "parameter {} at point {}",
                    tIdx,
                    minRayParam,
                    ray.at(minRayParam)));
    }

    for(int j = 0; j < numTris; ++j)
    {
      TriangleIndex localIdx = triSet[j];
      if(localIdx == idx) continue;

      if(primal::intersect(m_meshWrapper.trianglePositions(localIdx), ray, rayParam))
      {
        if(rayParam < minRayParam)
        {
          minRayParam = rayParam;
          tIdx = localIdx;

          QUEST_OCTREE_DEBUG_LOG_IF(
            DEBUG_BLOCK_1 == leafBlk || DEBUG_BLOCK_2 == leafBlk,
            fmt::format("... intersection for triangle w/ index {} at ray "
                        "parameter {} at point {}",
                        tIdx,
                        minRayParam,
                        ray.at(minRayParam)));
        }
      }
    }

    if(tIdx == MeshWrapper::NO_TRIANGLE)
    {
      continue;
    }

    // Inside when the dot product of the normal with this triangle is positive
    SpaceVector normal = (tIdx == idx)
      ? tri.normal()
      : m_meshWrapper.trianglePositions(tIdx).normal();

    return normal.dot(ray.direction()) > 0.;
  }

  SLIC_DEBUG("Could not determine inside/outside for point "
             << queryPt << " on block " << leafBlk);

  return false;  // query points on boundary might get here -- revisit this.
}

template <int DIM>
void InOutOctree<DIM>::updateSurfaceMeshVertices()
{
  // Create a map from old vertex indices to new vertex indices
  MeshVertexSet origVerts(m_meshWrapper.numMeshVertices());
  VertexIndexMap vertexIndexMap(&origVerts, MeshWrapper::NO_VERTEX);

  // Generate unique indices for new mesh vertices
  int uniqueVertexCounter = 0;
  for(int i = 0; i < origVerts.size(); ++i)
  {
    // Find the block and its indexed vertex in the octree
    BlockIndex leafBlock =
      this->findLeafBlock(m_meshWrapper.getMeshVertexPosition(i));
    SLIC_ASSERT((*this)[leafBlock].hasData());
    VertexIndex vInd = (*this)[leafBlock].dataIndex();

    // If the indexed vertex doesn't have a new id, give it one
    if(vertexIndexMap[vInd] == MeshWrapper::NO_VERTEX)
      vertexIndexMap[vInd] = uniqueVertexCounter++;

    // If this is not the indexed vertex of the block, set the new index
    if(vInd != i) vertexIndexMap[i] = vertexIndexMap[vInd];
  }

  // Use the index map to reindex the mesh verts and elements
  m_meshWrapper.reindexMesh(uniqueVertexCounter, vertexIndexMap);

  // Update the octree leaf vertex IDs to the new mesh IDs
  // and create the map from the new vertices to their octree blocks
  m_vertexToBlockMap = VertexBlockMap(&m_meshWrapper.vertexSet());
  for(int i = 0; i < m_meshWrapper.numMeshVertices(); ++i)
  {
    const SpacePt& pos = m_meshWrapper.vertexPosition(i);
    BlockIndex leafBlock = this->findLeafBlock(pos);
    SLIC_ASSERT(this->isLeaf(leafBlock) && (*this)[leafBlock].hasData());

    (*this)[leafBlock].setData(i);
    m_vertexToBlockMap[i] = leafBlock;
  }
}

template <int DIM>
bool InOutOctree<DIM>::allTrianglesIncidentInCommonVertex(
  const BlockIndex& leafBlock,
  DynamicGrayBlockData& leafData) const
{
  bool shareCommonVert = false;

  VertexIndex commonVert = leafData.vertexIndex();

  const int numTris = leafData.numTriangles();
  const auto& tris = leafData.triangles();

  if(blockIndexesVertex(commonVert, leafBlock))
  {
    // This is a leaf node containing the indexed vertex
    // Loop through the triangles and check that all are incident with this
    // vertex
    for(int i = 0; i < numTris; ++i)
    {
      if(!m_meshWrapper.incidentInVertex(
           m_meshWrapper.triangleVertexIndices(tris[i]),
           commonVert))
      {
        return false;
      }
    }
    shareCommonVert = true;
  }
  else
  {
    SLIC_ASSERT(numTris > 0);
    switch(numTris)
    {
    case 1:
      /// Choose an arbitrary vertex from this triangle
      commonVert = m_meshWrapper.triangleVertexIndices(tris[0])[0];
      shareCommonVert = true;
      break;
    case 2:
      /// Find a vertex that both triangles share
      shareCommonVert =
        m_meshWrapper.haveSharedVertex(tris[0], tris[1], commonVert);
      break;
    default:  // numTris >= 3
      /// Find a vertex that the first three triangles share
      shareCommonVert =
        m_meshWrapper.haveSharedVertex(tris[0], tris[1], tris[2], commonVert);

      /// Check that all other triangles have this vertex
      for(int i = 3; shareCommonVert && i < numTris; ++i)
      {
        if(!m_meshWrapper.incidentInVertex(
             m_meshWrapper.triangleVertexIndices(tris[i]),
             commonVert))
          shareCommonVert = false;
      }
      break;
    }

    if(shareCommonVert) leafData.setVertex(commonVert);
  }

  return shareCommonVert;
}

template <int DIM>
bool InOutOctree<DIM>::within(const SpacePt& pt) const
{
  if(this->boundingBox().contains(pt))
  {
    const BlockIndex block = this->findLeafBlock(pt);
    const InOutBlockData& data = (*this)[block];

    switch(data.color())
    {
    case InOutBlockData::Black:
      return true;
    case InOutBlockData::White:
      return false;
    case InOutBlockData::Gray:
      return withinGrayBlock(pt, block, data);
    case InOutBlockData::Undetermined:
      SLIC_ASSERT_MSG(false,
                      "Error -- All leaf blocks must have a color. "
                        << " The color of leafBlock " << block
                        << " was 'Undetermined' when querying point " << pt);
      break;
    }
  }

  return false;
}

template <int DIM>
void InOutOctree<DIM>::printOctreeStats() const
{
  detail::InOutOctreeStats<DIM> octreeStats(*this);
  SLIC_INFO(octreeStats.summaryStats());

#ifdef DUMP_VTK_MESH
  // Print out some debug meshes for vertex, triangle and/or blocks defined in
  // DEBUG_XXX macros
  if(m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED)
  {
    detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);

    if(DEBUG_VERT_IDX >= 0 && DEBUG_VERT_IDX < m_meshWrapper.numMeshVertices())
    {
      meshDumper.dumpLocalOctreeMeshesForVertex("debug_", DEBUG_VERT_IDX);
    }
    if(DEBUG_TRI_IDX >= 0 && DEBUG_TRI_IDX < m_meshWrapper.numMeshElements())
    {
      meshDumper.dumpLocalOctreeMeshesForTriangle("debug_", DEBUG_TRI_IDX);
    }

    if(DEBUG_BLOCK_1 != BlockIndex::invalid_index() &&
       this->hasBlock(DEBUG_BLOCK_1))
    {
      meshDumper.dumpLocalOctreeMeshesForBlock("debug_", DEBUG_BLOCK_1);
    }

    if(DEBUG_BLOCK_2 != BlockIndex::invalid_index() &&
       this->hasBlock(DEBUG_BLOCK_2))
    {
      meshDumper.dumpLocalOctreeMeshesForBlock("debug_", DEBUG_BLOCK_2);
    }
  }
#endif
}

template <int DIM>
void InOutOctree<DIM>::checkAllLeavesColoredAtLevel(int AXOM_DEBUG_PARAM(level)) const
{
#ifdef AXOM_DEBUG
  detail::InOutOctreeValidator<DIM> validator(*this);
  validator.checkAllLeavesColoredAtLevel(level);
#endif
}

template <int DIM>
void InOutOctree<DIM>::checkValid() const
{
#ifdef AXOM_DEBUG
  SLIC_DEBUG("Inside InOutOctree::checkValid() to verify state of "
             << (m_generationState >= INOUTOCTREE_ELEMENTS_INSERTED ? "PM" : "PR")
             << " octree.");

  detail::InOutOctreeValidator<DIM> validator(*this);
  validator.checkValid();

  SLIC_DEBUG("done.");
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpTriMeshVTK(const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpTriMeshVTK(name);

#else
  AXOM_DEBUG_VAR(name);  // avoids warning about unsued param
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpOctreeMeshVTK(const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpOctreeMeshVTK(name);

#else
  AXOM_DEBUG_VAR(name);  // avoids warning about unsued param
#endif
}

template <int DIM>
void InOutOctree<DIM>::dumpDifferentColoredNeighborsMeshVTK(
  const std::string& name) const
{
#ifdef DUMP_VTK_MESH

  detail::InOutOctreeMeshDumper<DIM> meshDumper(*this);
  meshDumper.dumpDifferentColoredNeighborsMeshVTK(name);

#else
  AXOM_DEBUG_VAR(name);  // avoids warning about unsued param
#endif
}

namespace detail
{
template <int DIM>
class InOutOctreeMeshDumper
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using TriangleIndexSet = typename InOutOctreeType::TriangleIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using SpacePt = typename InOutOctreeType::SpacePt;
  using GridPt = typename InOutOctreeType::GridPt;
  using VertexIndex = typename InOutOctreeType::VertexIndex;
  using TriangleIndex = typename InOutOctreeType::TriangleIndex;
  using TriVertIndices = typename InOutOctreeType::MeshWrapper::TriVertIndices;
  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpaceTriangle = typename InOutOctreeType::SpaceTriangle;

  using LeafVertMap = slam::Map<slam::Set<VertexIndex>, VertexIndex>;
  using LeafIntMap = slam::Map<slam::Set<axom::IndexType>, axom::IndexType>;
  using LeafGridPtMap = slam::Map<slam::Set<axom::IndexType>, GridPt>;

  using DebugMesh = mint::UnstructuredMesh<mint::MIXED_SHAPE>;

  using ColorsMap = std::map<InOutBlockData::LeafColor, int>;

  using GridPtHash = spin::PointHash<typename GridPt::CoordType>;
  using GridIntMap = std::unordered_map<GridPt, int, GridPtHash>;

public:
  InOutOctreeMeshDumper(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
  {
    // Create a small lookup table to map block colors to ints
    m_colorsMap[InOutBlockData::White] = -1;
    m_colorsMap[InOutBlockData::Gray] = 0;
    m_colorsMap[InOutBlockData::Black] = 1;
  }

  /**
   *  Generates a hexahedral VTK mesh with all neighboring blocks where one is
   * inside and the other is outside
   *  \note By construction, there should be no such pairs in a valid
   * InOutOctree mesh.
   */
  void dumpDifferentColoredNeighborsMeshVTK(const std::string& name) const
  {
    if(m_generationState < InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      SLIC_INFO("Need to generate octree colors before visualizing them.");
      return;
    }

    using LevelGridIntMap = slam::Map<slam::Set<>, GridIntMap>;
    LevelGridIntMap diffBlocks(&(m_octree.m_levels));

    int totalBlocks = 0;

    // Iterate through the octree leaves
    // looking for neighbor blocks with different labelings
    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      diffBlocks[lev] = GridIntMap(0, GridPtHash());

      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;
        if(data.isLeaf() && data.color() != InOutBlockData::Gray)
        {
          for(int i = 0; i < block.numFaceNeighbors(); ++i)
          {
            BlockIndex neighborBlk =
              m_octree.coveringLeafBlock(block.faceNeighbor(i));
            if(neighborBlk != BlockIndex::invalid_index())
            {
              const InOutBlockData& neighborData = m_octree[neighborBlk];
              switch(neighborData.color())
              {
              case InOutBlockData::Black:  // intentional fallthrough
              case InOutBlockData::White:
                if(data.color() != neighborData.color())
                {
                  diffBlocks[lev][block.pt()] = 1;
                  diffBlocks[neighborBlk.level()][neighborBlk.pt()] = 1;
                }
                break;
              case InOutBlockData::Gray:
              case InOutBlockData::Undetermined:
                break;
              }
            }
          }
        }
      }
      totalBlocks += diffBlocks[lev].size();
    }

    // Add all such blocks to a vector
    std::vector<BlockIndex> blocks;
    blocks.reserve(totalBlocks);

    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      for(auto&& blk : diffBlocks[lev])
      {
        blocks.emplace_back(BlockIndex(blk.first, lev));
      }
    }

    // Generate a VTK mesh with these blocks
    dumpOctreeMeshBlocks(name, blocks, false);
  }

  void dumpLocalOctreeMeshesForVertex(const std::string& name,
                                      VertexIndex vIdx) const
  {
    std::stringstream sstr;
    sstr << name << "vertex_" << vIdx;

    BlockIndex vertexBlock = m_octree.m_vertexToBlockMap[vIdx];

    // Dump a mesh for the vertex's containing block
    std::stringstream blockStr;
    blockStr << sstr.str() << "_block_" << vertexBlock.pt()[0] << "_"
             << vertexBlock.pt()[1] << "_" << vertexBlock.pt()[2];
    std::vector<BlockIndex> blocks;
    blocks.push_back(vertexBlock);
    dumpOctreeMeshBlocks(blockStr.str(), blocks, true);

    // Dump a mesh for the incident triangles
    std::stringstream triStr;
    triStr << sstr.str() << "_triangles";
    std::vector<TriangleIndex> tris;

    TriangleIndexSet triSet =
      m_octree.leafTriangles(vertexBlock, m_octree[vertexBlock]);
    for(int i = 0; i < triSet.size(); ++i)
    {
      TriangleIndex tIdx = triSet[i];
      TriVertIndices tv = m_octree.m_meshWrapper.triangleVertexIndices(tIdx);
      for(int j = 0; j < tv.size(); ++j)
      {
        if(tv[j] == vIdx) tris.push_back(tIdx);
      }
    }
    dumpTriangleMesh(triStr.str(), tris, true);
  }

  void dumpLocalOctreeMeshesForTriangle(const std::string& name,
                                        TriangleIndex tIdx) const
  {
    std::stringstream sstr;
    sstr << name << "triangle_" << tIdx;

    // Dump a triangle mesh with the single triangle
    std::vector<TriangleIndex> tris;
    tris.push_back(tIdx);
    dumpTriangleMesh(sstr.str(), tris, true);

    // Dump a hex mesh of all blocks that index this triangle
    std::vector<BlockIndex> blocks;
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        if(it->isLeaf() && it->hasData())
        {
          BlockIndex leafblk(it.pt(), lev);
          TriangleIndexSet triSet = m_octree.leafTriangles(leafblk, *it);

          bool found = false;
          for(int i = 0; !found && i < triSet.size(); ++i)
          {
            if(triSet[i] == tIdx)
            {
              blocks.push_back(leafblk);
              found = true;
            }
          }
        }
      }
    }
    sstr << "_blocks";
    dumpOctreeMeshBlocks(sstr.str(), blocks, true);
  }

  void dumpLocalOctreeMeshesForBlock(const std::string& name,
                                     const BlockIndex& block) const
  {
    // Dump a mesh with the single block
    std::stringstream blockStr;
    blockStr << name << "block_" << block.pt()[0] << "_" << block.pt()[1] << "_"
             << block.pt()[2];
    std::vector<BlockIndex> blocks;
    blocks.push_back(block);
    dumpOctreeMeshBlocks(blockStr.str(), blocks, true);

    // Dump a mesh with the indexed triangles for this block
    const InOutBlockData& blkData = m_octree[block];
    if(blkData.isLeaf() && blkData.hasData())
    {
      std::vector<TriangleIndex> tris;
      TriangleIndexSet triSet = m_octree.leafTriangles(block, blkData);
      for(int i = 0; i < triSet.size(); ++i)
      {
        tris.push_back(triSet[i]);
      }

      dumpTriangleMesh(blockStr.str() + "_triangles", tris, true);
    }
  }

  /** Generates a hexahedral VTK mesh with all octree blocks   */
  void dumpOctreeMeshVTK(const std::string& name) const
  {
    std::vector<BlockIndex> blocks;

    // Create an stl vector of all leaf blocks
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        if(it->isLeaf()) blocks.push_back(BlockIndex(it.pt(), lev));
      }
    }
    SLIC_INFO("Dump vtk:: Octree has " << blocks.size() << " leaves.");

    dumpOctreeMeshBlocks(name, blocks, false);
  }

  /** Generates a VTK mesh with all triangles in the mesh */
  void dumpTriMeshVTK(const std::string& name) const
  {
    const int numElts = m_octree.m_meshWrapper.numMeshElements();

    std::vector<TriangleIndex> tris;
    tris.reserve(numElts);

    for(int i = 0; i < numElts; ++i)
    {
      tris.push_back(i);
    }
    SLIC_INFO("Dump vtk:: Mesh has " << numElts << " triangles.");

    dumpTriangleMesh(name, tris, false);
  }

private:
  void dumpOctreeMeshBlocks(const std::string& name,
                            const std::vector<BlockIndex>& blocks,
                            bool shouldLogBlocks = false) const
  {
    if(blocks.empty()) return;

    // Dump an octree mesh containing all blocks
    std::stringstream fNameStr;
    fNameStr << name << ".vtk";

    DebugMesh* debugMesh = new DebugMesh(3, 8 * blocks.size(), blocks.size());
    const bool hasTriangles =
      (m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED);
    const bool hasColors =
      (m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED);

    // Allocate Slam Maps for the field data
    slam::PositionSet<> leafSet(blocks.size());

    LeafVertMap leafVertID(&leafSet);
    LeafVertMap leafVertID_unique(&leafSet);
    LeafIntMap leafTriCount(&leafSet);
    LeafIntMap leafColors(&leafSet);
    LeafIntMap leafLevel(&leafSet);
    LeafGridPtMap leafPoint(&leafSet);

    // Iterate through blocks -- and set the field data
    int leafCount = 0;
    for(auto it = blocks.begin(); it < blocks.end(); ++it)
    {
      const BlockIndex& block = *it;

      // Add the hex to the mesh
      addOctreeBlock(debugMesh, block, shouldLogBlocks);

      const InOutBlockData& leafData = m_octree[block];

      int vIdx = leafData.hasData() ? m_octree.leafVertex(block, leafData)
                                    : InOutOctreeType::MeshWrapper::NO_VERTEX;

      leafVertID[leafCount] = vIdx;
      leafLevel[leafCount] = block.level();
      leafPoint[leafCount] = block.pt();

      if(hasTriangles)
      {
        leafVertID_unique[leafCount] = m_octree.blockIndexesVertex(vIdx, block)
          ? vIdx
          : InOutOctreeType::MeshWrapper::NO_VERTEX;

        leafTriCount[leafCount] = leafData.hasData()
          ? m_octree.leafTriangles(block, leafData).size()
          : 0;
      }

      if(hasColors)
      {
        leafColors[leafCount] = m_colorsMap.at(leafData.color());
      }

      leafCount++;
    }

    // Add the fields to the mint mesh
    axom::IndexType* vertID = addIntField(debugMesh, "vertID");
    axom::IndexType* lLevel = addIntField(debugMesh, "level");

    axom::IndexType* blockCoord[3];
    blockCoord[0] = addIntField(debugMesh, "block_x");
    blockCoord[1] = addIntField(debugMesh, "block_y");
    blockCoord[2] = addIntField(debugMesh, "block_z");

    for(int i = 0; i < leafSet.size(); ++i)
    {
      vertID[i] = leafVertID[i];
      lLevel[i] = leafLevel[i];

      blockCoord[0][i] = leafPoint[i][0];
      blockCoord[1][i] = leafPoint[i][1];
      blockCoord[2][i] = leafPoint[i][2];
    }

    if(hasTriangles)
    {
      axom::IndexType* uniqVertID = addIntField(debugMesh, "uniqVertID");
      axom::IndexType* triCount = addIntField(debugMesh, "triCount");

      for(int i = 0; i < leafSet.size(); ++i)
      {
        uniqVertID[i] = leafVertID_unique[i];
        triCount[i] = leafTriCount[i];
      }
    }

    if(hasColors)
    {
      axom::IndexType* colors = addIntField(debugMesh, "colors");
      for(int i = 0; i < leafSet.size(); ++i) colors[i] = leafColors[i];
    }

    mint::write_vtk(debugMesh, fNameStr.str());

    delete debugMesh;
    debugMesh = nullptr;
  }

  void dumpTriangleMesh(const std::string& name,
                        const std::vector<TriangleIndex>& tris,
                        bool shouldLogTris = false) const
  {
    std::stringstream fNameStr;
    fNameStr << name << ".vtk";

    DebugMesh* debugMesh = new DebugMesh(3, 3 * tris.size(), tris.size());

    for(auto it = tris.begin(); it < tris.end(); ++it)
    {
      TriangleIndex tIdx = *it;
      addTriangle(debugMesh, tIdx, shouldLogTris);
    }

    // Add fields to the triangle mesh
    int numTris = tris.size();

    // Index of each triangle within the mesh
    axom::IndexType* triIdx = addIntField(debugMesh, "triangle_index");

    // Indices of the three boundary vertices of this triangle
    axom::IndexType* vertIdx[3];
    vertIdx[0] = addIntField(debugMesh, "vertex_index_0");
    vertIdx[1] = addIntField(debugMesh, "vertex_index_1");
    vertIdx[2] = addIntField(debugMesh, "vertex_index_2");

    for(int i = 0; i < numTris; ++i)
    {
      TriangleIndex tIdx = tris[i];
      triIdx[i] = tIdx;

      TriVertIndices tv = m_octree.m_meshWrapper.triangleVertexIndices(tIdx);
      vertIdx[0][i] = tv[0];
      vertIdx[1][i] = tv[1];
      vertIdx[2][i] = tv[2];
    }

    // other possible fields on triangles
    // -- number of blocks that index this triangle (blockCount)?
    // -- vertex field for number of triangles incident in the vertex (vtCount)?

    mint::write_vtk(debugMesh, fNameStr.str());

    delete debugMesh;
    debugMesh = nullptr;
  }

private:
  axom::IndexType* addIntField(DebugMesh* mesh, const std::string& name) const
  {
    axom::IndexType* fld =
      mesh->createField<axom::IndexType>(name, mint::CELL_CENTERED);
    SLIC_ASSERT(fld != nullptr);
    return fld;
  }

  double* addRealField(DebugMesh* mesh, const std::string& name) const
  {
    double* fld = mesh->createField<double>(name, mint::CELL_CENTERED);
    SLIC_ASSERT(fld != nullptr);
    return fld;
  }

  void addTriangle(DebugMesh* mesh,
                   const TriangleIndex& tIdx,
                   bool shouldLogTris) const
  {
    SpaceTriangle triPos = m_octree.m_meshWrapper.trianglePositions(tIdx);

    axom::IndexType vStart = mesh->getNumberOfNodes();
    mesh->appendNode(triPos[0][0], triPos[0][1], triPos[0][2]);
    mesh->appendNode(triPos[1][0], triPos[1][1], triPos[1][2]);
    mesh->appendNode(triPos[2][0], triPos[2][1], triPos[2][2]);

    axom::IndexType data[3];
    for(int i = 0; i < 3; ++i)
    {
      data[i] = vStart + i;
    }

    mesh->appendCell(data, mint::TRIANGLE);

    // Log the triangle info as primal code to simplify adding a test for this case
    if(shouldLogTris)
    {
      SLIC_INFO("// Triangle " << tIdx << "\n\t"
                               << "TriangleType tri("
                               << "PointType::make_point" << triPos[0] << ","
                               << "PointType::make_point" << triPos[1] << ","
                               << "PointType::make_point" << triPos[2] << ");");
    }
  }

  void addOctreeBlock(DebugMesh* mesh,
                      const BlockIndex& block,
                      bool shouldLogBlocks) const
  {
    GeometricBoundingBox blockBB = m_octree.blockBoundingBox(block);

    axom::IndexType vStart = mesh->getNumberOfNodes();

    const SpacePt& bMin = blockBB.getMin();
    const SpacePt& bMax = blockBB.getMax();

    mesh->appendNode(bMin[0], bMin[1], bMin[2]);
    mesh->appendNode(bMax[0], bMin[1], bMin[2]);
    mesh->appendNode(bMax[0], bMax[1], bMin[2]);
    mesh->appendNode(bMin[0], bMax[1], bMin[2]);

    mesh->appendNode(bMin[0], bMin[1], bMax[2]);
    mesh->appendNode(bMax[0], bMin[1], bMax[2]);
    mesh->appendNode(bMax[0], bMax[1], bMax[2]);
    mesh->appendNode(bMin[0], bMax[1], bMax[2]);

    axom::IndexType data[8];
    for(int i = 0; i < 8; ++i) data[i] = vStart + i;

    mesh->appendCell(data, mint::HEX);

    // Log bounding box info to simplify adding a test for this case
    if(shouldLogBlocks)
    {
      static int counter = 0;
      SLIC_INFO("// Block index " << block);
      SLIC_INFO("BoundingBoxType box"
                << ++counter << "(PointType::make_point" << bMin << ","
                << "PointType::make_point" << bMax << ");");
    }
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;

  ColorsMap m_colorsMap;
};

template <int DIM>
class InOutOctreeValidator
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using TriangleIndexSet = typename InOutOctreeType::TriangleIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using SpacePt = typename InOutOctreeType::SpacePt;
  using VertexIndex = typename InOutOctreeType::VertexIndex;
  using TriangleIndex = typename InOutOctreeType::TriangleIndex;
  using TriVertIndices = typename InOutOctreeType::MeshWrapper::TriVertIndices;
  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;

public:
  InOutOctreeValidator(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
  { }

  void checkAllLeavesColored() const
  {
    SLIC_DEBUG(
      "--Checking that all leaves have a color -- black, white and gray");

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      checkAllLeavesColoredAtLevel(lev);
    }
  }

  void checkAllLeavesColoredAtLevel(int level) const
  {
    const auto& levelLeafMap = m_octree.getOctreeLevel(level);
    auto itEnd = levelLeafMap.end();
    for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
    {
      if(!it->isLeaf()) continue;

      SLIC_ASSERT_MSG(it->isColored(),
                      "Error after coloring level "
                        << level << " leaf block " << BlockIndex(it.pt(), level)
                        << " was not colored.");
    }
  }

  void checkEachVertexIsIndexed() const
  {
    SLIC_DEBUG("--Checking that each vertex is in a leaf block of the tree.");

    const VertexIndex numVertices = m_octree.m_meshWrapper.numMeshVertices();
    for(VertexIndex i = 0; i < numVertices; ++i)
    {
      const SpacePt& pos = m_octree.m_meshWrapper.vertexPosition(i);
      BlockIndex vertBlock = m_octree.findLeafBlock(pos);
      const InOutBlockData& leafData = m_octree[vertBlock];

      VertexIndex vertInBlock = m_octree.leafVertex(vertBlock, leafData);
      AXOM_DEBUG_VAR(vertInBlock);

      // Check that we can find the leaf block indexing each vertex from its
      // position
      SLIC_ASSERT_MSG(
        leafData.hasData() && vertInBlock == i,
        " Vertex "
          << i << " at position " << pos << " \n\t was not indexed in block "
          << vertBlock << " with bounding box "
          << m_octree.blockBoundingBox(vertBlock) << " ( point is"
          << (m_octree.blockBoundingBox(vertBlock).contains(pos) ? "" : " NOT")
          << " contained in block )."
          << " \n\n *** \n Leaf data: " << leafData << " \n ***");

      // Check that our cached value of the vertex's block is accurate
      SLIC_ASSERT_MSG(
        vertBlock == m_octree.m_vertexToBlockMap[i],
        "Cached block for vertex "
          << i << " differs from its found block. "
          << "\n\t -- cached block " << m_octree.m_vertexToBlockMap[i]
          << "-- is leaf? "
          << (m_octree[m_octree.m_vertexToBlockMap[i]].isLeaf() ? "yes" : "no")
          << "\n\t -- actual block " << vertBlock << "-- is leaf? "
          << (m_octree[vertBlock].isLeaf() ? "yes" : "no")
          << "\n\t -- vertex set size: " << numVertices);
    }
  }

  void checkTrianglesReferencedInBoundaryVertexBlocks() const
  {
    SLIC_DEBUG(
      "--Checking that each triangle is referenced by the leaf blocks "
      "containing its vertices.");

    const axom::IndexType numTriangles = m_octree.m_meshWrapper.numMeshElements();
    for(axom::IndexType tIdx = 0; tIdx < numTriangles; ++tIdx)
    {
      TriVertIndices tvRel = m_octree.m_meshWrapper.triangleVertexIndices(tIdx);
      for(int j = 0; j < tvRel.size(); ++j)
      {
        VertexIndex vIdx = tvRel[j];
        BlockIndex vertBlock = m_octree.m_vertexToBlockMap[vIdx];
        const InOutBlockData& leafData = m_octree[vertBlock];

        // Check that this triangle is referenced here.
        bool foundTriangle = false;
        TriangleIndexSet leafTris = m_octree.leafTriangles(vertBlock, leafData);
        for(int k = 0; !foundTriangle && k < leafTris.size(); ++k)
        {
          if(leafTris[k] == tIdx) foundTriangle = true;
        }

        SLIC_ASSERT_MSG(foundTriangle,
                        "Did not find triangle with index "
                          << tIdx << " and vertices" << tvRel << " in block "
                          << vertBlock << " containing vertex " << vIdx
                          << " \n\n *** \n Leaf data: " << leafData
                          << " \n\t Triangles in block? " << leafTris
                          << " \n ***");
      }
    }
  }

  void checkBlockIndexingConsistency() const
  {
    // Check that internal blocks have no triangle / vertex
    //       and leaf blocks satisfy the conditions above
    SLIC_DEBUG(
      "--Checking that internal blocks have no data, and that leaves satisfy "
      "all PM conditions");

    const double bb_scale_factor = m_octree.m_boundingBoxScaleFactor;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(!data.isLeaf())
        {
          SLIC_ASSERT(!data.hasData());
        }
        else  // leaf block
        {
          if(data.hasData())
          {
            VertexIndex vIdx = m_octree.leafVertex(block, data);
            TriangleIndexSet triSet = m_octree.leafTriangles(block, data);
            for(int i = 0; i < triSet.size(); ++i)
            {
              TriangleIndex tIdx = triSet[i];

              // Check that vIdx is one of this triangle's vertices
              TriVertIndices tvRel =
                m_octree.m_meshWrapper.triangleVertexIndices(tIdx);

              SLIC_ASSERT_MSG(
                m_octree.m_meshWrapper.incidentInVertex(tvRel, vIdx),
                "All triangles in a gray block must be incident on a common "
                "vertex,"
                  << " but triangles " << tIdx << " with vertices " << tvRel
                  << " in block " << block << " is not incident in vertex "
                  << vIdx);

              // Check that this triangle intersects the bounding box of the
              // block
              GeometricBoundingBox blockBB = m_octree.blockBoundingBox(block);
              blockBB.expand(bb_scale_factor);
              SLIC_ASSERT_MSG(
                m_octree.blockIndexesElementVertex(tIdx, block) ||
                  intersect(m_octree.m_meshWrapper.trianglePositions(tIdx),
                            blockBB),
                "Triangle " << tIdx << " was indexed in block " << block
                            << " but it does not intersect the block."
                            << "\n\tBlock bounding box: " << blockBB
                            << "\n\tTriangle positions: "
                            << m_octree.m_meshWrapper.trianglePositions(tIdx)
                            << "\n\tTriangle vertex indices: " << tvRel
                            << "\n\tLeaf vertex is: " << vIdx
                            << "\n\tLeaf triangles: " << triSet << "("
                            << triSet.size() << ")");
            }
          }
        }
      }
    }
  }

  void checkNeighboringBlockColors() const
  {
    SLIC_DEBUG("--Checking that inside blocks do not neighbor outside blocks");
    for(int lev = m_octree.maxLeafLevel() - 1; lev >= 0; --lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);
      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const BlockIndex block(it.pt(), lev);
        const InOutBlockData& data = *it;

        if(data.isLeaf() && data.color() != InOutBlockData::Gray)
        {
          // Leaf does not yet have a color... try to find its color from
          // same-level face neighbors
          for(int i = 0; i < block.numFaceNeighbors(); ++i)
          {
            BlockIndex neighborBlk =
              m_octree.coveringLeafBlock(block.faceNeighbor(i));
            if(neighborBlk != BlockIndex::invalid_index())
            {
              const InOutBlockData& neighborData = m_octree[neighborBlk];
              switch(neighborData.color())
              {
              case InOutBlockData::Black:  // intentional
              // fallthrough
              case InOutBlockData::White:
                SLIC_CHECK_MSG(
                  data.color() == neighborData.color(),
                  "Problem at block "
                    << block << " with data " << data << " --- neighbor "
                    << neighborBlk << " with data "
                    << neighborData << ". Neighboring leaves that are not gray must have the same color.");
                break;
              case InOutBlockData::Gray:  // intentional
              // fallthrough
              case InOutBlockData::Undetermined:
                break;
              }
            }
          }
        }
      }
    }
  }

  void checkValid() const
  {
    // We are assumed to be valid before we insert the vertices
    if(m_generationState < InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
      return;

    // Iterate through the tree
    // Internal blocks should not have associated vertex data
    // Leaf block consistency depends on 'color'
    //      Black or White blocks should not have vertex data
    //      Gray blocks should have a vertex reference; it may or may not be
    // located within the block
    //          The sum of vertices located within a block should equal the
    // number of mesh vertices.
    //      Gray blocks should have one or more triangles.
    //          All triangles should be incident in a common vertex -- which
    // equals the indexed vertex, if it exists.

    if(m_generationState > InOutOctreeType::INOUTOCTREE_VERTICES_INSERTED)
    {
      checkEachVertexIsIndexed();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      checkTrianglesReferencedInBoundaryVertexBlocks();
      checkBlockIndexingConsistency();
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      checkAllLeavesColored();
      checkNeighboringBlockColors();
    }
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;
};

template <int DIM>
class InOutOctreeStats
{
public:
  using InOutOctreeType = InOutOctree<DIM>;
  using TriangleIndexSet = typename InOutOctreeType::TriangleIndexSet;

  using OctreeBaseType = typename InOutOctreeType::OctreeBaseType;
  using OctreeLevels = typename OctreeBaseType::OctreeLevels;
  using BlockIndex = typename OctreeBaseType::BlockIndex;

  using LeafCountMap = slam::Map<slam::Set<>, int>;
  using TriCountMap = slam::Map<slam::Set<>, int>;
  using CardinalityVTMap = slam::Map<slam::Set<>, int>;

  using LogHistogram = std::map<int, int>;
  using MinMaxRange = primal::BoundingBox<double, 1>;
  using LengthType = MinMaxRange::PointType;
  using LogRangeMap = std::map<int, MinMaxRange>;

  /** A simple struct to track totals within the octree levels */
  struct Totals
  {
    /** Default constructor to set everything to 0 */
    Totals()
      : blocks(0)
      , leaves(0)
      , leavesWithVert(0)
      , triangleRefCount(0)
      , whiteBlocks(0)
      , blackBlocks(0)
      , grayBlocks(0)
    { }

    int blocks;
    int leaves;
    int leavesWithVert;
    int triangleRefCount;
    int whiteBlocks;
    int blackBlocks;
    int grayBlocks;
  };

public:
  InOutOctreeStats(const InOutOctreeType& octree)
    : m_octree(octree)
    , m_generationState(m_octree.m_generationState)
    , m_levelBlocks(&m_octree.m_levels)
    , m_levelLeaves(&m_octree.m_levels)
    , m_levelLeavesWithVert(&m_octree.m_levels)
    , m_levelTriangleRefCount(&m_octree.m_levels)
    , m_levelWhiteBlockCount(&m_octree.m_levels)
    , m_levelBlackBlockCount(&m_octree.m_levels)
    , m_levelGrayBlockCount(&m_octree.m_levels)
  {
    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      m_triCount = TriCountMap(&m_octree.m_meshWrapper.elementSet());
    }

    // Iterate through blocks -- count the numbers of internal and leaf blocks
    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      const auto& levelLeafMap = m_octree.getOctreeLevel(lev);

      m_levelBlocks[lev] = levelLeafMap.numBlocks();
      m_levelLeaves[lev] = levelLeafMap.numLeafBlocks();
      m_levelLeavesWithVert[lev] = 0;
      m_levelTriangleRefCount[lev] = 0;
      m_levelWhiteBlockCount[lev] = 0;
      m_levelBlackBlockCount[lev] = 0;
      m_levelGrayBlockCount[lev] = 0;

      auto itEnd = levelLeafMap.end();
      for(auto it = levelLeafMap.begin(); it != itEnd; ++it)
      {
        const InOutBlockData& blockData = *it;
        BlockIndex block(it.pt(), lev);

        if(blockData.isLeaf())
        {
          if(blockData.hasData())
          {
            ++m_levelLeavesWithVert[lev];

            if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
            {
              m_levelTriangleRefCount[lev] +=
                m_octree.leafTriangles(block, blockData).size();

              BlockIndex blk(it.pt(), lev);
              TriangleIndexSet tris = m_octree.leafTriangles(blk, blockData);
              for(int i = 0; i < tris.size(); ++i)
              {
                ++m_triCount[tris[i]];
              }
            }
          }

          if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
          {
            switch(blockData.color())
            {
            case InOutBlockData::Black:
              ++m_levelBlackBlockCount[lev];
              break;
            case InOutBlockData::White:
              ++m_levelWhiteBlockCount[lev];
              break;
            case InOutBlockData::Gray:
              ++m_levelGrayBlockCount[lev];
              break;
            case InOutBlockData::Undetermined:
              break;
            }
          }
        }
      }

      m_totals.blocks += m_levelBlocks[lev];
      m_totals.leaves += m_levelLeaves[lev];
      m_totals.leavesWithVert += m_levelLeavesWithVert[lev];
      m_totals.triangleRefCount += m_levelTriangleRefCount[lev];
      m_totals.whiteBlocks += m_levelWhiteBlockCount[lev];
      m_totals.blackBlocks += m_levelBlackBlockCount[lev];
      m_totals.grayBlocks += m_levelGrayBlockCount[lev];
    }
  }

  /** Generates a string summarizing information about the leaves and blocks of
     the octree */
  std::string blockDataStats() const
  {
    std::stringstream sstr;

    for(int lev = 0; lev < m_octree.m_levels.size(); ++lev)
    {
      if(m_levelBlocks[lev] > 0)
      {
        int percentWithVert =
          integerPercentage(m_levelLeavesWithVert[lev], m_levelLeaves[lev]);

        sstr << "\t Level " << lev << " has " << m_levelBlocks[lev]
             << " blocks -- " << m_levelBlocks[lev] - m_levelLeaves[lev]
             << " internal; " << m_levelLeaves[lev] << " leaves "
             << " (" << percentWithVert << "% w/ vert); ";

        if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
        {
          sstr << " Leaves with colors -- B,W,G ==> "
               << m_levelBlackBlockCount[lev] << ","
               << m_levelWhiteBlockCount[lev] << ","
               << m_levelGrayBlockCount[lev] << " and "
               << m_levelTriangleRefCount[lev] << " triangle references.";
        }
        //sstr <<"Hash load factor: "
        //     << this->m_leavesLevelMap[ lev ].load_factor()
        //     << " -- max lf: " << this->m_leavesLevelMap[ lev
        // ].max_load_factor();
        sstr << "\n";
      }
    }

    return sstr.str();
  }

  /** Generates a string summarizing information about the mesh elements indexed
     by the octree */
  std::string meshDataStats() const
  {
    std::stringstream sstr;

    double meshNumTriangles = m_octree.m_meshWrapper.numMeshElements();

    int percentWithVert =
      integerPercentage(m_totals.leavesWithVert, m_totals.leaves);
    sstr << "  Mesh has " << m_octree.m_meshWrapper.numMeshVertices()
         << " vertices."
         << "\n  Octree has " << m_totals.blocks << " blocks; "
         << m_totals.blocks - m_totals.leaves << " internal; "
         << m_totals.leaves << " leaves "
         << " (" << percentWithVert << "% w/ vert); ";

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      sstr << " \n\t There were " << m_totals.triangleRefCount
           << " triangle references "
           << " (avg. " << (m_totals.triangleRefCount / meshNumTriangles)
           << " refs per triangle).";
    }

    return sstr.str();
  }

  std::string triangleCountHistogram() const
  {
    std::stringstream sstr;

    // Generate and output a histogram of the bucket counts on a lg-scale
    LogHistogram triCountHist;  // Create histogram of edge lengths (log
                                // scale)
    LogRangeMap triCountRange;

    int numElems = m_octree.m_meshWrapper.numMeshElements();

    for(int i = 0; i < numElems; ++i)
    {
      LengthType count(m_triCount[i]);
      int expBase2;
      std::frexp(m_triCount[i], &expBase2);
      triCountHist[expBase2]++;
      triCountRange[expBase2].addPoint(count);
    }

    std::stringstream triCountStr;
    triCountStr << "\tTriangle index count "
                << "(lg-arithmic bins for number of references per triangle):";
    for(auto it = triCountHist.begin(); it != triCountHist.end(); ++it)
    {
      triCountStr << "\n\t exp: " << it->first << "\t count: " << (it->second)
                  << "\tRange: " << triCountRange[it->first];
    }

    return triCountStr.str();
  }

  std::string vertexCardinalityHistogram() const
  {
    std::stringstream sstr;

    using TriVertIndices = typename InOutOctreeType::MeshWrapper::TriVertIndices;

    // Generate and output histogram of VT relation
    CardinalityVTMap cardVT(&m_octree.m_meshWrapper.vertexSet());

    int numElems = m_octree.m_meshWrapper.numMeshElements();
    for(int i = 0; i < numElems; ++i)
    {
      TriVertIndices tvRel = m_octree.m_meshWrapper.triangleVertexIndices(i);
      cardVT[tvRel[0]]++;
      cardVT[tvRel[1]]++;
      cardVT[tvRel[2]]++;
    }

    using LinHistogram = std::map<int, int>;
    LinHistogram vtCardHist;
    int numVerts = m_octree.m_meshWrapper.numMeshVertices();
    for(int i = 0; i < numVerts; ++i)
    {
      LengthType count(cardVT[i]);
      vtCardHist[cardVT[i]]++;
    }

    sstr << "\tCardinality VT relation histogram (linear): ";
    for(LinHistogram::const_iterator it = vtCardHist.begin();
        it != vtCardHist.end();
        ++it)
    {
      sstr << "\n\t exp: " << it->first << "\t count: " << (it->second);
    }

    return sstr.str();
  }

  std::string summaryStats() const
  {
    std::stringstream octreeStatsStr;

    octreeStatsStr << "*** "
                   << (m_generationState >=
                           InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED
                         ? "PM"
                         : "PR")
                   << " octree summary *** \n";

    octreeStatsStr << blockDataStats() << "\n" << meshDataStats();

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_LEAVES_COLORED)
    {
      octreeStatsStr << "\n\tColors B,W,G ==> " << m_totals.blackBlocks << ","
                     << m_totals.whiteBlocks << "," << m_totals.grayBlocks;
    }

    if(m_generationState >= InOutOctreeType::INOUTOCTREE_ELEMENTS_INSERTED)
    {
      octreeStatsStr << "\n"
                     << triangleCountHistogram() << "\n"
                     << vertexCardinalityHistogram();
    }

    return octreeStatsStr.str();
  }

private:
  int integerPercentage(double val, double size) const
  {
    return (size > 0) ? static_cast<int>((100. * val) / size) : 0;
  }

private:
  const InOutOctreeType& m_octree;
  typename InOutOctreeType::GenerationState m_generationState;

  LeafCountMap m_levelBlocks;
  LeafCountMap m_levelLeaves;
  LeafCountMap m_levelLeavesWithVert;
  LeafCountMap m_levelTriangleRefCount;

  LeafCountMap m_levelWhiteBlockCount;
  LeafCountMap m_levelBlackBlockCount;
  LeafCountMap m_levelGrayBlockCount;

  TriCountMap m_triCount;

  Totals m_totals;
};

}  // end namespace detail

}  // end namespace quest
}  // end namespace axom

#endif  // SPATIAL_OCTREE__HXX_
