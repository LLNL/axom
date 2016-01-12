
#ifndef OCTREE_HXX_
#define OCTREE_HXX_

#include "quest/BitTwiddle.hpp"
#include "quest/BoundingBox.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include "slic/slic.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"

#include "quest/Mesh.hpp"

#if defined(USE_CXX11)
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

#include <utility>      // for std::pair's operator==()


namespace quest
{

/**
 * \class
 * \brief Handles the non-geometric operations for our octree such as refinement,
 * finding the parents and children of a node and determining whether a leaf node exists
 *
 */
template<int DIM, typename LeafNodeType>
class TopologicalOctree
{
public:

  enum{ MAX_LEV = 30
      , NUM_CHILDREN = 1 << DIM
  };

  typedef quest::Point<int,DIM> GridPt;
  typedef quest::Vector<int,DIM> GridVec;

  typedef asctoolkit::slam::policies::CompileTimeSizeHolder<int, MAX_LEV> MAX_LEVEL_SIZE;
  typedef asctoolkit::slam::OrderedSet<MAX_LEVEL_SIZE> OctreeLevels;


#if defined(USE_CXX11)
  typedef std::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
#else
  typedef boost::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
#endif
  typedef typename MapType::iterator LevelMapIterator;

  typedef asctoolkit::slam::Map<MapType> LeafIndicesLevelMap;

  /**
   * \brief Inner class encapsulating the index of an octree <em>block</em>.
   *
   * Each block index is represented as a point on an integer grid (the minimum point of the block's extent)
   * at a given level of resolution.
   *
   * Each level of resolution is a regular grid with \f$ 2^{level} \f$
   * grid points along each dimension.  The <em>root</em> block (at level 0) covers the entire domain.
   * An octree block at level \f$ \ell \f$ has \f$ 2^{DIM} \f$ <em>children</em> at level \f$ \ell + 1 \f$
   * covering its domain.
   */
  class BlockIndex {
  public:
      typedef asctoolkit::slam::policies::CompileTimeSizeHolder<int, NUM_CHILDREN> OCTREE_CHILDREN_SIZE;
      typedef asctoolkit::slam::OrderedSet<OCTREE_CHILDREN_SIZE> ChildIndexSet;

  public:
      /**
       * \brief Default constructor
       */
      BlockIndex() : m_pt( GridPt() ), m_lev(0) {}

      /**
       * \brief Constructor from a point and a level
       */
      BlockIndex(const GridPt& pt, int level) : m_pt(pt), m_lev(level) {}

      /**
       * \brief Accessor for the BlockIndex instance's point
       * \returns const reference to the instance's point
       */
      const GridPt& pt() const  { return m_pt; }

      /**
       * \brief Accessor for the BlockIndex instance's point
       * \returns reference to the instance's point
       */
      GridPt& pt()        { return m_pt; }

      /**
       * \brief Accessor for the BlockIndex instance's level
       * \returns const reference to the instance's level
       */
      const int& level() const  { return m_lev; }

      /**
       * \brief Accessor for the BlockIndex instance's level
       * \returns reference to the instance's level
       */
      int& level()        { return m_lev; }

      /**
       * \brief The level of the block index's parent
       */
       int   parentLevel() const { return m_lev -1; }

      /**
       * \brief The level of the block index's child
       */
      int   childLevel()  const { return m_lev +1; }


      /**
       * \brief Returns the grid point of the block's parent
       */
      GridPt parentPt() const
      {
          return GridPt( m_pt.array() /2);
      }

        /**
         * \brief Returns the grid point of the block's child at index childIndex
         * \param [in] childIndex The index of the child whose grid point we are finding
         * \pre \f$ 0 \le childIndex < \f$ Octree::NUM_CHILDREN
         */
        GridPt childPt(int childIndex) const
        {
            SLIC_ASSERT( ChildIndexSet().isValidIndex(childIndex) );

            GridPt cPoint;

            // Child is at next level of resolution (multiply by two)
            // and offset according to whether corresponding bit
            // in child index is set or not
            for(int dim =0; dim< DIM; ++dim)
            {
                cPoint[dim] = (m_pt[dim] << 1)
                          + (childIndex & (1 << dim) ? 1 : 0);
            }

            return cPoint;
        }

        /**
         * \brief Returns the parent BlockIndex of this block
         * \note Returns an invalid BlockIndex if we attempt to find
         *       the parent of the root block
         */
        BlockIndex parent() const
        {
            return BlockIndex(parentPt(), parentLevel());
        }

        /**
         * \brief Returns the child BlockIndex of this block
         * \param [in] childIndex The index of the child whose grid point we are finding
         * \pre \f$ 0 \le childIndex < \f$ Octree::NUM_CHILDREN
         */
        BlockIndex child(int childIndex) const
        {
            return BlockIndex( childPt(childIndex), childLevel());
        }


        bool operator==(const BlockIndex& other) const {
            return (m_lev == other.m_lev) && (m_pt == other.m_pt);
        }

        bool operator!=(const BlockIndex& other) const {
            return !(*this==other);
        }

        /**
         * \brief Checks the validity of the index.
         * A block index is valid when its level is \f$ \ge 0 \f$
         * and each coordinate p[i] of its grid point is \f$ 0 \le p[i] < 2^{level()} \f$.
         *
         * \returns true if the block index is valid, else false
         */
        bool isValid() const
        {
            bool bValid = (m_lev >= 0);

            for(int i = 0; i < DIM; ++i)
                bValid = bValid && (m_pt[i] >=0) && (m_pt[i] < (1<< m_lev) );

            return bValid;
        }


    /**
     * \brief Helper function to generate an invalid block index.
     * \return A new BlockIndex instance blk
     * \post  blk.isValid() will return false
     */
      static BlockIndex invalid_index() { return BlockIndex( GridPt(), -1); }

  private:
      GridPt m_pt;
      int    m_lev;
  };


public:
  /**
   * \brief Default constructor.
   * Sets up an octree containing only the root block
   */
  TopologicalOctree()
    : m_leavesLevelMap(&m_levels)
  {
      BlockIndex rootBlock = root();
      m_leavesLevelMap[rootBlock.level()][rootBlock.pt()] = LeafNodeType();
  }

public:
   // \todo KW Convert these two functions to static class functions.
   //        This will require converting m_levels to a static set

  //@{

  /**
   * \brief Utility function to find the number of (possible) grid cells at a given level or resolution
   * \param [in] level The level or resolution.
   * \pre \f$ 0 \le lev < \text{TopologicalOctree::MAX_LEV}  \f$
   * \todo Convert this to a static class function.
   */
  GridPt maxGridCellAtLevel(int level) const
  {
    return GridPt(1<< m_levels[level] );
  }

  /**
   * Auxiliary function to return the root of the octree
   * \note The root block has no parent.
   *       Its parent is an invalid BlockIndex.
   *       I.e. octree.parent( octree.root()).isValid() = false.
   * \todo Convert this to a static class function
   */
  BlockIndex root() const { return BlockIndex(); }

  // @}

public:
  // @{
  // KW: The following four functions are probably not necessary any more
  //     Since their functionality is in the BlockIndex inner class.


  /**
   * \brief Finds the grid index and level of the current octree block's parent.
   * \note The root level is 0 and its children are at level 1
   * \note The root node has no parent.  The returned level will be '-1'
   * \param [in] pt The grid index of the block whose parent we want to find.
   * \param [in] level The level of the block whose parent we want to find.
   * \return The parent of the provided octree leaf.
   */
  BlockIndex parent(const GridPt & pt, int level) const
  {
      return BlockIndex(pt, level).parent();
  }

  /**
   * \brief Finds the BlockIndex of the given block's parent.
   * \param [in] block The block whose parent we want to find
   * \return The BlockIndex of the parent of the provided octree leaf.
   */
  BlockIndex parent(const BlockIndex& block) const
  {
      return block.parent();
  }


  /**
   * \brief Finds the BlockIndex of the given block's child
   * \param [in] pt The grid index of the block whose child we want to find.
   * \param [in] level The level of the block whose child we want to find.
   * \param [in] childIndex The index of the child to find
   * \pre \f$ 0 \le childIndex < 2^{DIM} \f$
   * \return The BlockIndex of the child of the provided octree leaf.
   */
  BlockIndex child(const GridPt & pt, int level, int childIndex) const
  {
      return BlockIndex(pt,level).child(childIndex);
  }

  /**
   * \brief Finds the BlockIndex of the given block's child.
   * \param [in] block The block whose child we want to find
   * \param [in] childIndex The index of the child to find
   * \pre \f$ 0 \le childIndex < 2^{DIM} \f$
   * \return The BlockIndex of the child of the provided octree leaf.
   */
  BlockIndex child(const BlockIndex& block, int childIndex) const
  {
      return block.child(childIndex);
  }

  // @}


public:
  /**
   * \brief Determine whether the octree contains a leaf block associated with grid point pt at level level
   * \param [in] pt The grid point to check
   * \param [in] level The level of the grid point
   * \returns true if the associated block is a leaf in the octree, false otherwise
   */
  bool isLeaf(const GridPt& pt, int level) const
  {
      const MapType& levelLeafMap = m_leavesLevelMap[ m_levels[level] ];
      return levelLeafMap.find(pt) != levelLeafMap.end();
  }

  /**
   * \brief Determine whether the octree contains a leaf block associated with grid point pt at level level
   * \param [in] block The BlockIndex of the tree to check
   * \returns true if the associated block is a leaf in the octree, false otherwise
   */
  bool isLeaf(const BlockIndex& block) const
  {
      return isLeaf(block.pt(), block.level());
  }


  /**
   * \brief Removes the given leaf block from the tree and adds its children
   * \pre leafBlock is a leaf block in the octree.  I.e. octree.leafExists(leafBlock)==true.
   * \note We might choose to simply mark leafBlock as a non-leaf rather than deleting
   * it from the tree.
   */
  void refineLeaf(const BlockIndex & leafBlock)
  {
    SLIC_ASSERT( isLeaf(leafBlock) );

    // 1. Find the leaf node
    // 2. Remove it from the tree (at the appropriate level
    MapType& currentNodeLevelMap = m_leavesLevelMap[ m_levels[leafBlock.level()] ];
    currentNodeLevelMap.erase(leafBlock.pt());

    // 3. Add its children to the tree
    MapType& childLevelMap = m_leavesLevelMap[ m_levels[leafBlock.childLevel()] ];

    typedef typename BlockIndex::ChildIndexSet ChildIndexSet;
    const int numChildren = ChildIndexSet().size();
    for(int childIdx=0; childIdx < numChildren; ++childIdx)
    {
        childLevelMap[ leafBlock.childPt(childIdx) ] = LeafNodeType();
    }

  }



private:
  DISABLE_COPY_AND_ASSIGNMENT(TopologicalOctree)

protected:
  OctreeLevels            m_levels;
  LeafIndicesLevelMap     m_leavesLevelMap;
};





/**
 * \class
 * \brief Adds spatial extents to a topological octree, allowing point location
 */
template<int DIM, typename LeafNodeType>
class SpatialOctree : public TopologicalOctree<DIM, LeafNodeType>
{
public:

    typedef quest::BoundingBox<double,DIM> GeometricBoundingBox;
    typedef quest::Point<double,DIM> SpacePt;
    typedef quest::Vector<double,DIM> SpaceVector;

    typedef typename TopologicalOctree<DIM, LeafNodeType>::GridPt GridPt;
    typedef typename TopologicalOctree<DIM, LeafNodeType>::MapType MapType;
    typedef typename TopologicalOctree<DIM, LeafNodeType>::BlockIndex BlockIndex;

//    typedef typename TopologicalOctree<DIM, LeafNodeType>::OctreeLevels OctreeLevels;

    typedef asctoolkit::slam::Map<SpaceVector> SpaceVectorLevelMap;

public:
    /**
     * \brief Construct a spatial octree from a spatial bounding box
     * \param [in] bb The spatial extent to be indexed by the octree
     */
    SpatialOctree(const GeometricBoundingBox& bb)
        : TopologicalOctree<DIM,LeafNodeType>()
        , m_deltaLevelMap(& this->m_levels)
        , m_boundingBox(bb)
    {
        SpaceVector bbRange = bb.range();
        for(int lev = 0; lev < this->m_levels.size(); ++lev)
        {
            m_deltaLevelMap[lev] = bbRange / static_cast<double>(1<<lev);
        }

    }

    /**
     * \brief Return the spatial bounding box of a grid cell at the given level or resolution
     */
    GeometricBoundingBox blockBoundingBox(const BlockIndex & block)
    {
        return blockBoundingBox( block.pt(), block.level() );
    }

    /**
     * \brief Return the spatial bounding box of a grid cell at the given level or resolution
     */
    GeometricBoundingBox blockBoundingBox(const GridPt & gridPt, int level)
    {
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];

        SpacePt lower(m_boundingBox.getMin());
        SpacePt upper(m_boundingBox.getMin());
        for(int i=0; i< DIM; ++i)
        {
            lower[i] +=  gridPt[i]   * deltaVec[i];
            upper[i] += (gridPt[i]+1)* deltaVec[i];
        }

        return GeometricBoundingBox( lower,upper );
    }


    /**
     * Returns the width of an octree block at level of resolution level
     */
    const SpaceVector& spacingAtLevel(int level) const
    {
        return m_deltaLevelMap[level];
    }

    /**
     * \brief Finds the index of the leaf block covering the query point pt
     * \param [in] pt The query point in space
     * \param [in] startingLevel Optional starting level for the query
     * \pre pt must be in the bounding box of the octree (i.e. boundingBox.contains(pt) == true )
     * \note The collection of leaves covers the bounding box, and the interiors of the leaves do not
     * intersect, so every point in the bounding box should be located in a unique leaf block.
     * \note We are assuming a half-open interval on the bounding boxes.
     * \return The block index (i.e. grid point and level) of the leaf block containing the query point
     */
    BlockIndex findLeafBlock(const SpacePt& pt, int startingLevel = 0) const
    {
        bool found = false;
        BlockIndex leafBlock = BlockIndex::invalid_index();

        SLIC_ASSERT( m_boundingBox.contains(pt) );

        for(int lev=startingLevel; !found && lev < this->m_levels.size(); ++lev)
        {
            GridPt gridPt = findGridCellAtLevel(pt, lev);
            found = this->isLeaf(gridPt, lev);

            if(found)
                leafBlock = BlockIndex(gridPt, lev);
        }

        SLIC_ASSERT_MSG(found, "Point " << pt << " not found "
                        <<"in a leaf block of the octree");

        return leafBlock;
    }

    /**
     * \brief Utility function to find the quantized level lev grid cell of Point pt
     * \param [in] pt The point at which we are querying.
     * \param [in] level The level or resolution.
     * \pre \f$ 0 \le lev < MAX_LEV == 32 \f$
     * \return The grid point of the block covering this point at this level
     * \todo KW: Should this function be protected? Is it generally useful?
     */
    GridPt findGridCellAtLevel(const SpacePt& pt, int level) const
    {
        SpaceVector ptVec(m_boundingBox.getMin(), pt);
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];

        // Find the octree block that covers us at this level
        // Note: we are assuming a half-open interval for all blocks
        //       that are not on the upper boundary of the domain
        GridPt quantizedPt = elementwiseQuantizedRatio( ptVec, deltaVec);

        // Fix for when we are on the upper boundary of the domain
        for(int i=0; i< DIM; ++i)
        {
            if( quantizedPt[i] == 1<<level )
                --quantizedPt[i];
        }

        return quantizedPt;
    }

private:
    /**
     * \brief Helper function to quantize to the integer grid
     */
    GridPt elementwiseQuantizedRatio(const SpaceVector& ptFromBBMin, const SpaceVector&  cellWidth) const
    {
        GridPt gridPt;
        for(int i=0; i< DIM; ++i)
            gridPt[i] = std::floor( ptFromBBMin[i] / cellWidth[i]);

        return gridPt;
    }

private:
  DISABLE_COPY_AND_ASSIGNMENT(SpatialOctree)

protected:
    SpaceVectorLevelMap     m_deltaLevelMap;
    GeometricBoundingBox    m_boundingBox;
};


namespace junkyard {
    // The code here was written quickly in an initial implementation.
    // It is currently being refactored into classes that better separate
    // the concerns:
    //      TopologicalOctree
    //      SpatialOctree
    //      Point and triangle mesh instantiations...

    /**
     * \class
     * \brief Simple encoding of all data we need for a leaf node
     *
     * \note This is a WIP implementation of an octree LeafNode
     *   that will be used during the octree construction.phase.
     *
     *   After getting a feel for the interface, I would like to convert this
     *   to using slam::DynamicSets (and use this as a motivating example to
     *   define/implement these).
     *
     *   Specifically, we will have a dynamic set of LeafNodes
     *   We will then have dynamic relations from LeafNodes
     *   to Mesh vertices and mesh triangles (cells)
     *   We will also have a map from leaf nodes to the color attribute (inside/outside/gray).
     *   These can be compacted / made static after construction for faster access
     *
     *   The hashmaps in the octree will be used to obtain the index of the corresponding leaf node
     *   which we will use as keys to the maps/relations.
     *
     */
    template<int DIM>
    struct LeafNode
    {
        typedef int VertexIndex;
        typedef int TriangleIndex;
        typedef std::vector<TriangleIndex> TriangleIndexVec;

        typedef Point<int, DIM> GridPt;

        enum { NO_VERTEX = -1 };

        enum LeafColor {
              UNDETERMINED_COLOR    // default state
            , BLACK                 // inside surface
            , WHITE                 // outside surface
            , GRAY                  // intersects geometry
            , INVALID_LEAF          // E.g. it became an internal node
        };

    public:

        LeafNode(GridPt gPt = GridPt(), int level=0)
            : m_gridPt(gPt)
            , m_level(level)
            , m_vert(NO_VERTEX)
            , m_color(UNDETERMINED_COLOR)
        {
        }

        // ------ Functions related to validity and refinement ------

        bool isValid() { return m_color != INVALID_LEAF; }

        static LeafNode make_invalid_leaf()
        {
              LeafNode ln;
              ln.m_color = INVALID_LEAF;
              return ln;
        }


        // ------ Functions related to vertices ------
        bool hasVertexId() const
        {
            return (m_vert != NO_VERTEX);
        }

        VertexIndex vertexIndex() const { return m_vert; }

        void addVertex(VertexIndex vIdx)
        {
            SLIC_ASSERT(!hasVertexId());

            m_vert = vIdx;
        }


        // ------ Functions related to triangles ------


        /**
         * The invariant of a leaf node in this octree
         * is that all triangles indexed by a leaf node
         * are incident in a common vertex.
         */
        bool allTrianglesShareAVertex(meshtk::Mesh* mesh)
        {
            SLIC_ERROR("Not implemented yet.");
            return false;
        }

        void addTriangle(TriangleIndex tIdx)
        {
            m_tris.push_back(tIdx);
        }

        /**
         * \brief Return the set of triangle indices referenced by this leaf node
         *
         * \note This will eventually be converted into a slam::set
         */
        const TriangleIndexVec& triangleIndices() const { return m_tris; }


    private:
        GridPt              m_gridPt;
        int                 m_level;
        VertexIndex         m_vert;
        TriangleIndexVec    m_tris;
        LeafColor           m_color;
    };


/**
 * \class
 * \brief A pointerless octree to aid in inside/outside point queries relative to a surface
 */
template<int DIM>
class Octree
{
public:

    enum{ MAX_LEV = 30
    };


    typedef quest::BoundingBox<double,DIM> GeometricBoundingBox;
    typedef quest::Point<double,DIM> SpacePt;
    typedef quest::Vector<double,DIM> SpaceVector;

    typedef Point<int,DIM> GridPt;

    typedef asctoolkit::slam::policies::CompileTimeSizeHolder<int, MAX_LEV> MAX_LEVEL_SIZE;
    typedef asctoolkit::slam::OrderedSet<MAX_LEVEL_SIZE> OctreeLevels;
    typedef asctoolkit::slam::Map<SpaceVector> SpaceVectorLevelMap;


    typedef int VertexIndex;
    typedef int TriangleIndex;

    typedef LeafNode<DIM> LeafNodeType;

    #if defined(USE_CXX11)
      typedef std::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
    #else
      typedef boost::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
    #endif
      typedef asctoolkit::slam::Map<MapType> LeafIndicesLevelMap;

public:
    Octree(const GeometricBoundingBox& bb)
        : m_boundingBox(bb)
        , m_deltaLevelMap(&m_levels)
        , m_leavesLevelMap(&m_levels)
    {
        SpaceVector bbRange = bb.range();
        for(int lev = 0; lev < m_levels.size(); ++lev)
        {
            m_deltaLevelMap[lev] = bbRange / static_cast<double>(1<<lev);
        }

        m_leavesLevelMap[0][GridPt()] = LeafNodeType(GridPt(), 0);
    }

    /**
     * \brief Utility function to find the quantized level lev grid cell of Point pt
     * \param [in] pt The point at which we are querying.
     * \param [in] level The level or resolution.
     * \pre \f$ 0 \le lev < MAX_LEV == 32 \f$
     */
    GridPt findGridCellAtLevel(const SpacePt& pt, int level)
    {
        SpaceVector ptVec(m_boundingBox.getMin(), pt);
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];

        return elementwiseQuantizedRatio( ptVec, deltaVec);
    }

    /**
     * \brief Find the spatial bounding box of a grid cell at the given level or resolution
     */
    GeometricBoundingBox leafCellBoundingBox(const GridPt & gridPt, int level)
    {
        const SpaceVector& deltaVec = m_deltaLevelMap[ level];
        SpaceVector lower(m_boundingBox.getMin());
        SpaceVector upper(m_boundingBox.getMin());
        for(int i=0; i< DIM; ++i)
        {
            lower[i] +=  gridPt[i]   * deltaVec[i];
            upper[i] += (gridPt[i]+1)* deltaVec[i];
        }

        return GeometricBoundingBox(lower,upper);
    }

    /**
     * \brief Utility function to find the number of (possible) grid cells at a given level or resolution
     * \param [in] The level or resolution.
     * \pre \f$ 0 \le lev < MAX_LEV == 32 \f$
     */
    GridPt maxGridCellAtLevel(int level) const
    {
        return GridPt(1<< m_levels[level] );
    }

    const SpaceVector& spacingAtLevel(int level) const
    {
        return m_deltaLevelMap[level];
    }

    void insertPoint(const SpacePt& pt, int level = 0)
    {
        // 1. Find leaf node
        LeafNodeType ln = findLeafNode(pt);

        // 2. Check if we can add the point directly
        //    e.g. it has no vertex, or the current vertex
        //    there can be melded with the one we are adding
        bool canAddVertex = false;
        bool hasVertex = ln.hasVertexId();
        if(hasVertex)
        {
            //canAddVertex = attemptMeldVertices( ln.vertexIndex());
            // needs an absolute and relative epsilon here!
        }
        else
        {
            canAddVertex = false;
        }

        // 3. If not,
        if(!canAddVertex )
        {
            // 3.a refine the leaf
           refineLeaf(ln);

           //  3.b. and try to add the vertex to the child
           insertPoint(pt, ln.level+1);
        }

    }

    void refineLeaf(const LeafNodeType & node)
    {
        // 1. Find the leaf node
        // 2. Remove it from the tree (at the appropriate level
        // 3. Add its children to the tree
        // 4. Add its geometry to the appropriate children
        //      e.g. its vertex and triangle indices

    }

    LeafNodeType findLeafNode(const GridPt& pt, int startingLevel = 0)
    {
        bool found = false;
        LeafNodeType foundNode;
        for(int lev=startingLevel; !found && lev < m_levels.size(); ++lev)
        {
            foundNode = findLeafNodeAtLevel(pt,lev);
            found = isValidLeafNode( foundNode );
        }

        return foundNode;
    }

    LeafNodeType findLeafNodeAtLevel(const GridPt& pt, int level) const
    {
        const MapType& leavesMap = m_leavesLevelMap[level];
        typename MapType::const_iterator leafIt;

        // Find the quantized grid point at this level of resolution
        GridPt gridPt = findGridCellAtLevel(pt, level);

        if( (leafIt = leavesMap.find(gridPt)) == leavesMap.end())
        {
            return LeafNodeType::make_invalid_leaf();
        }

        return leafIt->second;;
    }

    bool isValidLeafNode(const LeafNodeType& leaf)
    {
        return leaf.isValue();
    }

    void addTriangle(TriangleIndex triInd)
    {
        // 1. Find the bounding box of the triangle
        // 2. Descend the tree
        //    Starting at root, recursively check the children that overlap the bounding box
        //    For each such child, check if the associated leaf exists and if the triangle actually
        //    intersects the bounding box.
        // 3. Try to insert the triangle into each of these, with octree refinement, when necessary
    }


private:
    Octree(){}

    /**
     * \brief Helper function to quantize to the integer grid
     */
    GridPt elementwiseQuantizedRatio(const SpaceVector& num, const SpaceVector&  denom) const
    {
        GridPt gridPt;
        for(int i=0; i< GridPt::NDIMS; ++i)
            gridPt[i] = std::floor( num[i] / denom[i]);
        return gridPt;
    }

private:
    OctreeLevels            m_levels;
    GeometricBoundingBox    m_boundingBox;
    SpaceVectorLevelMap     m_deltaLevelMap;
    LeafIndicesLevelMap     m_leavesLevelMap;
};

} // end namespace junkyard


} // end namespace quest

#endif  // OCTREE_HXX_
