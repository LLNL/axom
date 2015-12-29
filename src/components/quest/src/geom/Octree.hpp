
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

  enum{ MAX_LEV = 30 };

  typedef quest::Point<int,DIM> GridPt;
  typedef quest::Vector<int,DIM> GridVec;

  typedef asctoolkit::slam::policies::CompileTimeSizeHolder<int, MAX_LEV> MAX_LEVEL_SIZE;
  typedef asctoolkit::slam::OrderedSet<MAX_LEVEL_SIZE> OctreeLevels;

#if defined(USE_CXX11)
  typedef std::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
#else
  typedef boost::unordered_map<GridPt, LeafNodeType, PointHash<int> > MapType;
#endif
  typedef asctoolkit::slam::Map<MapType> LeafIndicesLevelMap;

  typedef std::pair< const GridPt, int> BlockIndex;


public:
  TopologicalOctree()
    : m_leavesLevelMap(&m_levels)
  {
      m_leavesLevelMap[0][GridPt()] = LeafNodeType(); //GridPt(), 0);
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

  /**
   * \brief Determine whether the octree contains a leaf block associated with grid point @pt at level @level
   * \returns true if the associated block is a leaf in the octree, false otherwise
   */
  bool leafExists(const GridPt& pt, int level) const
  {
      const LeafIndicesLevelMap& levelLeafMap = m_leavesLevelMap[ m_levels[level] ];
      return levelLeafMap.find(pt) != levelLeafMap.end();
  }

  BlockIndex root() const {
      return std::make_pair( GridPt(), 0);
  }

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
      return std::make_pair( GridPt( pt.array() /2 ), m_levels[level -1]);
  }

  BlockIndex child(const GridPt & pt, int level, int childIndex) const
  {
      SLIC_ASSERT_MSG(0 <= childIndex && childIndex < (1<<DIM), "Child index must be between 0 and " << (1 << DIM) );

      const GridPt childPt = GridPt( pt.array() * 2 + childOffsetFromIndex(childIndex).array() );
      int childLevel = m_levels[level +1];

      return std::make_pair(childPt, childLevel);
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


private:
  /**
   * \brief Helper function to compute the child offset from an index
   * \param [in] childIndex The index of the child to find
   * \pre 0 <= @childIndex < 2^DIM
   */
  GridPt childOffsetFromIndex(int childIndex) const
  {
      GridPt childOffset;

      for(int dim=0; dim < DIM; ++dim)
      {
          // Use the convention that each bit in the provided childIndex parameter
          // determines whether we chose the left or right child in that dimension
          if(childIndex & (1 << dim))
          {
              ++childOffset[dim];
          }
      }

      return childOffset;

  }

private:
  OctreeLevels            m_levels;
  LeafIndicesLevelMap     m_leavesLevelMap;
};





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
        typedef quest::Point<int,DIM> GridPt;

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


} // end namespace quest

#endif  // OCTREE_HXX_
