
#ifndef OCTREE_BASE__HXX_
#define OCTREE_BASE__HXX_

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
class OctreeBase
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
  OctreeBase()
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
  DISABLE_COPY_AND_ASSIGNMENT(OctreeBase)

protected:
  OctreeLevels            m_levels;
  LeafIndicesLevelMap     m_leavesLevelMap;
};







} // end namespace quest

#endif  // OCTREE_BASE_HXX_
