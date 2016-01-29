
#ifndef OCTREE_BASE__HXX_
#define OCTREE_BASE__HXX_

#include <ostream>   // for ostream in print

#include "quest/BoundingBox.hpp"
#include "quest/MortonIndex.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/Mesh.hpp"

#include "slic/slic.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"

#if defined(USE_CXX11)
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

/**
 * \file
 * \brief Defines templated OctreeBase class and its inner class BlockIndex
 */


namespace quest
{

  /**
   * \brief Minimal implementation of a BlockDataType for an OctreeBase.
   * LeafData is default constructible and provides an isLeaf() and setInternal() functions.
   * \note This implementation uses ones-complement to differentiate between leaf and internal blocks.
   * This has the nice properties that (a) all internal LeafData have an id with a negative number,
   * (b) the LeafData id's can be zero-based (c) all LeafDatas -- representing leaf and internal blocks
   * can live in the same index space.
   */
  class LeafData
  {
  public:
      LeafData()
      {
          static int idGenerator = 0;
          m_id = idGenerator++;
      }

      LeafData(const LeafData& other): m_id(other.m_id) {}
      LeafData& operator=(const LeafData& other) { m_id = other.m_id; return *this;}

      /**
       * \brief Predicate to determine if this LeafData represents a leaf block in the tree
       */
      bool isLeaf() const { return m_id >= 0; }

      /**
       * \brief Sets the data associated with this leaf
       */
      void setData(int leafID) { m_id = leafID; }

      /**
       * Returns the normalized form of the id for this LeafData instance
       * \note The normalized form is a non-negative integer.
       */
      int getID() const { return isLeaf()? m_id : ~m_id; }

      /**
       * \brief Predicate to determine if this LeafData is associated with a leaf block
       * or internal block of the octree.
       * \return True, if the associated block is internal, false otherwise
       */
      void setInternal() { if(isLeaf()) { m_id = ~m_id; } }

      /**
       * Equality operator for comparing two LeafData instances
       */
      friend bool operator==(const LeafData& lhs, const LeafData& rhs )
      {
          return lhs.m_id == rhs.m_id;
      }

  protected:
      int m_id;
  };




/**
 * \class
 * \brief Handles the non-geometric operations for our octree such as refinement,
 * finding the parents and children of a node and determining whether a leaf node exists
 *
 * \note There are two concepts here.
 *      A set of nested dyadic integer grids -- A grid at level n has 2^n cells in each of the DIM dimensions
 *      and an adaptive octree defined by a subset of the grid cells from these nested grids.
 *      The former is a conceptual aid in the sense that it is implicitly encoded.
 *
 *      An octree block at level l is refined by adding its 2^DIM children blocks at level (l+1).
 *      The root of the octree covers the entire domains and has no parent.
 *      The leaf blocks of the octree have no children.
 *      The interior of the children do not overlap, and their union covers that of their parent block.
 *      Non-leaf blocks are referred to as 'internal'
 *
 * \note Requirements for BlockDataType: it must be default constructible and provide an isLeaf() predicate as well
 * as a setInternal() function that changes its state from representing a leaf block to representing an internal
 * block.
 */
template<int DIM, typename BlockDataType>
class OctreeBase
{
public:
  typedef int CoordType;
  typedef quest::Point<CoordType,DIM> GridPt;
  typedef quest::Vector<CoordType,DIM> GridVec;

  typedef asctoolkit::slam::policies::CompileTimeSizeHolder<CoordType, std::numeric_limits<CoordType>::digits> MAX_LEVEL_SIZE;
  typedef asctoolkit::slam::OrderedSet<MAX_LEVEL_SIZE> OctreeLevels;


#if defined(USE_CXX11)
  typedef std::unordered_map<GridPt, BlockDataType, PointHash<int> > MapType;
#else
  typedef boost::unordered_map<GridPt, BlockDataType, PointHash<int> > MapType;
#endif
  typedef typename MapType::iterator LevelMapIterator;
  typedef typename MapType::const_iterator LevelMapCIterator;

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
  private:

      enum  {
          /** The number of children of an octree block (\f$ 2^{DIM} \f$ in dimension DIM ) */
          NUM_CHILDREN = 1 << DIM
      };

      typedef asctoolkit::slam::policies::CompileTimeSizeHolder<int, NUM_CHILDREN> OCTREE_CHILDREN_SIZE;

  public:
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

        std::ostream& print(std::ostream& os) const
        {
            os << "{pt: " << m_pt
                <<"; level: " << m_lev <<"}";

            return os;
        }

        friend std::ostream& operator<<(std::ostream& os, const BlockIndex& block)
        {
            block.print(os);
            return os;
        }


    /**
     * \brief Helper function to generate an invalid block index.
     * \return A new BlockIndex instance blk
     * \post  blk.isValid() will return false
     */
      static BlockIndex invalid_index() { return BlockIndex( GridPt(), -1); }

      /**
       * \brief The number of children that an octree block can have
       */
      static int numChildren() { return ChildIndexSet().size(); }


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
      m_leavesLevelMap[rootBlock.level()][rootBlock.pt()] = BlockDataType();
  }

  /**
   * \brief The resolution octree level for leaf blocks of the octree
   */
  int maxLeafLevel() const { return m_levels.size(); }

  /**
   * \brief The resolution octree level for internal blocks of the octree
   */
  int maxInternalLevel() const { return m_levels.size()-1; }
public:
   // \todo KW Convert these two functions to static class functions.
   //        This will require converting m_levels to a static set

  //@{

  /**
   * \brief Utility function to find the number of (possible) grid cells at a given level or resolution
   * \param [in] level The level or resolution.
   * \pre \f$ 0 \le lev < \f$ maxLeafLevel()
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
      LevelMapCIterator blockIt = levelLeafMap.find(pt);

      return (blockIt != levelLeafMap.end() && blockIt->second.isLeaf() );
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
   * \brief Refines the given leaf block in the octree
   * Marks leafBlock as internal (non-leaf) and adds its children to the tree
   * \pre leafBlock is a valid leaf block in the octree.
   */
  void refineLeaf(const BlockIndex & leafBlock)
  {
    SLIC_ASSERT( isLeaf(leafBlock) );

    // Find the leaf node and set as internal
    MapType& currentNodeLevelMap = m_leavesLevelMap[ m_levels[leafBlock.level()] ];
    currentNodeLevelMap[ leafBlock.pt() ].setInternal();

    // Add its children to the tree
    MapType& childLevelMap = m_leavesLevelMap[ m_levels[leafBlock.childLevel()] ];
    typedef typename BlockIndex::ChildIndexSet ChildIndexSet;
    const int numChildren = ChildIndexSet().size();
    for(int childIdx=0; childIdx < numChildren; ++childIdx)
    {
        childLevelMap[ leafBlock.childPt(childIdx) ] = BlockDataType();
    }
  }

  /**
   * Accessor to the data associated with leafBlock
   * \param leafBlock A valid leaf block in the tree
   * \pre leafBlock is a valid leaf in the tree
   */
  BlockDataType& operator[](const BlockIndex& leafBlock)
  {
      SLIC_ASSERT(isLeaf(leafBlock));

      return m_leavesLevelMap[ m_levels[leafBlock.level()] ][ leafBlock.pt()];
  }

  /**
   * \brief const accessor to leaf data associated with leafBlock
   * \param leafBlock A valid leaf block in the octree
   * \pre leafBlock is a valid leaf in the tree
   */
  const BlockDataType& operator[](const BlockIndex& leafBlock) const
  {
      SLIC_ASSERT( isLeaf(leafBlock) );

      // Note: Using find() method on hashmap since operator[] is non-const
      const MapType& levelLeafMap = m_leavesLevelMap[ m_levels[leafBlock.level()] ];
      LevelMapCIterator blockIt = levelLeafMap.find(leafBlock.pt());

      return blockIt->second;
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(OctreeBase)

protected:
  OctreeLevels            m_levels;
  LeafIndicesLevelMap     m_leavesLevelMap;
};


} // end namespace quest

#endif  // OCTREE_BASE_HXX_
