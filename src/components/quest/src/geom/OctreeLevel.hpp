
#ifndef OCTREE_LEVEL__HXX_
#define OCTREE_LEVEL__HXX_

#include <ostream>   // for ostream in print

#include "quest/BoundingBox.hpp"
#include "quest/MortonIndex.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include "common/config.hpp"

#include "slic/slic.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"


#include "fmt/fmt.hpp"

#ifdef USE_CXX11
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

#include <boost/iterator/iterator_facade.hpp>

/**
 * \file
 * \brief Defines templated OctreeLevel class
 */

namespace quest
{

    /**
     * \brief Helper enumeration for status of a BlockIndex within an OctreeLevel instance
     */
    enum TreeBlock { BlockNotInTree, LeafBlock, InternalBlock};

    template<int DIM, typename BlockDataType>
    class OctreeLevel
    {
    public:
        typedef int CoordType;
        typedef quest::Point<CoordType,DIM> GridPt;
        typedef quest::Vector<CoordType,DIM> GridVec;

    private:
        /**
         * \class
         * \brief Private inner class to handle subindexing of block data within siblings
         */
      struct Brood {
          enum { NUM_CHILDREN = 1 << DIM };

          Brood(const GridPt& pt)
              : m_broodPt( pt.array() /2), m_idx(0)
          {
              for(int i=0; i< DIM; ++i)
              {
                  m_idx |= (pt[i]& 1) << i;
              }
          }

          const GridPt& pt() const { return m_broodPt; }
          const int& idx() const { return m_idx; }

          GridPt m_broodPt;
          int m_idx;
      };

    public:

        // A brood is a collection of sibling blocks that are created at the same time
        typedef BlockDataType BroodData[ Brood::NUM_CHILDREN ];


      #if defined(USE_CXX11)
        typedef std::unordered_map<GridPt, BroodData, PointHash<int> > MapType;
      #else
        typedef boost::unordered_map<GridPt, BroodData, PointHash<int> > MapType;
      #endif

        typedef typename MapType::iterator MapIter;
        typedef typename MapType::const_iterator ConstMapIter;

        template<typename OctreeLevel, typename InnerIterType, typename DataType> class BlockIterator;

        typedef BlockIterator<OctreeLevel, MapIter, BlockDataType> BlockIter;
        typedef BlockIterator<const OctreeLevel, ConstMapIter, const BlockDataType> ConstBlockIter;


    public:

        /**
         * \class
         * \brief An iterator type for the blocks of an octree level
         */
        template<typename OctreeLevel, typename InnerIterType, typename DataType>
        class BlockIterator : public boost::iterator_facade< BlockIterator<OctreeLevel, InnerIterType, DataType>
                                   , DataType
                                   , boost::forward_traversal_tag
                                   , DataType
                                   >
        {
        public:
          typedef BlockIterator<OctreeLevel, InnerIterType, DataType>              iter;

          BlockIterator(OctreeLevel* octLevel, bool begin = false)
              : m_octLevel(octLevel), m_data(ATK_NULLPTR), m_idx(0)
          {
              SLIC_ASSERT(octLevel != ATK_NULLPTR);

              m_endIter = m_octLevel->m_map.end();
              m_levelIter = begin ? m_octLevel->m_map.begin() : m_endIter;

              update();
          }

		 // valid for const access to the data
          const DataType& dereference() const { return *m_data; }

          const GridPt& pt() const { return m_pt; }
          
          // Use for non-const access to the data
          DataType& data() { return *m_data; }
          const DataType& data() const { return *m_data; }

          bool equal(const iter& other) const
          {
              return (m_octLevel == other.m_octLevel)       // point to same object
                   && (m_levelIter == other.m_levelIter)   // iterators are the same
                   && (m_idx == other.m_idx);               // brood indices are the same
          }

          void increment()
          {
              ++m_idx;

              if(m_idx == Brood::NUM_CHILDREN || m_octLevel->m_level == 0)
              {
                  ++m_levelIter;
                  m_idx = 0;
              }

              update();
          }
          void update()
          {
              // Test to see if we have finished iterating
              if(m_levelIter == m_endIter)
                  return;

              const GridPt& itPt = m_levelIter->first;
              for(int i=0; i<DIM; ++i)
                  m_pt[i] = (itPt[i]<<1) + ( m_idx & (1 << i)? 1 : 0);

              m_data = &m_levelIter->second[m_idx];
          }

        private:
          friend class boost::iterator_core_access;
          OctreeLevel*  m_octLevel;
          InnerIterType m_levelIter, m_endIter;
          GridPt        m_pt;
          DataType*     m_data;
          int           m_idx;
        };

    public:

        OctreeLevel(int level = -1): m_level(level){}

        //void setLevel(int level) { m_level = level; }


        CoordType maxCoord() const
        {
            return (1<< m_level) -1;
        }

        GridPt maxGridCell() const
        {
            return GridPt(maxCoord());
        }

        bool isLeaf(const GridPt& pt) const { return blockStatus(pt) == LeafBlock; }
        bool isInternal(const GridPt& pt) const { return blockStatus(pt) == InternalBlock; }
        bool hasBlock(const GridPt& pt) const
        {
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.pt());
            return blockIt != m_map.end();
        }

        bool inBounds(const GridPt& pt) const
        {
            const CoordType maxVal = maxCoord();
            for(int i=0; i< DIM; ++i)
                if( pt[i] < 0 || pt[i] > maxVal)
                    return false;
            return true;
        }

        BlockDataType& operator[](const GridPt& pt)
        {
            const Brood brood(pt);
            return m_map[ brood.pt() ][brood.idx()];
        }
        const BlockDataType& operator[](const GridPt& pt) const
        {
            SLIC_ASSERT_MSG(hasBlock(pt)
                            ,"(" << pt <<", "<< m_level << ") was not a block in the tree at level.");

            // Note: Using find() method on hashmap since operator[] is non-const
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.pt());
            return blockIt->second[brood.idx()];
        }

        /**
         * \brief Helper function to determine the status of an octre block within an octree level
         * \note This function is meant to help with implementing basic octree functionality
         *       and is not meant to be exposed in the public API
         * \param pt The grid point of the block index that we are testing
         * \param lev The level of the block index that we are testing
         */
        TreeBlock blockStatus(const GridPt & pt) const
        {
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.pt());

            return (blockIt == m_map.end())
                    ? BlockNotInTree
                    : (blockIt->second[brood.idx()].isLeaf())
                        ? LeafBlock
                        : InternalBlock;
        }

        BlockIter      begin()       { return BlockIter(this,true); }
        ConstBlockIter begin() const { return ConstBlockIter(this,true); }
        BlockIter      end()         { return BlockIter(this,false); }
        ConstBlockIter end()   const { return ConstBlockIter(this,false); }

        bool empty() const { return m_map.empty(); }

    //private:
    //  DISABLE_COPY_AND_ASSIGNMENT(OctreeLevel);

    private:
      int m_level;
      MapType m_map;
    };


} // end namespace quest

#endif  // OCTREE_BASE_HXX_
