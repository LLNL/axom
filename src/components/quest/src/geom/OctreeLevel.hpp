
#ifndef OCTREE_LEVEL__HXX_
#define OCTREE_LEVEL__HXX_

#include <ostream>   // for ostream in print

#include "quest/BoundingBox.hpp"
#include "quest/MortonIndex.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/Mesh.hpp"

#include "common/config.hpp"

#include "slic/slic.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/Map.hpp"

#ifdef USE_CXX11
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

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

      #if defined(USE_CXX11)
        typedef std::unordered_map<GridPt, BlockDataType, PointHash<int> > MapType;
      #else
        typedef boost::unordered_map<GridPt, BlockDataType, PointHash<int> > MapType;
      #endif

        typedef typename MapType::iterator MapIterator;
        typedef typename MapType::const_iterator ConstMapIterator;

    public:

        OctreeLevel(int level = -1): m_level(level) {}

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
            ConstMapIterator blockIt = m_map.find(pt);
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
            return m_map[ pt ];
        }
        const BlockDataType& operator[](const GridPt& pt) const
        {
            SLIC_ASSERT_MSG(hasBlock(pt)
                            ,"(" << pt <<", "<< m_level << ") was not a block in the tree at level.");

            // Note: Using find() method on hashmap since operator[] is non-const
            ConstMapIterator blockIt = m_map.find(pt);
            return blockIt->second;
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
            ConstMapIterator blockIt = m_map.find(pt);

            if(blockIt == m_map.end())
                return BlockNotInTree;

            return (blockIt->second.isLeaf()) ? LeafBlock: InternalBlock;
        }

        MapIterator      begin()       { return m_map.begin(); }
        ConstMapIterator begin() const { return m_map.begin(); }
        MapIterator      end()         { return m_map.end(); }
        ConstMapIterator end()   const { return m_map.end(); }

        bool empty() const { return m_map.empty(); }

    //private:
    //  DISABLE_COPY_AND_ASSIGNMENT(OctreeLevel);

    private:
      int m_level;
      MapType m_map;
    };


} // end namespace quest

#endif  // OCTREE_BASE_HXX_
