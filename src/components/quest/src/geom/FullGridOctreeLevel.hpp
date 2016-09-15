#ifndef FULL_GRID_OCTREE_LEVEL__HXX_
#define FULL_GRID_OCTREE_LEVEL__HXX_


#include "common/config.hpp"
#include "common/CommonTypes.hpp"

#include "quest/OctreeLevel.hpp"


namespace quest
{

    /**
     * \class
     * \brief A FullGridOctreeLevel is a representation of an OctreeLevel.
     *  that allocates space for all possible blocks at the given level.
     *
     *  It uses a Brood-based organization of the data, where the data for all octree
     *  siblings are stored contiguously, and uses a Morton-based order on the Broods
     *  in the level.  Associated with each block is a BlockDataType
     *  For efficiency, the data is associated with an entire brood, a collection of
     *  siblings that are created simultaneously. In dimension DIM, there are 2^DIM siblings in a brood.
     *
     *  \see OctreeLevel
     *  \see Brood
     */
    template<int DIM, typename MortonIndexType, typename BlockDataType>
    class FullGridOctreeLevel : public OctreeLevel<DIM,BlockDataType>
    {
    public:
        typedef OctreeLevel<DIM, BlockDataType> Base;
        typedef typename Base::GridPt                   GridPt;
        typedef typename Base::BroodData                BroodData;
        typedef typename Base::BlockIteratorHelper      BaseBlockIteratorHelper;
        typedef typename Base::ConstBlockIteratorHelper ConstBaseBlockIteratorHelper;

        template<typename OctreeLevelType, typename ParentType> class FullGridBlockIteratorHelper;

        typedef FullGridBlockIteratorHelper<FullGridOctreeLevel, BaseBlockIteratorHelper> FullGridBlockIterHelper;
        typedef FullGridBlockIteratorHelper<const FullGridOctreeLevel, ConstBaseBlockIteratorHelper> ConstFullGridBlockIterHelper;

        typedef Brood<GridPt, MortonIndexType> BroodType;

    public:

        /**
         * \brief Concrete instance of the BlockIteratorHelper class defined in the OctreeLevel base class.
         * \note ParenType must be either BlockIteratorHelper or ConstBlockIteratorHelper,
         *       both are defined in the OctreeLevel base class
         */
        template<typename OctreeLevelType, typename ParentType>
        class FullGridBlockIteratorHelper : public ParentType
        {
        public:
            typedef FullGridBlockIteratorHelper<OctreeLevelType, ParentType> FullGridBlockItType;
            typedef ParentType     BaseBlockItType;

            FullGridBlockIteratorHelper(OctreeLevelType* octLevel, bool begin)
                : m_octreeLevel(octLevel)
                , m_endIdx( octLevel->m_broodCapacity)
                , m_offset(0)
                , m_isLevelZero( octLevel->level() == 0)
            {
                m_currentIdx = begin? 0 : m_endIdx;

                // Advance the iterator to point to a valid Block
                if( begin && !data()->isBlock())
                    increment();
            }

            /** Increment to next block of the level */
            void increment()
            {
                // Note, must skip blocks that are not in the tree
                do {
                    ++m_offset;
                    if(m_offset == Base::BROOD_SIZE || m_isLevelZero)
                    {
                        ++m_currentIdx;
                        m_offset = 0;
                    }
                } while( m_currentIdx < m_endIdx && !data()->isBlock());
            }

            /** Access to point associated with the block pointed to by the iterator */
            GridPt pt() const
            {
                // Reconstruct the grid point from its brood representation
                typedef Mortonizer<typename GridPt::CoordType, MortonIndexType, GridPt::NDIMS> Mort;
                return Mort::demortonize( (m_currentIdx << DIM)  + m_offset);
            }

            /** Access to data associated with the block pointed to by the iterator */
            BlockDataType* data() { return &m_octreeLevel->m_data[m_currentIdx][m_offset]; }
            /** Const access to data associated with the block pointed to by the iterator */
            const BlockDataType* data() const { return &m_octreeLevel->m_data[m_currentIdx][m_offset]; }

            /** Determine if two BlockIterators are pointing to the same block */
            bool equal(const BaseBlockItType* other)
            {
                const FullGridBlockItType* pother = dynamic_cast<const FullGridBlockItType*>(other);

                return (pother != ATK_NULLPTR)
                     && (m_currentIdx == pother->m_currentIdx)   // iterators are the same
                     && (m_offset == pother->m_offset);               // brood indices are the same
            }
        private:
            OctreeLevelType* m_octreeLevel;
            MortonIndexType m_currentIdx, m_endIdx;
            int m_offset;
            bool m_isLevelZero;
        };

    public:

        /** \brief Default constructor for an octree level */
        FullGridOctreeLevel(int level = -1): Base(level), m_blockCount(0)
        {
            if(level < 0)
            {
                m_broodCapacity = 0;
                m_data = ATK_NULLPTR;
            }
            else
            {
                m_broodCapacity = 1;

                if(level > 0)
                {
                    const int rowsize = 1 << (level-1);
                    for(int i=0; i < DIM; ++i)
                        m_broodCapacity *= rowsize;
                }

                m_data = new BroodData[m_broodCapacity];
            }

            // Mark all blocks as non blocks
            for(int i=0; i<m_broodCapacity; ++i)
            {
                BroodData& bd = m_data[i];
                for(int j=0; j< Base::BROOD_SIZE; ++j)
                    bd[j].setNonBlock();
            }
        }

        ~FullGridOctreeLevel()
        {
            if(m_data != ATK_NULLPTR)
            {
                delete [] m_data;
                m_data = ATK_NULLPTR;
            }

            m_broodCapacity = 0;
            m_blockCount = 0;
        }

        /**
         * \brief Factory function to return a GridBlockIterHelper for this level
         *  \param begin A boolean to determine if this is to be a begin (true) or end (false) iterator
         */
        BaseBlockIteratorHelper* getIteratorHelper(bool begin)
        {
            return new FullGridBlockIterHelper(this, begin);
        }

        /**
         * \brief Factory function to return a ConstGridBlockIterHelper for this level
         *  \param begin A boolean to determine if this is to be a begin (true) or end (false) iterator
         */
        ConstBaseBlockIteratorHelper* getIteratorHelper(bool begin) const
        {
            return new ConstFullGridBlockIterHelper(this, begin);
        }


        /**
         * \brief Predicate to check whether the block associated with the given GridPt pt is in the current level
         */
        bool hasBlock(const GridPt& pt) const
        {
            const BroodType brood(pt);
            return m_data[brood.base()][brood.offset()].isBlock();
        }

        /**
         * \brief Adds all children of the given grid point to the octree level
         * \param [in] pt The gridPoint associated with the parent of the children that are being added
         * \pre pt must be in bounds for the level
         * \sa inBounds()
         */
        void addAllChildren(const GridPt& pt)
        {
            SLIC_ASSERT_MSG(this->inBounds(pt)
                           , "Problem while inserting children of point " << pt
                           << " into octree level " << this->m_level
                           << ". Point was out of bounds -- "
                           << "each coordinate must be between 0 and " << this->maxCoord() << ".");

            getBroodData(pt) = BroodData();

            // Handle level 0 -- only add the root, mark its 'siblings' as non-blocks
            if( this->level() == 0)
            {
                for(int j=1; j< Base::BROOD_SIZE; ++j)
                    m_data[0][j].setNonBlock();
                ++m_blockCount;
            }
            else
            {
                m_blockCount += Base::BROOD_SIZE;
            }
        }


        /** \brief Accessor for the data associated with pt */
        BlockDataType& operator[](const GridPt& pt)
        {
            const BroodType brood(pt);
            return m_data[brood.base()][brood.offset()];
        }

        /** \brief Const accessor for the data associated with pt */
        const BlockDataType& operator[](const GridPt& pt) const
        {
            const BroodType brood(pt);
            return m_data[brood.base()][brood.offset()];
        }

        /** \brief Access the data associated with the entire brood */
        BroodData& getBroodData(const GridPt& pt)
        {
            return m_data[BroodType::MortonizerType::mortonize(pt)];
        }

        /** \brief Const access to data associated with the entire brood */
        const BroodData& getBroodData(const GridPt& pt) const
        {
            return m_data[BroodType::MortonizerType::mortonize(pt)];
        }


        /** \brief Predicate to check if there are any blocks in this octree level */
        bool empty() const { return m_blockCount == 0; }

        /** \brief Returns the number of blocks (internal and leaf) in the level */
        int numBlocks() const { return m_blockCount; }

        /** \brief Returns the number of internal blocks in the level */
        int numInternalBlocks() const { return numBlocks() - numLeafBlocks(); }

        /** \brief Returns the number of leaf blocks in the level */
        int numLeafBlocks() const
        {
            if(empty())
                return 0;

            int count = 0;
            for(int i=0; i< m_broodCapacity; ++i)
            {
                const BroodData& bd  = m_data[i];
                for(int j=0; j< Base::BROOD_SIZE; ++j)
                {
                    if(bd[j].isLeaf())
                         ++count;
                }
            }
            return count;
        }


        /**
         * \brief Helper function to determine the status of an octree block within this octree level
         * \param pt The grid point of the block index that we are testing
         * \return The status of the grid point pt (e.g. LeafBlock, InternalBlock, ...)
         */
        TreeBlockStatus blockStatus(const GridPt & pt) const
        {
            if(! this->inBounds(pt))
                return BlockNotInTree;

            const BroodType brood(pt);
            const BlockDataType& blockData =  m_data[brood.base()][brood.offset()];

            return blockData.isBlock()
                    ? (blockData.isLeaf() ? LeafBlock : InternalBlock )
                    : BlockNotInTree;
        }

    //private:
    //  DISABLE_COPY_AND_ASSIGNMENT(OctreeLevel);

    private:
      BroodData* m_data;
      int m_broodCapacity;
      int m_blockCount;
    };

} // end namespace quest

#endif  // FULL_GRID_OCTREE_LEVEL__HXX_
