#ifndef MORTON_OCTREE_LEVEL__HXX_
#define MORTON_OCTREE_LEVEL__HXX_


#include "common/config.hpp"
#include "common/CommonTypes.hpp"

#include "quest/OctreeLevel.hpp"

#ifdef USE_CXX11
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

namespace quest
{
    /**
     * \class
     * \brief A MortonOctreeLevel is a representation of a sparse OctreeLevel.
     *  It is a concrete implementation of an OctreeLevel
     *
     *  It associates data with its Octree block using a hash map
     *  whose key type is a MortonIndex type (an integer type),
     *  and whose value type is a BlockDataType.
     *  For efficiency, the data is associated with an entire brood, a collection of
     *  siblings that are created simultaneously. In dimension DIM, there are 2^DIM siblings in a brood.
     *
     *  \see OctreeLevel
     */
    template<int DIM, typename MortonIndexType, typename BlockDataType>
    class MortonOctreeLevel : public OctreeLevel<DIM,BlockDataType>
    {
    public:
      typedef OctreeLevel<DIM, BlockDataType> Base;
      typedef typename Base::GridPt GridPt;
      typedef typename Base::BroodData BroodData;
      typedef typename Base::BlockIteratorHelper BaseBlockIteratorHelper;
      typedef typename Base::ConstBlockIteratorHelper ConstBaseBlockIteratorHelper;

      #if defined(USE_CXX11)
          typedef std::unordered_map<MortonIndexType, BroodData> MapType;
    #else
          typedef boost::unordered_map<MortonIndexType, BroodData> MapType;
      #endif

        typedef typename MapType::iterator       MapIter;
        typedef typename MapType::const_iterator ConstMapIter;

        template<typename OctreeLevelType, typename InnerIterType, typename ParentType> class MortonBlockIteratorHelper;

        typedef MortonBlockIteratorHelper<MortonOctreeLevel, MapIter, BaseBlockIteratorHelper> MortonBlockIterHelper;
        typedef MortonBlockIteratorHelper<const MortonOctreeLevel, ConstMapIter, ConstBaseBlockIteratorHelper> ConstMortonBlockIterHelper;

    protected:
        /**
         * \class
         * \brief Private helper class to handle subindexing of block data within octree siblings
         * \note A brood is a collection of siblings that are generated simultaneously.
         * \note This class converts a grid point at the given level into a brood index of the point.
         *       The base brood is the MortonIndex of the grid point's octree parent
         *       and its offset index is obtained by interleaving the least significant bit of its coordinates.
         */
      struct Brood {
          enum { BROOD_BITMASK = Base::BROOD_SIZE -1 };

          typedef Mortonizer<typename GridPt::CoordType, MortonIndexType, GridPt::NDIMS> MortonizerType;

          /**
           * \brief Constructor for a brood offset relative to the given grid point pt
           * \param [in] pt The grid point within the octree level
           */
          Brood(const GridPt& pt)
          {
              m_broodIdx = MortonizerType::mortonize(pt);
              m_offset = m_broodIdx & BROOD_BITMASK;
              m_broodIdx >>= DIM;
          }

          /** \brief Accessor for the base point of the entire brood */
          const MortonIndexType& base() const { return m_broodIdx; }

          /** \brief Accessor for the index of the point within the brood */
          const int& offset() const { return m_offset; }

      private:
          MortonIndexType m_broodIdx;   /** MortonIndex of the base point of all blocks within the brood */
          int m_offset;                 /** Index of the block within the brood. Value is in [0, 2^DIM) */
      };

    public:

       /**
        * \brief Concrete instance of the BlockIteratorHelper class defined in the OctreeLevel base class.
        */
        template<typename OctreeLevelType, typename InnerIterType, typename ParentType>
        class MortonBlockIteratorHelper : public ParentType
        {
        public:
            typedef MortonBlockIteratorHelper<OctreeLevelType, InnerIterType, ParentType> MortonBlockItType;
            typedef ParentType     BaseBlockItType;

            MortonBlockIteratorHelper(OctreeLevelType* octLevel, bool begin)
                : m_endIter( octLevel->m_map.end() )
                , m_offset(0)
                , m_isLevelZero( octLevel->level() == 0)
            {
                m_currentIter = begin ? octLevel->m_map.begin() : m_endIter;
            }

            /** Increment to next block of the level */
            void increment()
            {
                ++m_offset;

                if(m_offset == Base::BROOD_SIZE || m_isLevelZero)
                {
                    ++m_currentIter;
                    m_offset = 0;
                }
            }

            /** Access to point associated with the block pointed to by the iterator */
            GridPt pt() const
            {
                // Reconstruct the grid point from its brood representation
                typedef Mortonizer<typename GridPt::CoordType, MortonIndexType, GridPt::NDIMS> Mort;
                return Mort::demortonize( (m_currentIter->first << DIM)  + m_offset);
            }

            /** Access to data associated with the block pointed to by the iterator */
            BlockDataType* data() { return &m_currentIter->second[m_offset]; }
            /** Const access to data associated with the block pointed to by the iterator */
            const BlockDataType* data() const { return &m_currentIter->second[m_offset]; }

            /** Determine if two BlockIterators are pointing to the same block */
            bool equal(const BaseBlockItType* other)
            {
                const MortonBlockItType* pother = dynamic_cast<const MortonBlockItType*>(other);

                return (pother != ATK_NULLPTR)
                     && (m_currentIter == pother->m_currentIter)   // iterators are the same
                     && (m_offset == pother->m_offset);               // brood indices are the same
            }
        private:
            InnerIterType m_currentIter, m_endIter;
            int m_offset;
            bool m_isLevelZero;
        };

    public:

        /** \brief Default constructor for an octree level */
        MortonOctreeLevel(int level = -1): Base(level){}

        /** Factory function to return an octree block iterator helper */
        BaseBlockIteratorHelper* getIteratorHelper(bool begin)
        {
            return new MortonBlockIterHelper(this, begin);
        }

        /** Factory function to return a const octree block iterator helper */
        ConstBaseBlockIteratorHelper* getIteratorHelper(bool begin) const
        {
            return new ConstMortonBlockIterHelper(this, begin);
        }


        /**
         * \brief Predicate to check whether the block associated with the given GridPt pt is in the current level
         */
        bool hasBlock(const GridPt& pt) const
        {
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.base());
            return blockIt != m_map.end();
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

            m_map[ Brood::MortonizerType::mortonize(pt) ];  // Adds children, if not already present, using default BlockDataType() constructor
        }



        /** \brief Accessor for the data associated with pt */
        BlockDataType& operator[](const GridPt& pt)
        {
            const Brood brood(pt);
            return m_map[brood.base()][brood.offset()];
        }

        /** \brief Const accessor for the data associated with pt */
        const BlockDataType& operator[](const GridPt& pt) const
        {
            SLIC_ASSERT_MSG(hasBlock(pt)
                            ,"(" << pt <<", "<< this->m_level << ") was not a block in the tree at level.");

            // Note: Using find() method on hashmap since operator[] is non-const
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.base());
            return blockIt->second[brood.offset()];
        }

        /** \brief Access the data associated with the entire brood */
        BroodData& getBroodData(const GridPt& pt)
        {
            return m_map[ Brood::MortonizerType::mortonize(pt)];
        }

        /** \brief Const access to data associated with the entire brood */
        const BroodData& getBroodData(const GridPt& pt) const {
            SLIC_ASSERT_MSG(hasBlock(pt)
                            ,"(" << pt <<", "<< this->m_level << ") was not a block in the tree at level.");

            // Note: Using find() method on hashmap since operator[] is non-const
            ConstMapIter blockIt = m_map.find( Brood::MortonizerType::mortonize(pt));
            return blockIt->second;
        }

        /** \brief Predicate to check if there are any blocks in this octree level */
        bool empty() const { return m_map.empty(); }

        /** \brief Returns the number of blocks (internal and leaf) in the level */
        int numBlocks() const
        {
            if(empty())
                return 0;
            return (this->m_level == 0)? 1 : (m_map.size() * Base::BROOD_SIZE);
        }

        /** \brief Returns the number of internal blocks in the level */
        int numInternalBlocks() const
        {
            if(empty())
                return 0;

            const int leaves = numLeafBlocks();
            return (this->m_level == 0 && leaves == 0)
                    ? 1
                    : numBlocks() - numLeafBlocks(); }

        /** \brief Returns the number of leaf blocks in the level */
        int numLeafBlocks() const
        {
            if(empty())
                return 0;

            int count = 0;
            for(ConstMapIter it = m_map.begin(), itEnd = m_map.end(); it != itEnd; ++it)
            {
                const BroodData& bd  = it->second;
                for(int i=0; i< Base::BROOD_SIZE; ++i)
                {
                    if(bd[i].isLeaf())
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
            const Brood brood(pt);
            ConstMapIter blockIt = m_map.find(brood.base());

            return (blockIt == m_map.end())
                    ? BlockNotInTree
                    : (blockIt->second[brood.offset()].isLeaf())
                        ? LeafBlock
                        : InternalBlock;
        }

    //private:
    //  DISABLE_COPY_AND_ASSIGNMENT(OctreeLevel);

    private:
      MapType m_map;
    };




} // end namespace quest

#endif  // MORTON_OCTREE_LEVEL__HXX_
