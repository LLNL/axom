
#ifndef OCTREE_LEVEL__HXX_
#define OCTREE_LEVEL__HXX_


#include "common/config.hpp"

#include "fmt/fmt.hpp"
#include "slic/slic.hpp"

#include "quest/MortonIndex.hpp"
#include "quest/NumericArray.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"


#ifdef USE_CXX11
#include <unordered_map>
#else
#include "boost/unordered_map.hpp"
#endif

#include <boost/iterator/iterator_facade.hpp>

/**
 * \file
 * \brief Defines templated OctreeLevel class
 * An OctreeLevel associates data with the integer points on a sparse grid.
 */

namespace quest
{

    /**
     * \brief Helper enumeration for status of a BlockIndex within an OctreeLevel instance
     */
    enum TreeBlockStatus {
          BlockNotInTree   // Status of blocks that are not in the tree
        , LeafBlock        // Status of blocks that are leaves in the tree
        , InternalBlock    // Status of blocks that are internal to the tree
    };

    /**
     * \class
     * \brief A class to represent a sparse level of blocks within an octree.
     * Each block is associated with an integer grid point whose coordinates
     * have values between 0 and 2^L (where L = this->level() is the encoded level).
     * The OctreeLevel associates data of (templated) type BlockDataType with each such block.
     * \note BlockDataType must define a predicate function with the signature: bool isLeaf() const;
     */
    template<int DIM, typename BlockDataType>
    class OctreeLevel
    {
    public:
        typedef int CoordType;
        typedef quest::Point<CoordType,DIM> GridPt;
        typedef quest::Vector<CoordType,DIM> GridVec;

        enum { BROOD_SIZE = 1 << DIM };

    public:

      // A brood is a collection of sibling blocks that are generated simultaneously
      typedef quest::NumericArray< BlockDataType, BROOD_SIZE> BroodData;

      template<typename OctreeLevel, typename IterHelper, typename DataType> class BlockIterator;

    protected:

      class BlockIteratorHelper
      {
      public:
          virtual ~BlockIteratorHelper() {}
          virtual void increment() = 0;
          virtual bool equal(const BlockIteratorHelper* other) = 0;
          virtual GridPt pt() const =0;
          virtual BlockDataType* data() = 0;
          virtual const BlockDataType* data() const = 0;
      };

      class ConstBlockIteratorHelper
      {
      public:
          virtual ~ConstBlockIteratorHelper() {}
          virtual void increment() = 0;
          virtual bool equal(const ConstBlockIteratorHelper* other) = 0;
          virtual GridPt pt() const =0;
          virtual const BlockDataType* data() const = 0;
      };

    public:
      typedef BlockIterator<OctreeLevel, BlockIteratorHelper, BlockDataType> BlockIter;
      typedef BlockIterator<const OctreeLevel, ConstBlockIteratorHelper, const BlockDataType> ConstBlockIter;

    public:
      OctreeLevel(int level = -1): m_level(level){}

      virtual ~OctreeLevel() {};

      /**
       * \brief Returns the maximum coordinate value in the level
       * \note This is (2^L -1), where L is the current level
       */
      CoordType maxCoord() const
      {
          return (1<< m_level) -1;
      }

      /**
       * \brief Returns a GridPt whose coordinates are set to maxCoord
       * \sa maxCoord()
       */
      GridPt maxGridCell() const
      {
          return GridPt(maxCoord());
      }

      int level() const { return m_level; }

      /**
       * \brief Predicate to check whether the block associated with the given GridPt pt is an allowed block ih the level
       * \param [in] pt The gridpoint of the block to check
       * \note pt is inBounds if each of its coordinates is a non-negative integer less than maxCoord()
       * \sa maxCoord()
       */
      bool inBounds(const GridPt& pt) const
      {
          const CoordType maxVal = maxCoord();
          for(int i=0; i< DIM; ++i)
              if( pt[i] < 0 || pt[i] > maxVal)
                  return false;
          return true;
      }


    public:

      /**
       * \class
       * \brief An iterator type for the blocks of an octree level
       */
      template<typename OctreeLevel, typename IterHelper, typename DataType>
      class BlockIterator : public boost::iterator_facade< BlockIterator<OctreeLevel, IterHelper, DataType>
                                 , DataType
                                 , boost::forward_traversal_tag
                                 , DataType
                                 >
      {
      public:
        typedef BlockIterator<OctreeLevel, IterHelper, DataType>  iter;

        BlockIterator(OctreeLevel* octLevel, bool begin = false)
            : m_octLevel(octLevel)
        {
            SLIC_ASSERT(octLevel != ATK_NULLPTR);
            m_iterHelper = octLevel->getIteratorHelper(octLevel,begin);
        }

        ~BlockIterator()
        {
            if(m_iterHelper != ATK_NULLPTR)
            {
                delete m_iterHelper;
                m_iterHelper = ATK_NULLPTR;
            }
        }


        /**
         * \brief A const dereference function used for
         * \note Only valid for constant access to data associated with an octree block
         * \note For non-const access on a non-const accessor, use the data() function
         */
        const DataType& dereference() const { return *m_iterHelper->data(); }

        /**
         * \brief Const accessor for the iterator's current grid point
         */
        GridPt pt() const { return m_iterHelper->pt(); }

        /**
         * \brief Non-const accessor for data associated with the iterator's current grid point
         */
        DataType& data() { return *m_iterHelper->data(); }

        /**
         * \brief Const accessor for data associated with the iterator's current grid point
         */
        const DataType& data() const { return *m_iterHelper->data(); }

        /**
         * \brief Equality test against another iterator
         * \param other The other iterator
         * \return true, if the two iterators are equal, false otherwise
         */
        bool equal(const iter& other) const
        {
            return (m_octLevel == other.m_octLevel)       // point to same object
                    && m_iterHelper->equal(other.m_iterHelper);
        }

        /** \brief Increment the iterator to the next point */
        void increment()
        {
            m_iterHelper->increment();
        }

      private:
        friend class boost::iterator_core_access;
        OctreeLevel*  m_octLevel;             /** Pointer to the iterator's container class */
        IterHelper* m_iterHelper;

      };


    public:

      /**
       * \brief Predicate to check whether the block associated with the given GridPt pt is a leaf block
       */
      bool isLeaf(const GridPt& pt) const { return this->blockStatus(pt) == LeafBlock; }

      /**
       * \brief Predicate to check whether the block associated with the given GridPt pt is an internal block
       */
      bool isInternal(const GridPt& pt) const { return this->blockStatus(pt) == InternalBlock; }


      /** \brief Begin iterator to points and data in tree level */
      BlockIter      begin()       { return BlockIter(this,true); }

      /** \brief Const begin iterator to points and data in tree level */
      ConstBlockIter begin() const { return ConstBlockIter(this,true); }

      /** \brief End iterator to points and data in tree level */
      BlockIter      end()         { return BlockIter(this,false); }

      /** \brief Const end iterator to points and data in tree level */
      ConstBlockIter end()   const { return ConstBlockIter(this,false); }

      virtual TreeBlockStatus blockStatus(const GridPt & pt) const = 0;

      virtual bool empty() const =0;
      virtual bool hasBlock(const GridPt& pt) const =0;
      virtual void addAllChildren(const GridPt& pt) = 0;

      virtual const BlockDataType& operator[](const GridPt& pt) const = 0;
      virtual       BlockDataType& operator[](const GridPt& pt)       = 0;

      virtual BlockIteratorHelper* getIteratorHelper(OctreeLevel*, bool) = 0;
      virtual ConstBlockIteratorHelper* getIteratorHelper(const OctreeLevel*, bool) const = 0;


    protected:
      int m_level;
    };

    template<int DIM, typename BlockDataType>
    class GridPointOctreeLevel : public OctreeLevel<DIM,BlockDataType>
    {
    public:
      typedef OctreeLevel<DIM, BlockDataType> Base;
      typedef typename Base::GridPt GridPt;
      typedef typename Base::BroodData BroodData;
      typedef typename Base::BlockIteratorHelper BaseBlockIteratorHelper;
      typedef typename Base::ConstBlockIteratorHelper ConstBaseBlockIteratorHelper;

      #if defined(USE_CXX11)
        typedef std::unordered_map<GridPt, BroodData, PointHash<int> > MapType;
      #else
        typedef boost::unordered_map<GridPt, BroodData, PointHash<int> > MapType;
      #endif

        typedef typename MapType::iterator       MapIter;
        typedef typename MapType::const_iterator ConstMapIter;

        template<typename OctreeLevelType, typename InnerIterType, typename ParentType> class GridBlockIteratorHelper;

        typedef GridBlockIteratorHelper<GridPointOctreeLevel, MapIter, BaseBlockIteratorHelper> GridBlockIterHelper;
        typedef GridBlockIteratorHelper<const GridPointOctreeLevel, ConstMapIter, ConstBaseBlockIteratorHelper> ConstGridBlockIterHelper;

    protected:
        /**
         * \class
         * \brief Private inner class to handle subindexing of block data within octree siblings
         * \note A brood is a collection of siblings that are generated simultaneously.
         * \note This class converts a grid point at the given level into a brood index of the point.
         *       The base brood point is that of the grid point's octree parent
         *       and its offset index is obtained by interleaving the least significant bit of its coordinates.
         */
      struct Brood {

          /**
           * \brief Constructor for a brood offset relative to the given grid point pt
           * \param [in] pt The grid point within the octree level
           */
          Brood(const GridPt& pt)
              : m_broodPt( pt.array() /2), m_offset(0)
          {
              for(int i=0; i< DIM; ++i)
              {
                  m_offset |= (pt[i]& 1) << i; // interleave the least significant bits
              }
          }

          /** \brief Accessor for the base point of the entire brood */
          const GridPt& base() const { return m_broodPt; }

          /** \brief Accessor for the index of the point within the brood */
          const int& offset() const { return m_offset; }

      private:
          GridPt m_broodPt;  /** Base point of all blocks within the brood */
          int m_offset;         /** Index of the block within the brood. Value is in [0, 2^DIM) */
      };

    public:

        template<typename OctreeLevelType, typename InnerIterType, typename ParentType>
        class GridBlockIteratorHelper : public ParentType
        {
        public:
            typedef GridBlockIteratorHelper<OctreeLevelType, InnerIterType, ParentType> GridBlockItType;
            typedef ParentType     BaseBlockItType;

            GridBlockIteratorHelper(OctreeLevelType* octLevel, bool begin)
                : m_endIter( octLevel->m_map.end() )
                , m_offset(0)
            {
                m_currentIter = begin ? octLevel->m_map.begin() : m_endIter;
            }

            void increment()
            {
                ++m_offset;

                if(m_offset == Base::BROOD_SIZE)
                {
                    ++m_currentIter;
                    m_offset = 0;
                }
            }

            GridPt pt() const
            {
                // Reconstruct the grid point from its brood representation
                GridPt itPt = m_currentIter->first;
                for(int i=0; i<DIM; ++i)
                    itPt[i] = (itPt[i]<<1) + ( m_offset & (1 << i)? 1 : 0);

                return itPt;
            }

            BlockDataType* data() { return &m_currentIter->second[m_offset]; }
            const BlockDataType* data() const { return &m_currentIter->second[m_offset]; }

            bool equal(const BaseBlockItType* other)
            {
                const GridBlockItType* pother = dynamic_cast<const GridBlockItType*>(other);

                return (pother != ATK_NULLPTR)
                     && (m_currentIter == pother->m_currentIter)   // iterators are the same
                     && (m_offset == pother->m_offset);               // brood indices are the same
            }
        private:
            InnerIterType m_currentIter, m_endIter;
            int m_offset;
        };

    public:

        /**
         * \brief Default constructor for an octree level
         */
        GridPointOctreeLevel(int level = -1): Base(level){}


        BaseBlockIteratorHelper* getIteratorHelper(Base* octLevel, bool begin)
        {
            return new GridBlockIterHelper(static_cast<GridPointOctreeLevel*>(octLevel), begin);
        }

        ConstBaseBlockIteratorHelper* getIteratorHelper(const Base* octLevel, bool begin) const
        {
            return new ConstGridBlockIterHelper(static_cast<const GridPointOctreeLevel*>(octLevel), begin);
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

            m_map[pt];  // Adds children, if not already present, using default BlockDataType() constructor
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



        /** \brief Predicate to check if there are any blocks in this octree level */
        bool empty() const { return m_map.empty(); }

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


    template<int DIM, typename BlockDataType>
    class MortonOctreeLevel : public OctreeLevel<DIM,BlockDataType>
    {
    public:
      typedef OctreeLevel<DIM, BlockDataType> Base;
      typedef typename Base::GridPt GridPt;
      typedef typename Base::BroodData BroodData;
      typedef typename Base::BlockIteratorHelper BaseBlockIteratorHelper;
      typedef typename Base::ConstBlockIteratorHelper ConstBaseBlockIteratorHelper;

      #if defined(USE_CXX11)
          typedef std::unordered_map<MortonIndex, BroodData> MapType;
    #else
          typedef boost::unordered_map<MortonIndex, BroodData> MapType;
      #endif

        typedef typename MapType::iterator       MapIter;
        typedef typename MapType::const_iterator ConstMapIter;

        template<typename OctreeLevelType, typename InnerIterType, typename ParentType> class MortonBlockIteratorHelper;

        typedef MortonBlockIteratorHelper<MortonOctreeLevel, MapIter, BaseBlockIteratorHelper> MortonBlockIterHelper;
        typedef MortonBlockIteratorHelper<const MortonOctreeLevel, ConstMapIter, ConstBaseBlockIteratorHelper> ConstMortonBlockIterHelper;

    protected:
        /**
         * \class
         * \brief Private inner class to handle subindexing of block data within octree siblings
         * \note A brood is a collection of siblings that are generated simultaneously.
         * \note This class converts a grid point at the given level into a brood index of the point.
         *       The base brood point is that of the grid point's octree parent
         *       and its offset index is obtained by interleaving the least significant bit of its coordinates.
         */
      struct Brood {
          enum { BROOD_BITMASK = Base::BROOD_SIZE -1 };

          typedef Mortonizer<typename GridPt::CoordType, GridPt::NDIMS> MortonizerType;

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
          const MortonIndex& base() const { return m_broodIdx; }

          /** \brief Accessor for the index of the point within the brood */
          const int& offset() const { return m_offset; }

      private:
          MortonIndex m_broodIdx;  /** Base point of all blocks within the brood */
          int m_offset;         /** Index of the block within the brood. Value is in [0, 2^DIM) */
      };

    public:

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

            void increment()
            {
                ++m_offset;

                if(m_offset == Base::BROOD_SIZE || m_isLevelZero)
                {
                    ++m_currentIter;
                    m_offset = 0;
                }
            }

            GridPt pt() const
            {
                // Reconstruct the grid point from its brood representation
                typedef Mortonizer<typename GridPt::CoordType, GridPt::NDIMS> Mort;
                return Mort::demortonize( (m_currentIter->first << DIM)  + m_offset);
            }

            BlockDataType* data() { return &m_currentIter->second[m_offset]; }
            const BlockDataType* data() const { return &m_currentIter->second[m_offset]; }

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

        /**
         * \brief Default constructor for an octree level
         */
        MortonOctreeLevel(int level = -1): Base(level){}


        BaseBlockIteratorHelper* getIteratorHelper(Base* octLevel, bool begin)
        {
            return new MortonBlockIterHelper(static_cast<MortonOctreeLevel*>(octLevel), begin);
        }

        ConstBaseBlockIteratorHelper* getIteratorHelper(const Base* octLevel, bool begin) const
        {
            return new ConstMortonBlockIterHelper(static_cast<const MortonOctreeLevel*>(octLevel), begin);
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



        /** \brief Predicate to check if there are any blocks in this octree level */
        bool empty() const { return m_map.empty(); }

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

#endif  // OCTREE_LEVEL__HXX_
