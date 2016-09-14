
#ifndef OCTREE_LEVEL__HXX_
#define OCTREE_LEVEL__HXX_


#include "common/config.hpp"
#include "common/CommonTypes.hpp"

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
 * OctreeLevel is an abstract base class.
 * This file also defines two concrete instantiations:
 * * GridPointOctreeLevel uses a GridPoint as a hash table key for its octree blocks
 * * MortonOctreeLevel uses a Morton index (of the given bit width) as a hash key
 *   for its octree blocks.
 */

namespace quest
{

    /**
     * \brief Helper enumeration for status of a BlockIndex within an OctreeLevel instance
     */
    enum TreeBlockStatus {
          BlockNotInTree   /** Status of blocks that are not in the tree */
        , LeafBlock        /** Status of blocks that are leaves in the tree */
        , InternalBlock    /** Status of blocks that are internal to the tree */
    };

    /**
     * \class
     * \brief An abstract base class to represent a sparse level of blocks within an octree.
     * Each block is associated with an integer grid point whose coordinates
     * have values between 0 and 2^L (where L = this->level() is the encoded level).
     * The OctreeLevel associates data of (templated) type BlockDataType with each such block.
     * \note For efficiency, the data is stored within a brood (collection of siblings that are created simultaneously).
     * \note BlockDataType must define a predicate function with the signature: bool isLeaf() const;
     */
    template<int DIM, typename BlockDataType>
    class OctreeLevel
    {
    public:
        /** The coordinate type of a block in the octree */
        typedef int CoordType;

        /**
         * \brief A type for the grid points of the octree.
         * \note CoordType must be an integral type
         */
        typedef quest::Point<CoordType,DIM> GridPt;

        enum { BROOD_SIZE = 1 << DIM };

        /** A brood is a collection of sibling blocks that are generated simultaneously */
        typedef quest::NumericArray< BlockDataType, BROOD_SIZE> BroodData;

        /** Predeclare the BlockIterator type */
        template<typename OctreeLevel, typename IterHelper, typename DataType> class BlockIterator;

    protected:

      /**
       * \brief A virtual base class to help with iteration of an OctreeLevel's blocks
       */
      class BlockIteratorHelper
      {
      public:
          /** Virtual destructor */
          virtual ~BlockIteratorHelper() {}

          /** \brief A function to increment to the next Block in the level */
          virtual void increment() = 0;

          /** \brief Predicate to determine if two BlockIteratorHelpers are equivalent */
          virtual bool equal(const BlockIteratorHelper* other) = 0;

          /** Accessor for the point associated with the current octree block */
          virtual GridPt pt() const =0;

          /** Accessor for the data associated with the current octree block */
          virtual BlockDataType* data() = 0;

          /** Const accessor for the data associated with the current octree block */
          virtual const BlockDataType* data() const = 0;
      };

      /**
       * \brief A virtual base class to help with constant iteration of an OctreeLevel's blocks
       */
      class ConstBlockIteratorHelper
      {
      public:
          /** Virtual desctructor */
          virtual ~ConstBlockIteratorHelper() {}

          /** \brief A function to increment to the next Block in the level */
          virtual void increment() = 0;

          /** \brief Predicate to determine if two BlockIteratorHelpers are equivalent */
          virtual bool equal(const ConstBlockIteratorHelper* other) = 0;

          /** Accessor for the point associated with the current octree block */
          virtual GridPt pt() const =0;

          /** Const accessor for the data associated with the current octree block */
          virtual const BlockDataType* data() const = 0;
      };

    public:
      typedef BlockIterator<OctreeLevel, BlockIteratorHelper, BlockDataType> BlockIter;
      typedef BlockIterator<const OctreeLevel, ConstBlockIteratorHelper, const BlockDataType> ConstBlockIter;

    public:
      /** \brief Constructor of an OctreeLevel at level lev */
      OctreeLevel(int level = -1): m_level(level){}

      /** \brief Virtual destructor of an OctreeLevel */
      virtual ~OctreeLevel() {};

      /**
       * \brief Returns the maximum coordinate value in the level
       * \note This is (2^l -1), where L is the current level
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

      /** Accessor for the instance's level */
      int level() const { return m_level; }

      /**
       * \brief Predicate to check whether the block associated with the given GridPt pt is an allowed block in the level
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
       * \note Uses a helper class to manage the polymorphic OctreeLevel's iteration
       *       The helper defines the following functions: increment(), equal() update(), pt() and data()
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
            m_iterHelper = octLevel->getIteratorHelper(begin); // factory function
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
         * \brief A const dereference function used for accessing the data
         * \note Only valid for constant access to data associated with an octree block
         * \note For non-const access on a non-const iterator, use the data() function
         * \sa data()
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
        IterHelper* m_iterHelper;             /** Instance of iterator helper class */

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

      /** \brief Virtual function to check the status of a block (e.g. Leaf, Internal, NotInTree) */
      virtual TreeBlockStatus blockStatus(const GridPt & pt) const = 0;

      /** \brief Virtual predicate to determine if the OctreeLevel is empty */
      virtual bool empty() const =0;

      /** \brief Virtual predicate to determine if the OctreeLevel has a block with the given grid point pt */
      virtual bool hasBlock(const GridPt& pt) const =0;

      /** \brief Virtual function to add all children of the given grid point pt to the OctreeLevel */
      virtual void addAllChildren(const GridPt& pt) = 0;

      /** \brief Virtual const accessor for the data associated with grid point pt */
      virtual const BlockDataType& operator[](const GridPt& pt) const = 0;
      /** \brief Virtual accessor for the data associated with grid point pt */
      virtual       BlockDataType& operator[](const GridPt& pt)       = 0;

      /** \brief Virtual accessor for the data associated with all children of the given grid point (i.e. the brood) */
      virtual BroodData& getBroodData(const GridPt& pt) =0;
      /** \brief Virtual const accessor for the data associated with all children of the given grid point (i.e. the brood) */
      virtual const BroodData& getBroodData(const GridPt& pt) const =0;

      /**
       *  \brief Virtual factory function to create an iterator helper
       *  \param A boolean to determine if the iterator should be a begin iterator (true) or an end iterator (false)
       */
      virtual BlockIteratorHelper* getIteratorHelper(bool) = 0;
      /**
       * \brief Virtual factory function to create a const iterator helper
       * \param A boolean to determine if the iterator should be a begin iterator (true) or an end iterator (false)
       */
      virtual ConstBlockIteratorHelper* getIteratorHelper(bool) const = 0;


      /** \brief Virtual function to compute the number of blocks (internal and leaf) in the level */
      virtual int numBlocks() const =0;

      /** \brief Virtual function to compute the number of internal blocks in the level */
      virtual int numInternalBlocks() const =0;

      /** \brief Virtual function to compute the number of leaf blocks in the level */
      virtual int numLeafBlocks() const =0;


    protected:
      int m_level;
    };

    /**
     * \class
     * \brief A GridPointOctreeLevel is a representation of a sparse OctreeLevel.
     *  It is a concrete implementation of an OctreeLevel
     *
     *  It associates data with its Octree block using a hash map
     *  whose key type is an integer GridPoint (using a Morton-based hash function),
     *  and whose value type is a BlockDataType.
     *  For efficiency, the data is associated with an entire brood, a collection of
     *  siblings that are created simultaneously. In dimension DIM, there are 2^DIM siblings in a brood.
     *
     *  \see OctreeLevel
     */
    template<int DIM, typename BlockDataType>
    class GridPointOctreeLevel : public OctreeLevel<DIM,BlockDataType>
    {
    public:
      typedef OctreeLevel<DIM, BlockDataType> Base;
      typedef typename Base::GridPt                   GridPt;
      typedef typename Base::BroodData                BroodData;
      typedef typename Base::BlockIteratorHelper      BaseBlockIteratorHelper;
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
         * \brief Private helper class to handle subindexing of block data within octree siblings
         * \note A brood is a collection of siblings that are generated simultaneously.
         * \note This class converts a grid point at the given level into a brood index of the point.
         *       The base brood point has the coordinates of the grid point's octree parent
         *       and its offset index is obtained by interleaving the least significant bit of its coordinates
         *       in each dimension.
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
          int m_offset;      /** Index of the block within the brood. Value is in [0, 2^DIM) */
      };

    public:

        /**
         * \brief Concrete instance of the BlockIteratorHelper class defined in the OctreeLevel base class.
         */
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

            /** \brief Advances the iterator to the next block */
            void increment()
            {
                ++m_offset;

                if(m_offset == Base::BROOD_SIZE)
                {
                    ++m_currentIter;
                    m_offset = 0;
                }
            }

            /** \brief Accessor for the iterator's current point */
            GridPt pt() const
            {
                // Reconstruct the grid point from its brood representation
                GridPt itPt = m_currentIter->first;
                for(int i=0; i<DIM; ++i)
                    itPt[i] = (itPt[i]<<1) + ( m_offset & (1 << i)? 1 : 0);

                return itPt;
            }

            /** \brief Accessor for the iterator's current data */
            BlockDataType* data() { return &m_currentIter->second[m_offset]; }
            /** \brief Const accessor for the iterator's current data */
            const BlockDataType* data() const { return &m_currentIter->second[m_offset]; }

            /** \brief Predicate to determine if two BlockIterators are the same */
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

        /** \brief Default constructor for an octree level */
        GridPointOctreeLevel(int level = -1): Base(level){}


        /**
         * \brief Factory function to return a GridBlockIterHelper for this level
         *  \param begin A boolean to determine if this is to be a begin (true) or end (false) iterator
         */
        BaseBlockIteratorHelper* getIteratorHelper(bool begin)
        {
            return new GridBlockIterHelper(this, begin);
        }

        /**
         * \brief Factory function to return a ConstGridBlockIterHelper for this level
         *  \param begin A boolean to determine if this is to be a begin (true) or end (false) iterator
         */
        ConstBaseBlockIteratorHelper* getIteratorHelper(bool begin) const
        {
            return new ConstGridBlockIterHelper(this, begin);
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

        /** \brief Access the data associated with the entire brood */
        BroodData& getBroodData(const GridPt& pt) { return m_map[ pt]; }

        /** \brief Const access to data associated with the entire brood */
        const BroodData& getBroodData(const GridPt& pt) const {
            SLIC_ASSERT_MSG(hasBlock(pt)
                            ,"(" << pt <<", "<< this->m_level << ") was not a block in the tree at level.");

            // Note: Using find() method on hashmap since operator[] is non-const
            ConstMapIter blockIt = m_map.find( pt);
            return blockIt->second;
        }


        /** \brief Predicate to check if there are any blocks in this octree level */
        bool empty() const { return m_map.empty(); }

        /** \brief Returns the number of blocks (internal and leaf) in the level */
        int numBlocks() const { return m_map.size() * Base::BROOD_SIZE; }

        /** \brief Returns the number of internal blocks in the level */
        int numInternalBlocks() const { return numBlocks() - numLeafBlocks(); }

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
        int numBlocks() const { return m_map.size() * Base::BROOD_SIZE; }

        /** \brief Returns the number of internal blocks in the level */
        int numInternalBlocks() const { return numBlocks() - numLeafBlocks(); }

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

#endif  // OCTREE_LEVEL__HXX_
