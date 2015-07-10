/**
 * \file OrderedSet.hpp
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef MESHAPI_ORDERED_SET_H_
#define MESHAPI_ORDERED_SET_H_

#include <cstddef>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include "common/CommonTypes.hpp" // for ATK_NULLPTR
#include "slic/slic.hpp"

#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

namespace asctoolkit {
namespace meshapi {

namespace policies {



    /**
     * \name OrderedSet_Size_Policies
     * \brief A few default policies for the size of an OrderedSet
     */

    /// \{

    /**
     * \brief A policy class for the size of a set whose value can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeSizeHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType();

        RuntimeSizeHolder(IntType sz = DEFAULT_VALUE) : m_sz(sz) {}

        inline IntType size() const { return m_sz;}
        inline IntType& size()            { return m_sz;}

        inline bool isEmpty() const       { return m_sz == IntType(); }
    private:
        IntType m_sz;
    };

    /**
     * \brief A policy class for a compile-time known set size
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeSizeHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeSizeHolder(){}


        CompileTimeSizeHolder(IntType val) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeSizeHolder -- tried to initialize a compile time size policy with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline IntType size() const { return INT_VAL;}
        inline IntType isEmpty() const { return INT_VAL == IntType();}
    };

    /**
      * \brief A policy class for an empty set (no size)
      */
     template<typename IntType>
     struct NoSize {
         static const IntType DEFAULT_VALUE = IntType();

         NoSize(){}


         NoSize(IntType val) {
                SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                               , "MeshAPI::NoSize policy-- tried to initialize a NoSize set with value with value ("
                              << val <<" ) but should always be zero.");
         }

         inline IntType size() const { return DEFAULT_VALUE;}
         inline IntType isEmpty() const { return true;}
     };

    /// \}

    /**
     * \name OrderedSet_Offset_Policies
     * \brief A few default policies for the offset of an OrderedSet
     */

    /// \{


    /**
     * \brief A policy class for the offset in a set.  The offset can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeOffsetHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType();


        RuntimeOffsetHolder(IntType off) : m_off(off) {}

        inline IntType offset() const { return m_off;}
        inline IntType& offset()            { return m_off;}

        inline bool hasOffset() const       { return m_off == IntType(); }
    private:
        IntType m_off;
    };


    /**
     * \brief A policy class for a compile-time known set offset
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeOffsetHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeOffsetHolder(IntType val) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeOffsetHolder -- tried to initialize a compile time offset with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline IntType offset() const { return INT_VAL;}
        inline IntType hasOffset() const { return INT_VAL == IntType();}
    };

    /**
     * \brief A policy class for when we have no offset
     */
    template<typename IntType>
    struct NoOffset {
        static const IntType DEFAULT_VALUE = IntType();

        NoOffset() {}
        NoOffset(IntType val) {
               SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                              , "MeshAPI::NoOffset policy -- tried to initialize a NoOffset policy with ("
                              << val <<", but should always be 0");
        }

        inline IntType offset() const { return DEFAULT_VALUE;}
        inline IntType hasOffset() const { return false;}
    };

    /// \}

    /**
     * \name OrderedSet_Stride_Policies
     * \brief A few default policies for the stride of an OrderedSet
     */

    /// \{


    /**
     * \brief A policy class for the stride in a set.  The stride can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeStrideHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType(1);

        RuntimeStrideHolder(IntType stride = DEFAULT_VALUE) : m_stride(stride) {}

        inline IntType stride() const { return m_stride;}
        inline IntType& stride()            { return m_stride;}

        inline bool hasStride() const       { return m_stride != IntType(); }
    private:
        IntType m_stride;
    };


    /**
     * \brief A policy class for a compile-time known stride
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeStrideHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeStrideHolder() {}
        CompileTimeStrideHolder(IntType val) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeStrideHolder -- tried to initialize a compile time stride with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline IntType stride() const { return INT_VAL;}
        inline IntType hasStride() const { return INT_VAL != IntType();}
    };

    /**
     * \brief A policy class for a set with stride one (i.e. the default stride)
     */
    template<typename IntType>
    struct StrideOne {
        static const IntType DEFAULT_VALUE = IntType(1);

        StrideOne() {}
        /**
         * This constructor only exists to allow the derived class to not have to specialize for when the stride is known at compile time
         */
        StrideOne(IntType val) {
               SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                              , "MeshAPI::StrideOne policy -- tried to initialize a stride-one StridePolicy with value ("
                              << val <<"), but should always be 1.");
        }

        inline IntType stride() const { return DEFAULT_VALUE;}
        inline IntType hasStride() const { return DEFAULT_VALUE != IntType();}
    };



    /// \}

    /**
     * \name OrderedSet_Indirection_Policies
     * \brief A few default policies for the indirection of an OrderedSet
     */

    /// \{


    /**
     * \brief A policy class for sets with no indirection
     */
    template<typename PositionType, typename ElementType>
    struct NoIndirection
    {
        inline ElementType indirection(PositionType pos) const { return static_cast<ElementType>(pos); }
        bool hasIndirection() const { return false;}
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct ArrayIndirection
    {
        ArrayIndirection(ElementType* buf = ATK_NULLPTR) : m_arrBuf(buf) {}

        ElementType*& data() { return m_arrBuf;}

        inline ElementType indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection()
                           , "MeshAPI::Set:ArrayIndirection -- Tried to dereference a null array in an array based indirection set.");
            return m_arrBuf[pos];
        }
        bool hasIndirection() const { return m_arrBuf == ATK_NULLPTR;}
    private:
        ElementType* m_arrBuf;
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct STLVectorIndirection
    {
        typedef std::vector<ElementType> VectorType;

        STLVectorIndirection(VectorType* buf = ATK_NULLPTR) : m_vecBuf(buf) {}

        VectorType& data() { return *m_vecBuf;}

        inline ElementType indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection(), "MeshAPI::Set:STLVectorIndirection -- Tried to dereference a null vector in a vector based indirection set.");
            SLIC_ASSERT_MSG( pos < m_vecBuf->size(), "MeshAPI::Set:STLVectorIndirection -- Tried to access an out of bounds element at position "
                    << pos << " in vector with only " << m_vecBuf->size() << " elements.");

            return (*m_vecBuf)[pos];
        }
        bool hasIndirection() const { return m_vecBuf == ATK_NULLPTR;}
    private:
        VectorType* m_vecBuf;
    };

    /// \}

    /**
     * \name OrderedSet_Subsetting_Policies
     * \brief A few default policies for the subsetting of an OrderedSet
     */

    /// \{

    struct NoSubset
    {
        static const NullSet s_nullSet;

        bool isSubset() const { return false; }
        const Set* parentSet() const { return &s_nullSet; }
    };

    struct VirtualParentSubset
    {
        static NullSet s_nullSet;

        VirtualParentSubset() : m_parentSet(&s_nullSet) {}

        bool isSubset() const { return *m_parentSet != s_nullSet; }
        const Set* parentSet() const { return m_parentSet; }
              Set*& parentSet()       { return m_parentSet; }

    private:
        Set* m_parentSet;
    };

    template<typename ParentSetType>
    struct ConcreteParentSubset
    {

        bool isSubset() const { return m_parentSet == ATK_NULLPTR; }
        const ParentSetType* parentSet() const { return m_parentSet; }
              ParentSetType*& parentSet()       { return m_parentSet; }

    private:
        ParentSetType* m_parentSet;
    };


    /// \}
}




/**
 * \class OrderedSet
 * \brief Models a set whose elements can be defined as a strided offsets of the position, possibly with a level of indirection.
 * \details Specifically, the element at position pos can be defined as:  static_cast<ElementType>( indirection[ pos * stride + offset ] )
 */
    template< typename SizePolicy           = policies::RuntimeSizeHolder<Set::PositionType>
            , typename OffsetPolicy         = policies::NoOffset<Set::PositionType>
            , typename StridePolicy         = policies::StrideOne<Set::PositionType>
            , typename IndirectionPolicy    = policies::NoIndirection<Set::PositionType, Set::ElementType>
            , typename SubsettingPolicy     = policies::NoSubset
            >
    struct OrderedSet: public Set, SizePolicy, OffsetPolicy, StridePolicy, IndirectionPolicy, SubsettingPolicy
  {
  public:

    typedef Set::IndexType                        IndexType;
    typedef Set::PositionType                     PositionType;
    typedef IndexType                             ElementType;

    typedef SizePolicy SizePolicyType;
    typedef OffsetPolicy OffsetPolicyType;
    typedef StridePolicy StridePolicyType;
    typedef IndirectionPolicy IndirectionPolicyType;

    typedef boost::counting_iterator<ElementType> iterator;
    typedef std::pair<iterator,iterator>          iterator_pair;

  public:
    OrderedSet(PositionType size    = SizePolicy::DefaultValue()
             , PositionType offset  = OffsetPolicy::DefaultValue()
             , PositionType stride  = StridePolicy::DefaultValue()
             // Note: constructor does not yet take an indirection type pointer...
             //, const Set* parentSet = &s_nullSet
            )
        : SizePolicy(size)
        , OffsetPolicy(offset)
        , StridePolicy(stride)
        //, SubsettingPolicy(parentSet)
      {}

public:
    /**
     * \brief Given a position in the Set, return a position in the larger index space
     */
    inline ElementType operator[](PositionType pos) const
    {
            verifyPosition(pos);
            return IndirectionPolicy::indirection( pos * StridePolicy::stride() + OffsetPolicy::offset() );
    }

    inline ElementType         at(PositionType pos)         const { return operator[](pos); }


    inline PositionType        size()  const { return SizePolicy::size(); }
    inline bool                isEmpty() const { return SizePolicy::isEmpty(); }

  public:   // Functions related to iteration
    iterator            begin() const { return iterator( OffsetPolicy::offset() ); }
    iterator            end()   const { return iterator( size() * StridePolicy::stride() + OffsetPolicy::offset() ); }
    iterator_pair       range() const { return std::make_pair(begin(), end()); }

    /// HACK: This function needs to be implemented
    bool                isValid(bool verboseOutput = false) const
    {
        if(verboseOutput)
            return true;
        return true;
    }

    bool isSubset() const { return SubsettingPolicy::isSubset(); }

  private:

    inline void         verifyPosition(PositionType pos)       const
    {
      SLIC_ASSERT_MSG( pos >= 0 && pos < size()
          , "MeshAPI::OrderedSet -- requested out-of-range element at position "
          << pos << ", but set only has " << size() << " elements." );
    }

  private:
    /// NOTE: All data for OrderedSet is associated with parent policy classes
  };


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_ORDERED_SET_H_
