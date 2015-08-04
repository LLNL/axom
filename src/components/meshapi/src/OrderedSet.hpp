/**
 * \file
 *
 * \brief Basic API for an ordered set of entities in a simulation
 * \note We are actually storing (ordered) multisets, since elements can be repeated an arbitrary number of times (e.g. for indirection sets)
 *
 */

#ifndef MESHAPI_ORDERED_SET_H_
#define MESHAPI_ORDERED_SET_H_

#include <cstddef>
#include <vector>


#ifndef MESHAPI_USE_COUNTING_ITERATOR
//    #define MESHAPI_USE_COUNTING_ITERATOR
#endif

#ifdef MESHAPI_USE_COUNTING_ITERATOR
    #include <boost/iterator/counting_iterator.hpp>
#else
    #include <boost/iterator/iterator_facade.hpp>
    #include <boost/utility/enable_if.hpp>
    #include <boost/type_traits.hpp>
#endif

#include "common/CommonTypes.hpp" // for ATK_NULLPTR
#include "slic/slic.hpp"

#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

#include "meshapi/SizePolicies.hpp"
#include "meshapi/OffsetPolicies.hpp"
#include "meshapi/StridePolicies.hpp"
#include "meshapi/IndirectionPolicies.hpp"
#include "meshapi/SubsettingPolicies.hpp"


namespace asctoolkit {
namespace meshapi {


/**
 * \class
 * \brief Models a set whose elements can be defined as a strided offsets of the position, possibly with a level of indirection.
 * \details Specifically, the element at position pos can be defined as:  static_cast<ElementType>( indirection[ pos * stride + offset ] )
 */
    template< typename SizePolicy          = policies::RuntimeSizeHolder<Set::PositionType>
            , typename OffsetPolicy        = policies::ZeroOffset<Set::PositionType>
            , typename StridePolicy        = policies::StrideOne<Set::PositionType>
            , typename IndirectionPolicy   = policies::NoIndirection<Set::PositionType, Set::ElementType>
            , typename SubsettingPolicy    = policies::NoSubset
            >
    struct OrderedSet: public Set
                            , SizePolicy
                            , OffsetPolicy
                            , StridePolicy
                            , IndirectionPolicy
                            , SubsettingPolicy
  {
  public:

    typedef Set::IndexType                              IndexType;
    typedef Set::PositionType                           PositionType;
    typedef IndexType                                   ElementType;

    typedef SizePolicy          SizePolicyType;
    typedef OffsetPolicy        OffsetPolicyType;
    typedef StridePolicy        StridePolicyType;
    typedef IndirectionPolicy   IndirectionPolicyType;
    typedef SubsettingPolicy    SubsettingPolicyType;

#ifdef MESHAPI_USE_COUNTING_ITERATOR
    typedef boost::counting_iterator<ElementType> iterator;
    typedef std::pair<iterator,iterator>          iterator_pair;

    typedef const iterator const_iterator;
    typedef std::pair<const_iterator,const_iterator>  const_iterator_pair;
#else
    template<typename OrderedSetType> class OrderedSetIterator;

    typedef OrderedSetIterator<const OrderedSet> const_iterator;
    typedef std::pair<const_iterator,const_iterator>          const_iterator_pair;

    typedef const_iterator iterator;
    typedef const_iterator_pair iterator_pair;
#endif

  public:
    OrderedSet(PositionType size    = SizePolicyType::DEFAULT_VALUE
             , PositionType offset  = OffsetPolicyType::DEFAULT_VALUE
             , PositionType stride  = StridePolicyType::DEFAULT_VALUE
             // Note: constructor does not yet take an indirection type pointer...
             //, const Set* parentSet = &s_nullSet
            )
        : SizePolicyType(size)
        , OffsetPolicyType(offset)
        , StridePolicyType(stride)
        //, SubsettingPolicyType(parentSet)
      {}


public:
    // define an iterator type -- w/ stride and indirection
    template<typename OrderedSet>
    class OrderedSetIterator : public boost::iterator_facade< OrderedSetIterator<OrderedSet>
                                                            , typename OrderedSet::ElementType
                                                            , std::random_access_iterator_tag
                                                            , typename OrderedSet::ElementType
                                                            , typename OrderedSet::PositionType
                                                            >
    {
    public:
        typedef OrderedSetIterator<OrderedSet> iter;
        typedef typename OrderedSet::ElementType ElementType;
        typedef typename OrderedSet::PositionType PositionType;

        typedef typename OrderedSet::IndirectionPolicyType IndirectionType;
        typedef typename OrderedSet::StridePolicyType StrideType;
    public:
        OrderedSetIterator(PositionType pos, const OrderedSet* oSet)
            : m_pos(pos), m_orderedSet(oSet) {}


        const ElementType & dereference()    const {
            // Note: Since we return a reference to the pointed-to value, we need different functions
            //       for OrderedSets with indirection buffers than with those that have no indirection
            typedef policies::NoIndirection<PositionType,ElementType> NoIndirectionType;
            return indirection( HasIndirection< not boost::is_same<IndirectionType, NoIndirectionType>::value >(), 0);
        }


        bool equal(const iter& other) const { return (m_orderedSet == other.m_orderedSet) && (m_pos == other.m_pos); }
        void increment() { advance(1); }
        void decrement() { advance(-1); }
        void advance(PositionType n) { m_pos += n*stride(); }
        const PositionType distance_to(const iter& other) const { return (other.m_pos - m_pos)/ stride(); }

    private:
        inline const PositionType stride() const { return m_orderedSet->StrideType::stride();}

        template<bool> class HasIndirection {};

        template<typename T>
        inline const ElementType& indirection(HasIndirection<true>, T) const { return m_orderedSet->IndirectionType::indirection(m_pos); }

        template<typename T>
        inline const ElementType& indirection(HasIndirection<false>, T) const { return m_pos; }

    private:
        friend class boost::iterator_core_access;

        PositionType m_pos;
        const OrderedSet* m_orderedSet;
    };

public:
    /**
     * \brief Given a position in the Set, return a position in the larger index space
     */
    inline typename IndirectionPolicy::IndirectionResult operator[](PositionType pos) const
    {
            verifyPosition(pos);
            return IndirectionPolicy::indirection( pos * StridePolicyType::stride() + OffsetPolicyType::offset() );
    }

    inline ElementType         at(PositionType pos)         const { return operator[](pos); }


    inline PositionType        size()  const { return SizePolicyType::size(); }
    inline bool                empty() const { return SizePolicyType::empty(); }

  public:   // Functions related to iteration

#ifdef MESHAPI_USE_COUNTING_ITERATOR
    const_iterator          begin() const { return iterator( OffsetPolicyType::offset() ); }
    const_iterator          end()   const { return iterator( size() + OffsetPolicyType::offset() );}
    const_iterator_pair     range() const { return std::make_pair(begin(), end()); }

    iterator          begin() { return iterator( OffsetPolicyType::offset() ); }
    iterator          end()   { return iterator( size() + OffsetPolicyType::offset() );}
    iterator_pair     range() { return std::make_pair(begin(), end()); }

#else
    const_iterator          begin() const { return const_iterator( OffsetPolicyType::offset(), this); }
    const_iterator          end()   const { return const_iterator( size() * StridePolicyType::stride() + OffsetPolicyType::offset(), this);}
    const_iterator_pair     range() const { return std::make_pair(begin(), end()); }
#endif

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
