/**
 * \file OrderedSet.h
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef MESHAPI_ORDERED_SET_H_
#define MESHAPI_ORDERED_SET_H_

#include <cstddef>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include "slic/slic.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

namespace asctoolkit {
namespace meshapi {


/**
 * \class RangeSet
 * \brief Models a set whose elements belong to a contiguous range \f$ \in [lowerIndex,upperIndex) \f$
 * \details The ElementType here needs to be computable as offsets (of PositionType) from the lowerIndex
 *          Examples include: signed and unsigned integral types
 */
  class RangeSet : public Set
  {
  public:
    typedef Set::IndexType                        IndexType;
    typedef Set::PositionType                     PositionType;
    typedef IndexType                             ElementType;


    typedef boost::counting_iterator<ElementType> iterator;
    typedef std::pair<iterator,iterator>          iterator_pair;

    static const NullSet s_nullSet;

  public:
    RangeSet(ElementType size = ElementType(), const Set* parentSet = &s_nullSet )
        : m_lowerIdx(ElementType()), m_upperIdx(size), m_parentSet(parentSet) {}

    RangeSet(ElementType lowerIndex, ElementType upperIndex, const Set* parentSet = &s_nullSet )
        : m_lowerIdx(lowerIndex), m_upperIdx(upperIndex), m_parentSet(parentSet) {}

    PositionType        size()  const { return (m_upperIdx - m_lowerIdx); }

    iterator            begin() const { return iterator(m_lowerIdx); }
    iterator            end()   const { return iterator(m_upperIdx); }
    iterator_pair       range() const { return std::make_pair(begin(), end()); }

    /**
     * \brief Given a position in the Set, return a position in the larger index space
     */
    ElementType operator[](PositionType pos) const { verifyPosition(pos); return pos + m_lowerIdx; }
    ElementType         at(PositionType pos)         const { return operator[](pos); }

    bool                isValid(bool verboseOutput = false) const;

    bool                isSubset() const { return *m_parentSet != s_nullSet; }
    const Set *         parentSet() const { return m_parentSet; }

    bool                isEmpty() const { return size() == PositionType(); }

  private:
    inline void         verifyPosition(PositionType pos)       const
    {
      SLIC_ASSERT_MSG( pos < size()
          , "MeshAPI::RangeSet -- requested out of range element at position "
          << pos << ", but set only has " << size() << " elements." );
    }

  private:
    ElementType m_lowerIdx;
    ElementType m_upperIdx;
    const Set * m_parentSet;

  };


#if 0
/**
 * \brief Two OrderedSets are equal if they have the same cardinality
 * \note Two sets of different types are (currently) considered to be unequal
 */
  inline bool operator==(RangeSet const& firstSet, RangeSet const& otherSet) { return firstSet.size() == otherSet.size(); }

/**
 * \brief Two OrderedSets are equal if they have the same cardinality
 * \note Two sets of different types are (currently) considered to be unequal
 */
  inline bool operator!=(RangeSet const& firstSet, RangeSet const& otherSet) { return !(firstSet==otherSet); }
#endif

} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_ORDERED_SET_H_
