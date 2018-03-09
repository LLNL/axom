/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file OrderedSet.h
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef SLAM_ORDERED_SET_H_
#define SLAM_ORDERED_SET_H_

#include <cstddef>
#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include "Utilities.hpp"
#include "Set.hpp"
#include "NullSet.hpp"

namespace axom
{
namespace slam
{


template<bool IS_SUBRANGE = true>
class RangeSet : public Set
{
public:
  typedef Set::IndexType SetIndex;
  typedef Set::SizeType SizeType;
  typedef Set::PositionType SetPosition;


  typedef boost::counting_iterator<SetIndex> iterator;
  typedef std::pair<iterator,iterator> iterator_pair;

  static const NullSet s_nullSet;

public:
  RangeSet(SizeType size = SizeType(),
           const Set* parentSet = &s_nullSet ) : m_size(size),
    m_parentSet(parentSet) {}

  SizeType size()  const { return m_size; }

  iterator  begin() const { return iterator(0); }
  iterator  end()   const { return iterator(m_size); }
  iterator_pair  range() const { return std::make_pair(begin(), end()); }

  SetIndex operator[](SetPosition pos) const {
    verifyPosition(pos); return pos;
  }
  SetIndex     at(SetPosition pos)         const { return operator[](pos); }

  bool isValid(bool verboseOutput = false) const;

  bool isSubset() const { return *m_parentSet != s_nullSet; }
  const Set* parentSet() const { return m_parentSet; }


private:
  inline void  verifyPosition(SetPosition pos)       const
  {
    SLIC_ASSERT_MSG(
      pos < static_cast<SetPosition>( size() )
      , "SLAM::RangeSet -- requested out of range element at position "
      << pos << ", but set only has " << size() <<
      " elements." );
  }

private:
  SizeType m_size;
  const Set* m_parentSet;

};


#if 0
/**
 * \brief Two OrderedSets are equal if they have the same cardinality
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator==(RangeSet const& firstSet, RangeSet const& otherSet) {
  return firstSet.size() == otherSet.size();
}

/**
 * \brief Two OrderedSets are equal if they have the same cardinality
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator!=(RangeSet const& firstSet, RangeSet const& otherSet) {
  return !(firstSet==otherSet);
}
#endif

} // end namespace slam
} // end namespace axom

#endif //  SLAM_ORDERED_SET_H_
