/**
 * \file IndirectionSet.h
 *
 * \brief Basic API for a set of entities in a simulation
 */

#ifndef MESHAPI_INDIRECTION_SET_H_
#define MESHAPI_INDIRECTION_SET_H_

#include <cstddef>
#include <vector>

#include "slic/slic.hpp"
#include "meshapi/Utilities.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

namespace asctoolkit {
namespace meshapi {


/**
 * \class IndirectionSet
 *
 * \brief An indexed set (a tuple) of entities in a simulation
 *
 * A container class for a set of entities in a simulation. Each entity has an index.
 *
 * Below is an initial implementation for a set with explicit indexes (encoded here using a vector).
 */
class IndirectionSet : public Set
{
public:

  typedef Set::IndexType                            IndexType;
  typedef Set::PositionType                         PositionType;
  typedef IndexType                                 ElementType;

  typedef std::vector<ElementType>                  ArrType;
  typedef ArrType::iterator                         iterator;
  typedef std::pair<iterator, iterator>             iterator_pair;

  typedef ArrType::const_iterator                   const_iterator;
  typedef std::pair<const_iterator, const_iterator> const_iterator_pair;


  static const NullSet s_nullSet;

public:
  IndirectionSet (const Set* parentSet = &s_nullSet) : m_parentSet(parentSet) {}
  ~IndirectionSet () {}

  /**
   * \brief Random access to the entities of the set.  Checked in debug builds.
   * \param pos The index of the desired element
   * \return A reference to the encoded index
   * \pre pos must be less than the number of elements in the set ( size() )
   */
  ElementType&       operator [](PositionType pos)        { verifyPosition(pos); return m_entities[pos]; }

  /**
   * \brief Random access to the entities of the set (const version).  Checked in debug builds.
   * \param pos The index of the desired element
   * \return A const reference to the encoded index
   * \pre pos must be less than the number of elements in the set ( size() )
   */
  ElementType const& operator [](PositionType pos) const { verifyPosition(pos); return m_entities[pos]; }

  /**
   * \brief Random access to the entities of the set -- implementation of Set's virtual function at()
   * \param pos The index of the desired element
   * \return The element at position pos (via the indirection map)
   * \pre pos must be less than the number of elements in the set ( size() )
   */
  ElementType                 at(PositionType pos) const { return operator[](pos); }

  /**
   * \brief Get the number of entities in the set
   * \return The number of entities in the set.
   */
  PositionType                size() const { return m_entities.size(); }

  /**
   * \return An iterator to the beginning of the entities
   */
  iterator                    begin()            { return m_entities.begin(); }

  /**
   * \return A const iterator to the beginning of the entities
   */
  const_iterator              begin() const { return m_entities.begin(); }

  /**
   * \return An iterator to the end of the entities.
   */
  iterator                    end()              { return m_entities.end(); }

  /**
   * \return A const iterator to the end of the entities
   */
  const_iterator              end() const { return m_entities.end(); }

  /**
   * \return A pair of begin/end iterators
   */
  iterator_pair               range()        { return std::make_pair(begin(), end()); }

  /**
   * \return A pair of begin/end iterators
   */
  const_iterator_pair         range() const { return std::make_pair(begin(), end()); }


  /**
   * \brief Checks that the set is valid
   */
  bool                        isValid(bool verboseOutput = false) const;

  /**
   * \brief Determines if the Set is a Subset of another set.
   * \return true if the set is a subset of another set, otherwise false.
   */
  bool                        isSubset() const { return *m_parentSet == s_nullSet; }

  /**
   * \brief Checks whether the set has any elements
   */
  bool                        isEmpty() const { return m_entities.empty(); }


  /**
   * \return A pointer to the parent set.  Parent is a NullSet if isSubset() is false
   */
  const Set*                  parentSet()    const { return m_parentSet; }

public:
  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying set data
   * \note We will have to figure out a good way to limit this access to situations where it makes sense.
   */

  /// \{

  ArrType &       data()        { return m_entities; }
  const ArrType & data() const { return m_entities; }

  /// \}

private:
  void verifyPosition(PositionType pos) const
  {
    SLIC_ASSERT_MSG( pos < size()
      , "MeshAPI::IndirectionSet -- requested out of range element at position "
      << pos << ", but set only has " << size() << " elements." );
  }
private:

  ArrType m_entities;
  const Set*      m_parentSet;
};

#if 0
/**
 * \brief Two IndirectionSets are equal if they have the same cardinality and all indexes are the same (and in the same order)
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator==(IndirectionSet const& firstSet, IndirectionSet const& otherSet)
{
  if(firstSet.size() != otherSet.size() )
    return false;
  for(IndexType idx = 0; idx< firstSet.size(); ++idx)
  {
    if( firstSet[idx] != otherSet[idx] )
      return false;
  }
  return true;
}

/**
 * \brief Two IndirectionSets are equal if they have the same cardinality and all indexes are the same (and in the same order)
 * \note Two sets of different types are (currently) considered to be unequal
 */
inline bool operator!=(IndirectionSet const& firstSet, IndirectionSet const& otherSet) { return !(firstSet==otherSet); }
#endif



} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_INDIRECTION_SET_H_
