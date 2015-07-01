/**
 * \file StaticVariableRelation.h
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to an arbitrary number of entities from the second set
 *
 */

#ifndef MESHAPI_STATIC_VARIABLE_RELATION_HPP_
#define MESHAPI_STATIC_VARIABLE_RELATION_HPP_

#include <vector>

//#include <iostream>

#include "slic/slic.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

class StaticVariableRelation : public Relation
{
private:
  /**
   * A small helper class to allow double subscripting on the relation
   */
  class SubscriptProxy {
  public:
    SubscriptProxy(RelationVecConstIterator it, SetPosition size) : m_iter(it), m_size(size) {}
    SetPosition const& operator[](SetPosition index) const
    {
      SLIC_ASSERT_MSG( index < m_size, "Inner array access out of bounds."
        << "\n\tPresented value: " << index
        << "\n\tMax allowed value: " << static_cast<int>(m_size - 1));
      return m_iter[index];
    }
  private:
    RelationVecConstIterator m_iter;
    SetPosition m_size;
  };
public:
  typedef Relation::SetPosition                                         SetPosition;

  typedef std::vector<SetPosition>                                      RelationVec;
  typedef RelationVec::iterator                                         RelationVecIterator;
  typedef std::pair<RelationVecIterator,RelationVecIterator>            RelationVecIteratorPair;

  typedef RelationVec::const_iterator                                   RelationVecConstIterator;
  typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>  RelationVecConstIteratorPair;



public:
  StaticVariableRelation (Set* fromSet = &s_nullSet, Set* toSet = &s_nullSet);
  virtual ~StaticVariableRelation(){}
  /**
   * \note TODO: swap this out for data in the datastore
   */
  void                      bindRelationData(RelationVec const& beginsVec, RelationVec const& toOffsets);

  RelationVecConstIterator  begin(SetPosition fromSetIndex)       const
  {
    verifyPosition(fromSetIndex);
    return m_toSetIndicesVec.begin() + toSetBeginIndex(fromSetIndex);
  }

  RelationVecConstIterator end(SetPosition fromSetIndex)         const
  {
    verifyPosition(fromSetIndex);
    return m_toSetIndicesVec.begin() + toSetEndIndex(fromSetIndex);
  }

  RelationVecConstIteratorPair range(SetPosition fromSetIndex)   const
  {
    return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
  }

  SubscriptProxy const operator[](SetPosition fromSetElt) const
  {
    return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
  }

  SetPosition size(SetPosition fromSetIndex)                  const
  {
    verifyPosition(fromSetIndex);
    return toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex);
  }

  bool isValid(bool verboseOutput = false) const;

public:

  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying relation data
   * \note We will have to figure out a good way to limit this access to situations where it makes sense.
   */

  /// \{

  /**
   * \brief Helper function to access the underlying relation data
   * \note The relation currently 'owns' the underlying vector.
   *       This will be changing soon, and we will only have a reference/pointer to the data.
   */
  RelationVec &       fromSetBeginsData()       { return m_fromSetBeginsVec; }

  /**
   * \brief Helper function to access the underlying relation data
   * \note The relation currently 'owns' the underlying vector.
   *       This will be changing soon, and we will only have a reference/pointer to the data.
   */
  const RelationVec & fromSetBeginsData() const { return m_fromSetBeginsVec; }

  /**
   * \brief Helper function to access the underlying relation data
   * \note The relation currently 'owns' the underlying vector.
   *       This will be changing soon, and we will only have a reference/pointer to the data.
   */
  RelationVec &       toSetPositionsData()       { return m_toSetIndicesVec; }

  /**
   * \brief Helper function to access the underlying relation data
   * \note The relation currently 'owns' the underlying vector.
   *       This will be changing soon, and we will only have a reference/pointer to the data.
   */
  const RelationVec & toSetPositionsData() const { return m_toSetIndicesVec; }

  /// \}

private:
  inline void         verifyPosition(SetPosition fromSetIndex)       const { SLIC_ASSERT( fromSetIndex <  m_fromSet->size()  ); }
  inline SetPosition  toSetBeginIndex(SetPosition fromSetIndex)   const { return m_fromSetBeginsVec[fromSetIndex]; }
  inline SetPosition  toSetEndIndex(SetPosition fromSetIndex)     const { return m_fromSetBeginsVec[fromSetIndex + 1]; }



private:

  Set* m_fromSet;
  Set* m_toSet;

  RelationVec m_fromSetBeginsVec;             // vector of size m_fromSet.size() + 1 that points into the to_set vectors
  RelationVec m_toSetIndicesVec;              // vector of toSet entries
};


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_VARIABLE_RELATION_HPP_
