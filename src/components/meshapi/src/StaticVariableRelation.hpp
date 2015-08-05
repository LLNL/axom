/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file StaticVariableRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to an arbitrary number of entities from the second set
 *
 */

#ifndef MESHAPI_STATIC_VARIABLE_RELATION_HPP_
#define MESHAPI_STATIC_VARIABLE_RELATION_HPP_

#ifndef MESHAPI_STATIC_VARIABLE_RELATION_ITERATOR_USE_PROXY
//  #define MESHAPI_STATIC_VARIABLE_RELATION_ITERATOR_USE_PROXY
#endif


#include <vector>

//#include <iostream>

#include "slic/slic.hpp"
#include "meshapi/OrderedSet.hpp"
#include "meshapi/NullSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

  class StaticVariableRelation : public Relation
  {

#ifdef MESHAPI_STATIC_VARIABLE_RELATION_ITERATOR_USE_PROXY
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
#endif

  public:
    typedef Relation::SetPosition                                         SetPosition;

    typedef std::vector<SetPosition>                                      RelationVec;
    typedef RelationVec::iterator                                         RelationVecIterator;
    typedef std::pair<RelationVecIterator,RelationVecIterator>            RelationVecIteratorPair;

    typedef RelationVec::const_iterator                                   RelationVecConstIterator;
    typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>  RelationVecConstIteratorPair;

    typedef OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>      // TODO: change this to a compile time size if/when parent is compile time
                      , policies::RuntimeOffsetHolder<Set::PositionType>
                      , policies::StrideOne<Set::PositionType>
                      , policies::STLVectorIndirection<Set::PositionType, Set::ElementType> > RelationSet;


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

#ifdef MESHAPI_STATIC_VARIABLE_RELATION_ITERATOR_USE_PROXY
    const SubscriptProxy operator[](SetPosition fromSetElt) const
    {
      return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
    }
#else
    /**
     * This function returns the OrderedSet of all elements in the toSet related to 'fromSetElt' in the fromSet.
     */
    const RelationSet operator[](SetPosition fromSetElt) const
    {
        // Note -- we need a better way to initialize an indirection set
        RelationSet rel(size(fromSetElt),toSetBeginIndex(fromSetElt) );
        rel.data() = &m_toSetIndicesVec;

        return rel;
    }
#endif

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
    inline void         verifyPosition(SetPosition fromSetIndex)    const { SLIC_ASSERT( fromSetIndex <  m_fromSet->size()  ); }
    inline SetPosition  toSetBeginIndex(SetPosition fromSetIndex)   const { return m_fromSetBeginsVec[fromSetIndex]; }
    inline SetPosition  toSetEndIndex(SetPosition fromSetIndex)     const { return m_fromSetBeginsVec[fromSetIndex + 1]; }



  private:

    Set* m_fromSet;
    Set* m_toSet;

    RelationVec m_fromSetBeginsVec;           // vector of size m_fromSet.size() + 1 that points into the to_set vectors
    RelationVec m_toSetIndicesVec;            // vector of toSet entries
  };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_VARIABLE_RELATION_HPP_
