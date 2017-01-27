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

#ifndef SLAM_STATIC_VARIABLE_RELATION_HPP_
#define SLAM_STATIC_VARIABLE_RELATION_HPP_


#include <vector>

//#include <iostream>

#include "axom/config.hpp"   // for AXOM_USE_BOOST

#include "slic/slic.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/NullSet.hpp"
#include "slam/Relation.hpp"


namespace axom {
namespace slam    {

  class StaticVariableRelation : public Relation
  {
  public:
    typedef Relation::SetPosition                                         SetPosition;

    typedef std::vector<SetPosition>                                      RelationVec;

#ifdef AXOM_USE_BOOST
    typedef RelationVec::iterator                                         RelationVecIterator;
    typedef std::pair<RelationVecIterator,RelationVecIterator>            RelationVecIteratorPair;

    typedef RelationVec::const_iterator                                   RelationVecConstIterator;
    typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>  RelationVecConstIteratorPair;
#endif // AXOM_USE_BOOST

    typedef OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>,      // TODO: change this to a compile time size if/when parent is compile time
        policies::RuntimeOffsetHolder<Set::PositionType>,
        policies::StrideOne<Set::PositionType>,
        policies::STLVectorIndirection<Set::PositionType, Set::ElementType> > RelationSet;

  public:
    StaticVariableRelation (Set* fromSet = &s_nullSet, Set* toSet = &s_nullSet);
    ~StaticVariableRelation(){}
    /**
     * \note TODO: swap this out for data in the datastore
     */
    void                      bindRelationData(RelationVec const& beginsVec, RelationVec const& toOffsets);

 #ifdef AXOM_USE_BOOST
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
  #endif // AXOM_USE_BOOST

    /**
     * This function returns the OrderedSet of all elements in the toSet related to 'fromSetElt' in the fromSet.
     */
    const RelationSet operator[](SetPosition fromSetElt) const
    {
      typedef RelationSet::SetBuilder SetBuilder;
      return SetBuilder()
             .size( elemSize(fromSetElt))
             .offset( toSetBeginIndex(fromSetElt) )
             .data( &m_toSetIndicesVec)
      ;
    }

    SetPosition size(SetPosition fromSetIndex)                  const
    {
      return elemSize(fromSetIndex);
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
    inline SetPosition elemSize(SetPosition fromSetIndex) const
    {
      verifyPosition(fromSetIndex);
      return toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex);
    }

    inline void         verifyPosition(SetPosition AXOM_DEBUG_PARAM(fromSetIndex))    const
    {
      SLIC_ASSERT( fromSetIndex >= 0 && fromSetIndex <  m_fromSet->size()  );
    }
    inline SetPosition  toSetBeginIndex(SetPosition fromSetIndex)   const { return m_fromSetBeginsVec[fromSetIndex]; }
    inline SetPosition  toSetEndIndex(SetPosition fromSetIndex)     const { return m_fromSetBeginsVec[fromSetIndex + 1]; }



  private:

    Set* m_fromSet;
    Set* m_toSet;

    RelationVec m_fromSetBeginsVec;           // vector of size m_fromSet.size() + 1 that points into the to_set vectors
    RelationVec m_toSetIndicesVec;            // vector of toSet entries
  };


} // end namespace slam
} // end namespace axom

#endif // SLAM_STATIC_VARIABLE_RELATION_HPP_
