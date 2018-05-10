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
 * \file DynamicConstantRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities from
 * the first set can be related to a constant number of entities from the
 * second set. This relation is dynamic; the related entities can change at
 * runtime.
 */

#ifndef SLAM_DYNAMIC_CONSTANT_RELATION_HPP_
#define SLAM_DYNAMIC_CONSTANT_RELATION_HPP_

#include "axom/config.hpp"

#include "slic/slic.hpp"
#include "slam/Set.hpp"
#include "slam/Relation.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/DynamicSet.hpp"
#include "slam/CardinalityPolicies.hpp"

#include <vector>

namespace axom
{
namespace slam
{

/**
 * \class DynamicConstantRelation
 * \brief  A relation class with constant cardinality that supports
 * adding/removing set relations.
 *
 * \detail An entry is considered valid if its set's entry is valid
 * and at least one of its relation is valid (ie. not equal to INVALID_INDEX).
 *
 * \note The current implementation fixes the value of INVALID_INDEX to -1.
 * A future update will allow users to set the value of INVALID_INDEX to a
 * more convenient value, when necessary.
 */
template<typename CardinalityPolicy>
class DynamicConstantRelation : public /*Relation,*/ CardinalityPolicy
{
public:
  enum
  {
    INVALID_INDEX = -1
  };

  typedef Relation::SetPosition SetPosition;
  typedef std::vector<SetPosition>                  RelationVec;

  typedef DynamicSet<>                              FromSetType;
  typedef DynamicSet<>                              ToSetType;

  typedef typename CardinalityPolicy::
    RelationalOperatorSizeType BeginsSizePolicy;

  typedef axom::slam::policies::
    STLVectorIndirection<SetPosition,SetPosition>   STLIndirection;
  typedef OrderedSet<
      BeginsSizePolicy,
      policies::RuntimeOffset<SetPosition>,
      policies::StrideOne<SetPosition>,
      STLIndirection >                              RelationSet;

#ifdef AXOM_USE_CXX11
  typedef typename RelationSet::iterator RelationIterator;
  typedef typename RelationSet::iterator_pair RelationIteratorPair;

  typedef typename RelationSet::const_iterator RelationConstIterator;
  typedef typename RelationSet::const_iterator_pair RelationConstIteratorPair;
#endif // AXOM_USE_CXX11

public:

  /**
   * \brief Default constructor with empty set for toSet and fromSet
   */
  DynamicConstantRelation () :
    m_fromSet( EmptySetTraits<FromSetType>::emptySet()  ),
    m_toSet(  EmptySetTraits<ToSetType>::emptySet() )
  {
    m_relationCardinality = CardinalityPolicy::size( 0 );
  }

  /**
   * \brief construct a DynamicConstantRelation from the given fromSet to toSet
   */
  DynamicConstantRelation (FromSetType* fromSet, ToSetType* toSet)
    : CardinalityPolicy(
      EmptySetTraits<Set>::isEmpty(fromSet) ? 0 : fromSet->size() ),
    m_fromSet(fromSet ),
    m_toSet( toSet )
  {
    m_relationCardinality = CardinalityPolicy::size( 0 );
    m_relationsVec.resize(m_relationCardinality*fromSet->size(), INVALID_INDEX);
  };

  ~DynamicConstantRelation(){};


public:
#ifdef AXOM_USE_CXX11
  /// \name DynamicConstantRelation iterator interface
  /// @{

  RelationIterator          begin(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].begin();
  }

  RelationConstIterator     begin(SetPosition fromSetInd ) const
  {
    return (*this)[fromSetInd].begin();
  }

  RelationIterator          end(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].end();
  }

  RelationConstIterator     end(SetPosition fromSetInd)    const
  {
    return (*this)[fromSetInd].end();
  }

  RelationIteratorPair      range(SetPosition fromSetInd)
  {
    return (*this)[fromSetInd].range();
  }

  RelationConstIteratorPair range(SetPosition fromSetInd)  const
  {
    return (*this)[fromSetInd].range();
  }

  /// @}

#endif // AXOM_USE_CXX11

public:

  RelationSet const at(SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return operator[](fromSetIndex);
  }

  RelationSet const operator[](SetPosition fromSetIndex) const
  {
    typedef typename RelationSet::SetBuilder SetBuilder;
    return SetBuilder()
           //.size( CardinalityPolicy::size(fromSetIndex) )
           .size( m_relationCardinality )
           //.offset( CardinalityPolicy::offset( fromSetIndex) )
           .offset( fromSetIndex * m_relationCardinality )
           .data( &m_relationsVec);
  }

  SetPosition size(SetPosition fromSetIndex ) const
  {
    verifyPosition(fromSetIndex);
    return m_relationCardinality;
  }

  SetPosition size() const
  {
    return m_relationsVec.size() / m_relationCardinality;
  }

  SetPosition numberOfValidEntries() const
  {
    SetPosition nvalid = 0;
    const int N = size();
    for( int i=0 ; i< N ; ++i)
    {
      nvalid += isValidEntry(i);
    }
    return nvalid;
  }

  /**
   * \brief return if an entry is valid or not.
   * \detailed an entry is considered valid if it has at least one valid value
   */
  bool isValidEntry(SetPosition idx) const
  {
    if( idx >= 0 && idx < (int)m_relationsVec.size()/m_relationCardinality )
    {
      for (int i = 0 ; i < m_relationCardinality ; ++i)
      {
        if (m_relationsVec[idx * m_relationCardinality + i] != INVALID_INDEX)
          return true;
      }
    }
    return false;
  }

  bool isValid(bool verboseOutput = false) const;

public:   // Modifying functions

  void insert(SetPosition fromSetIndex, SetPosition toSetIndex)
  {
    expandSizeIfNeeded(fromSetIndex+1);

    verifyPosition( fromSetIndex );

    //find the first invalid place to put it
    for(int i=0 ; i<m_relationCardinality ; ++i)
    {
      SetPosition idx = m_relationCardinality*fromSetIndex + i;
      if( m_relationsVec[idx] == INVALID_INDEX )
      {
        m_relationsVec[idx] = toSetIndex;
        return;
      }
    }

    //The entry was not inserted
    SLIC_WARNING(
      "Relation from "
      << fromSetIndex <<" to " << toSetIndex
      << " was not inserted because the entry is full.");
  }

  /**
   * \brief Function to modify the value at offset \a offset of the
   * fromSet index \fromSetIndex to the value \a toSetIndex
   *
   * \note This is a temporary stopgap function until operator[]
   *  allows us to modify values.
   *
   * \detail This should be replaced with operator[] which returns a non-const
   * RelationSet so users can do relation[fromSetIndex][offset] = toSetIndex;
   */
  void modify(SetPosition fromSetIndex,
              SetPosition offset,
              SetPosition toSetIndex)
  {
    expandSizeIfNeeded( fromSetIndex + 1);
    m_relationsVec[m_relationCardinality * fromSetIndex + offset] = toSetIndex;
  }

  /**
   * \brief Mark all values in entry \a fromSetIndex as invalid.
   */
  void remove(SetPosition fromSetIndex)
  {
    if(!isValidEntry(fromSetIndex))
      return;

    for(int i=0 ; i<m_relationCardinality ; ++i)
    {
      m_relationsVec[m_relationCardinality*fromSetIndex + i] = INVALID_INDEX;
    }
  }

public:
  /**
   * \brief Access the set of positions in the 'toSet'
   * associated with the given position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet' whose relation data
   * (in the 'toSet') we are requesting
   */
  RelationVec &       data()       { return m_relationsVec; }

  /**
   * \brief Access the set of positions in the 'toSet' associated with the given
   * position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet' whose relation data
   * (in the 'toSet') we are requesting
   */
  const RelationVec & data() const { return m_relationsVec; }

private:

  void expandSizeIfNeeded(SetPosition s)
  {
    if(s > (int)m_relationsVec.size()/m_relationCardinality)
    {
      m_relationsVec.resize(s*m_relationCardinality, INVALID_INDEX);
    }
  }

  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(fromSetIndex))
  const
  {
    SLIC_ASSERT_MSG(
      fromSetIndex >= 0 &&
      fromSetIndex < static_cast<SetPosition>(m_fromSet->size() ),
      "Index " << fromSetIndex
               << " out of range [0," << m_fromSet->size() <<  ")");
  }


private:

  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  int m_relationCardinality;
  RelationVec m_relationsVec;
};


/**
 * \brief Checks whether the relation is valid.
 */
template<typename CardinalityPolicy>
bool DynamicConstantRelation<CardinalityPolicy>::isValid(bool verboseOutput)
const
{
  std::stringstream errSstr;

  bool setsAreValid = true;
  bool relationdataIsValid = true;

  // Check if the sets are valid
  bool isFromSetNull = (m_fromSet == AXOM_NULLPTR);
  bool isToSetNull   = (m_toSet == AXOM_NULLPTR);

  if(isFromSetNull || isToSetNull)
  {
    if(verboseOutput)
    {
      errSstr
        << "\n\t Static relations require both the fromSet "
        << "and toSet to be non-null"
        << "\n\t -- fromSet was " << (isFromSetNull ? "" : " not ") << "null"
        << "\n\t -- toSet was " << (isToSetNull ? "" : " not ") << "null";
    }

    setsAreValid = false;
  }


  // Check the sizes of fromSet matches relationVec
  if( setsAreValid )
  {
    if( m_fromSet->size() * m_relationCardinality !=
        (int)m_relationsVec.size() )
    {
      if(verboseOutput)
      {
        errSstr << "\n\t Size of relationVec does not match toSet size. "
                << "\n\t -- fromSet size is " << m_fromSet->size() << ","
                << "\n\t -- m_relationsVec size is " << m_relationsVec.size()
                << ".";
      }
      setsAreValid =  false;
    }
  }

  // Check if the relation data is valid
  if(setsAreValid)
  {
    // Check that invalid set entry points to invalid relation.
    // Note: the reverse can be valid. ie. valid set entry may have invalid
    // relation entry.
    for(SetPosition pos = 0 ; pos < (int) m_fromSet->size() ; ++pos)
    {
      if( m_fromSet->at(pos) == FromSetType::INVALID_ENTRY &&
          isValidEntry(pos) )
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* invalid entries in fromSet has a valid Relation."
                  << "\n\t-- at index: " << pos;
        }
        relationdataIsValid = false;
      }
    }

    // Check that all relation indices are in range for m_toSet
    for(SetPosition pos = 0 ; pos < (int) m_relationsVec.size() ; ++pos)
    {
      if ( m_relationsVec[pos] != INVALID_INDEX &&
           !m_toSet->isValidEntry( m_relationsVec[pos]) )
      {
        if(verboseOutput)
        {
          errSstr
            << "\n\t* Relation index out of range or invalid."
            << "\n\t-- position " << pos / m_relationCardinality
            << "-" <<  pos % m_relationCardinality
            << " with value: " << m_relationsVec[pos]
            << " needs to be in range [0," << m_toSet->size()
            << ") and index a valid entry.";
        }
        relationdataIsValid = false;
      }
    }
  }

  // We are done.  Output the messages if applicable and return
  bool bValid = setsAreValid && relationdataIsValid;

  if(verboseOutput && !bValid)
  {
    SLIC_DEBUG( errSstr.str() );
  }

  return bValid;
}

} // end namespace slam
} // end namespace axom

#endif // SLAM_DYNAMIC_CONSTANT_RELATION_HPP_
