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
 * \file DynamicVariableRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities from
 * the first set can be related to an arbitrary number of entities from the
 * second set. This relation is dynamic; the related entities can change
 * at runtime.
 */

#ifndef SLAM_DYNAMIC_VARIABLE_RELATION_HPP_
#define SLAM_DYNAMIC_VARIABLE_RELATION_HPP_

#include "axom/config.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/Relation.hpp"

#include <vector>
//#include <iostream>

namespace axom
{
namespace slam
{

class DynamicVariableRelation : public Relation
{
public:
  typedef Relation::SetPosition SetPosition;

  typedef std::vector<SetPosition>
    RelationVec;
  typedef RelationVec::iterator RelationVecIterator;
  typedef std::pair<RelationVecIterator,RelationVecIterator>
    RelationVecIteratorPair;

  typedef RelationVec::const_iterator RelationVecConstIterator;
  typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>
    RelationVecConstIteratorPair;

  typedef std::vector< RelationVec>
    RelationsContainer;
  typedef RelationsContainer::const_iterator RelationsContainerCIt;
  typedef RelationsContainer::iterator RelationsContainerIt;

public:
  DynamicVariableRelation (Set* fromSet = &s_nullSet,
                           Set* toSet = &s_nullSet);
  ~DynamicVariableRelation(){}

public:
#ifdef AXOM_USE_CXX11
  /// \name DynamicVariableRelation iterator interface
  /// @{
  RelationVecConstIterator begin(SetPosition fromSetIndex)       const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).begin();
  }

  RelationVecConstIterator end(SetPosition fromSetIndex)         const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).end();
  }

  RelationVecConstIteratorPair range(SetPosition fromSetIndex)   const
  {
    return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
  }
  /// @}

#endif // AXOM_USE_CXX11


  RelationVec const& operator[](SetPosition fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

  SetPosition size(SetPosition fromSetIndex)                  const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).size();
  }

  bool isValid(bool verboseOutput = false) const;


public:     // Modifying functions

  void insert(SetPosition fromSetIndex, SetPosition toSetIndex)
  {
    verifyPosition(fromSetIndex);
    m_relationsVec[fromSetIndex].push_back(toSetIndex);

  }

  RelationVec& operator[](SetPosition fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

public:
  /**
   * \name DirectDataAccess
   * \brief Accessor functions to get the underlying relation data for each
   *  element

   * \note We will have to figure out a good way
   * to limit this access to situations where it makes sense.
   */

  /// \{

  /**
   * \brief Access the set of positions in the 'toSet'
   * associated with the given position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet'
   * whose relation data (in the 'toSet') we are requesting
   */
  RelationVec &       data(SetPosition fromSetPos)
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /**
   * \brief Access the set of positions in the 'toSet'
   * associated with the given position in 'fromSet'
   *
   * \param fromSetPos The position within the 'fromSet'
   * whose relation data (in the 'toSet') we are requesting
   */
  const RelationVec & data(SetPosition fromSetPos) const
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /// \}

private:
  inline void
  verifyPosition(SetPosition AXOM_DEBUG_PARAM(fromSetIndex))        const
  {
    SLIC_ASSERT_MSG(
      fromSetIndex >= 0 &&
      fromSetIndex < static_cast<SetPosition>(m_fromSet->size() ),
      "Index " << fromSetIndex
               << " out of range [0," << m_fromSet->size() << ")");
  }

  inline RelationVec &      fromSetRelationsVec(SetPosition fromSetIndex)
  {
    return m_relationsVec[fromSetIndex];
  }
  inline RelationVec const& fromSetRelationsVec(SetPosition fromSetIndex) const
  {
    return m_relationsVec[fromSetIndex];
  }



private:

  Set* m_fromSet;
  Set* m_toSet;

  RelationsContainer m_relationsVec;
};


} // end namespace slam
} // end namespace axom

#endif // SLAM_DYNAMIC_VARIABLE_RELATION_HPP_
