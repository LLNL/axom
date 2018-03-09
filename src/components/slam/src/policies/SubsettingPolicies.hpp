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
 * \file SubsettingPolicies.hpp
 *
 * \brief Subsetting policies for SLAM
 *
 * Subsetting policies encompass the type and availability of a set's parent
 * A valid subset policy must support the following interface:
 *   * [required]
 *   * isSubset(): bool -- returns whether the set is a subset of another set
 *   * parentSet() : ParentSetType -- returns a pointer to the parent set.
 *                                     AXOM_NULLPTR when isSubset() is false
 *   * isValid() : bool -- indicates whether the Subsetting policy of the set
 *     is valid
 *   * [optional]
 *   * operator(): IntType -- alternate accessor for indirection
 */

#ifndef SLAM_POLICIES_SUBSET_H_
#define SLAM_POLICIES_SUBSET_H_

#include "axom/Macros.hpp"

#include <set>

namespace axom
{
namespace slam
{
namespace policies
{

/**
 * \name OrderedSet_Subsetting_Policies
 * \brief A few default policies for the subsetting of an OrderedSet
 */

/// \{

struct NoSubset
{
  static const NullSet s_nullSet;
  typedef const Set ParentSetType;

  NoSubset() {}

  // This empty .ctor is here to satisfy the SubsettingPolicy API
  NoSubset(ParentSetType*) {}

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool                  isSubset() const { return false; }
  const ParentSetType* parentSet() const { return AXOM_NULLPTR; }

  template<typename OrderedSetIt>
  bool                  isValid(OrderedSetIt, OrderedSetIt, bool) const
  {
    return true;
  }
};

struct VirtualParentSubset
{
  static NullSet s_nullSet;

  typedef Set ParentSetType;

  VirtualParentSubset(ParentSetType* parSet = &s_nullSet)
    : m_parentSet(parSet) {}

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool        isSubset() const { return *m_parentSet != s_nullSet; }
  const Set* parentSet() const { return m_parentSet; }
  Set*&       parentSet() { return m_parentSet; }

  template<typename OrderedSetIt>
  bool        isValid(OrderedSetIt beg, OrderedSetIt end,
                      bool verboseOutput = false) const
  {
    // We allow parent sets to be null (i.e. the subset feature is deactivated)
    if( !isSubset() || m_parentSet == AXOM_NULLPTR)
      return true;

    // Next, check if child is empty -- null set is a subset of all sets
    bool childIsEmpty = (beg == end);
    if( childIsEmpty )
      return true;

    // Next, since child has at least one element, the parent cannot be empty
    if(verboseOutput)
    {
      bool bValid = ( m_parentSet->size() > 0);
      AXOM_DEBUG_VAR(bValid);
      SLIC_CHECK_MSG(
        bValid,
        "VirtualParentSubset -- if we are a subset and input set is "
        <<"non-empty then parent set must be non-empty");
    }

    // At this point, parent and child are both non-null
    //  -- ensure that all elts of child are in parent
    std::set<typename Set::ElementType> pSet;
    for(typename Set::PositionType pos = 0 ; pos < m_parentSet->size() ; ++pos)
      pSet.insert( m_parentSet->at(pos));
    for( ; beg != end ; ++beg)
    {
      if( pSet.find(*beg) == pSet.end())
      {
        // For now, we only warn about the first failure '
        // Note: in the future, we might want to show all problems.
        SLIC_CHECK_MSG(
          verboseOutput,
          "VirtualParentSubset :: parent set does not contain element "
          << *beg << " so child cannot be a subset of parent");
        return false;
      }
    }
    return true;
  }

private:
  Set* m_parentSet;
};

template<typename TheParentSetType>
struct ConcreteParentSubset
{
  typedef const TheParentSetType ParentSetType;

  ConcreteParentSubset(ParentSetType* parSet = AXOM_NULLPTR)
    : m_parentSet(parSet) {}

  /**
   * \brief Checks whether the set containing this policy class is a subset
   */
  bool                  isSubset()  const
  {
    return m_parentSet != AXOM_NULLPTR;
  }

  ParentSetType* const& parentSet() const { return m_parentSet; }
  ParentSetType* &      parentSet()       { return m_parentSet; }


  template<typename OrderedSetIt>
  bool isValid( OrderedSetIt beg, OrderedSetIt end,
                bool verboseOutput = false) const
  {
    // We allow parent sets to be null (i.e. the subset feature is deactivated)
    if( !isSubset() )
      return true;

    // Next, check if child is empty -- null set is a subset of all sets
    bool childIsEmpty = (beg == end);
    if( childIsEmpty )
      return true;

    // Next, since child has at least one element, the parent cannot be empty
    if(verboseOutput)
    {
      bool bValid = (m_parentSet->size() > 0);
      AXOM_DEBUG_VAR(bValid);
      SLIC_CHECK_MSG(
        bValid,
        "VirtualParentSubset -- if input set is non-empty "
        << " parent set must be non-empty");
    }

    // At this point, parent and child are both non-null
    std::set<typename Set::ElementType> pSet;
    for(typename ParentSetType::PositionType pos = 0 ;
        pos < m_parentSet->size() ; ++pos)
      pSet.insert( (*m_parentSet)[pos] );
    for( ; beg != end ; ++beg)
    {
      if( pSet.find(*beg) == pSet.end())
      {
        // For now, we only warn about the first failure '
        // Note: in the future, we might want to show all problems.
        SLIC_CHECK_MSG(
          verboseOutput,
          "ConcreteParentSubset :: parent set does not contain element "
          << *beg << " so child cannot be a subset of parent");
        return false;
      }
    }
    return true;
  }

private:
  ParentSetType* m_parentSet;
};

/// \}

} // end namespace policies
} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICIES_SUBSET_H_
