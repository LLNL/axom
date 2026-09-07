// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file DynamicVariableRelation.hpp
 *
 * \brief Editable connectivity with a variable number of entries per from-set position.
 */

#include "axom/config.hpp"
#include "axom/slic.hpp"

#include "axom/slam/policies/PolicyTraits.hpp"

#include "axom/slam/Set.hpp"
#include "axom/slam/Utilities.hpp"

#include <vector>
#include <sstream>
#include <iterator>

namespace axom::slam
{
/**
 * \class DynamicVariableRelation
 * \brief Store and edit a collection of to-set positions for each from-set position.
 *
 * For a vertex-to-cell relation, each vertex can have a different number of
 * incident cells. insert(from, to) appends a to-set position. operator[] and
 * data() provide direct access to the corresponding std::vector.
 *
 * The relation owns its vectors and borrows both sets. The from-set size must
 * remain equal to its size at construction. Referenced sets must outlive the
 * relation. Changes to a vector can invalidate its iterators and references.
 */
template <typename FirstSetType = slam::Set<>, typename SecondSetType = slam::Set<>>
class DynamicVariableRelation
{
public:
  using FromSetType = FirstSetType;
  using ToSetType = SecondSetType;

  using FromPositionType = typename FromSetType::PositionType;
  using ToPositionType = typename ToSetType::PositionType;
  using FlatPositionType = detail::default_flat_position_t<FromPositionType, ToPositionType>;

  using RelationVec = std::vector<ToPositionType>;
  using RelationVecIterator = typename RelationVec::iterator;
  using RelationVecIteratorPair = std::pair<RelationVecIterator, RelationVecIterator>;
  using RelationVecConstIterator = typename RelationVec::const_iterator;
  using RelationVecConstIteratorPair = std::pair<RelationVecConstIterator, RelationVecConstIterator>;
  using RelationsContainer = std::vector<RelationVec>;
  using RelationsContainerCIt = typename RelationsContainer::const_iterator;
  using RelationsContainerIt = typename RelationsContainer::iterator;

public:
  DynamicVariableRelation(FirstSetType* fromSet = policies::EmptySetTraits<FirstSetType>::emptySet(),
                          SecondSetType* toSet = policies::EmptySetTraits<SecondSetType>::emptySet())
    : m_fromSet(fromSet)
    , m_toSet(toSet)
  {
    if(m_fromSet)
    {
      m_relationsVec.resize(m_fromSet->size());
    }
  }

  ~DynamicVariableRelation() { }

public:
  /// \name DynamicVariableRelation iterator interface
  /// @{
  RelationVecConstIterator begin(FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).begin();
  }

  RelationVecConstIterator end(FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return fromSetRelationsVec(fromSetIndex).end();
  }

  RelationVecConstIteratorPair range(FromPositionType fromSetIndex) const
  {
    return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
  }
  /// @}

  /// \brief Return the to-set positions associated with fromSetIndex.
  /// \pre 0 <= fromSetIndex < fromSetSize()
  RelationVec const& operator[](FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

  /// \brief Return the number of entries associated with fromSetIndex.
  FlatPositionType size(FromPositionType fromSetIndex) const
  {
    verifyPosition(fromSetIndex);
    return static_cast<FlatPositionType>(fromSetRelationsVec(fromSetIndex).size());
  }

  /// \brief Sum the per-element cardinalities to obtain the total number of entries.
  FlatPositionType totalSize() const
  {
    FlatPositionType sz = 0;
    for(auto& vec : m_relationsVec)
    {
      sz += static_cast<FlatPositionType>(vec.size());
    }
    return sz;
  }

  bool hasFromSet() const { return !policies::EmptySetTraits<FromSetType>::isEmpty(m_fromSet); }
  FromSetType* fromSet() { return m_fromSet; }
  const FromSetType* fromSet() const { return m_fromSet; }

  bool hasToSet() const { return !policies::EmptySetTraits<ToSetType>::isEmpty(m_toSet); }
  ToSetType* toSet() { return m_toSet; }
  const ToSetType* toSet() const { return m_toSet; }

  FromPositionType fromSetSize() const
  {
    return static_cast<FromPositionType>(m_relationsVec.size());
  }

  ToPositionType toSetSize() const { return m_toSet->size(); }

  bool isValid(bool verboseOutput = false) const;

public:  // Modifying functions
  /// \brief Append a to-set position to the entries associated with fromSetIndex.
  /// \pre Both positions identify elements in their respective sets.
  void insert(FromPositionType fromSetIndex, ToPositionType toSetIndex)
  {
    verifyPosition(fromSetIndex);
    m_relationsVec[fromSetIndex].push_back(toSetIndex);
  }

  RelationVec& operator[](FromPositionType fromSetIndex)
  {
    verifyPosition(fromSetIndex);
    return m_relationsVec[fromSetIndex];
  }

public:
  /**
   * \name Direct data access
   * \brief Access the vector of to-set positions for a from-set position.
   * \note Writes must preserve valid to-set positions.
   */

  /// \{

  /**
   * \brief Access the to-set positions associated with fromSetPos.
   *
   * \param fromSetPos Position in the from-set.
   */
  RelationVec& data(FromPositionType fromSetPos)
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /**
   * \brief Access the to-set positions associated with fromSetPos.
   *
   * \param fromSetPos Position in the from-set.
   */
  const RelationVec& data(FromPositionType fromSetPos) const
  {
    verifyPosition(fromSetPos);
    return m_relationsVec[fromSetPos];
  }

  /// \}

private:
  inline void verifyPosition(FromPositionType AXOM_DEBUG_PARAM(fromSetIndex)) const
  {
    SLIC_ASSERT_MSG(
      fromSetIndex >= 0 && fromSetIndex < static_cast<FromPositionType>(m_fromSet->size()),
      "Index " << fromSetIndex << " out of range [0," << m_fromSet->size() << ")");
  }

  inline RelationVec& fromSetRelationsVec(FromPositionType fromSetIndex)
  {
    return m_relationsVec[fromSetIndex];
  }
  inline RelationVec const& fromSetRelationsVec(FromPositionType fromSetIndex) const
  {
    return m_relationsVec[fromSetIndex];
  }

private:
  FromSetType* m_fromSet;
  ToSetType* m_toSet;

  RelationsContainer m_relationsVec;
};

template <typename FirstSetType, typename SecondSetType>
bool DynamicVariableRelation<FirstSetType, SecondSetType>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream sstr;

  if(!hasFromSet() || !hasToSet())
  {
    if(!m_relationsVec.empty())
    {
      if(verboseOutput)
      {
        sstr << "\n\t* relations vector was not empty "
             << " -- fromSet was " << (!hasFromSet() ? "" : " not ") << "null"
             << " , toSet was " << (!hasToSet() ? "" : " not ") << "null";
      }

      bValid = false;
    }
  }
  else
  {
    if(verboseOutput)
    {
      sstr << "\n\t* Neither set was null";
    }

    // Check that the the relations vector has the right size
    // (should be same as fromSet's size() )
    if(static_cast<FromPositionType>(m_relationsVec.size()) != m_fromSet->size())
    {
      if(verboseOutput)
      {
        sstr << "\n\t* relations vector has the wrong size."
             << "\n\t-- from set size is: " << m_fromSet->size()
             << "\n\t-- expected relation size: " << m_fromSet->size()
             << "\n\t-- actual size: " << m_relationsVec.size();
      }
      bValid = false;
    }

    // Check that all elements of the relations vector point to
    // valid  set elements in the toSet
    for(FromPositionType fromIdx = 0; fromIdx < m_fromSet->size(); ++fromIdx)
    {
      FromPositionType idx = fromIdx;
      for(RelationVecConstIterator rIt = begin(idx), rEnd = end(idx); rIt < rEnd; ++rIt)
      {
        if(*rIt >= m_toSet->size())
        {
          if(verboseOutput)
          {
            sstr << "\n\t* relation for element " << m_fromSet->at(fromIdx)
                 << " of fromSet had an out-of-range element.-- value "
                 << std::distance(begin(idx), rIt) << " was " << *rIt
                 << ". Max possible value should be " << m_toSet->size() << ".";
          }
          bValid = false;
        }
      }
    }
  }

  if(verboseOutput)
  {
    std::stringstream sstr2;
    sstr2 << "\n*** Detailed results of isValid on the relation.\n";
    if(bValid)
    {
      sstr2 << "(dynamic,variable) Relation was valid." << std::endl;
    }
    else
    {
      sstr2 << "Relation was NOT valid.\n" << sstr.str() << std::endl;
    }

    if(m_fromSet)
    {
      sstr2 << "\n** fromSet has size " << m_fromSet->size() << ": ";
    }
    if(m_toSet)
    {
      sstr2 << "\n** toSet has size " << m_toSet->size() << ": ";
    }

    if(m_relationsVec.empty())
    {
      sstr2 << "\n** relations vec is empty:";
    }
    else
    {
      FlatPositionType overallCount = 0;
      sstr2 << "\n** relations vec elements:";

      for(FromPositionType fromIdx = 0; fromIdx < m_fromSet->size(); ++fromIdx)
      {
        FromPositionType idx = fromIdx;
        sstr2 << "\n\t" << m_fromSet->at(fromIdx) << " (" << size(idx) << "):\t";
        std::copy(begin(idx), end(idx), std::ostream_iterator<ToPositionType>(sstr2, " "));
        overallCount += size(idx);
      }
      sstr2 << "\n\n\tOverall size of relation" << overallCount << std::endl;

      SLIC_INFO(sstr2.str());
    }
  }

  return bValid;
}

}  // end namespace axom::slam
