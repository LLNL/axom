/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "StaticVariableRelation.hpp"

#include <sstream>
#include <iostream>
#include <iterator>

namespace axom {
namespace slam {

  StaticVariableRelation::StaticVariableRelation (Set* fromSet, Set* toSet)
      : m_fromSet(fromSet), m_toSet(toSet)
  {}

  void StaticVariableRelation::bindRelationData(RelationVec const& beginsVec, RelationVec const& toOffsets)
  {
    m_fromSetBeginsVec.clear();
    m_fromSetBeginsVec.reserve(beginsVec.size());
    std::copy(beginsVec.begin(), beginsVec.end(), std::back_inserter(m_fromSetBeginsVec));

    m_toSetIndicesVec.clear();
    m_toSetIndicesVec.reserve(toOffsets.size());
    std::copy(toOffsets.begin(), toOffsets.end(), std::back_inserter(m_toSetIndicesVec));
  }

  bool StaticVariableRelation::isValid(bool verboseOutput) const
  {
    bool bValid = true;

    std::stringstream errSstr;


    if( *m_fromSet == s_nullSet || *m_toSet == s_nullSet)
    {
      if(!m_fromSetBeginsVec.empty())
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* fromSetBeginsVec was not empty "
                  << " -- fromSet was " << (*m_fromSet == s_nullSet ? "" : " not ") << "null"
                  << " , toSet was " << (*m_toSet == s_nullSet ? "" : " not ") << "null";
        }

        bValid = false;
      }
      if(!m_toSetIndicesVec.empty())
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* toSetIndicesVec was not empty "
                  << " -- fromSet was " << (*m_fromSet == s_nullSet ? "" : " not ") << "null"
                  << " , toSet was " << (*m_toSet == s_nullSet ? "" : " not ") << "null";
        }

        bValid = false;
      }
    }
    else
    {
      if(verboseOutput)
        errSstr << "\n\t* Neither set was null";

      // Check that the fromSet vector has the correct size
      if( static_cast<SetPosition>(m_fromSetBeginsVec.size()) != (m_fromSet->size() + 1) )
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* fromSetBeginsVec was the wrong size."
                  << " -- expected " << (m_fromSet->size() + 1)
                  << ", actual " << m_fromSetBeginsVec.size() << ")";
        }

        bValid = false;
      }

      // Check that no element of the fromSetBegins vector points outside of the toSetIndices vector
      for(RelationVecConstIterator it = m_fromSetBeginsVec.begin(), itEnd = m_fromSetBeginsVec.end(); it != itEnd; ++it)
      {
        if( *it > static_cast<SetPosition>(m_toSetIndicesVec.size()) )
        {
          if(verboseOutput)
          {
            errSstr << "\n\t* fromSetBeginsVec had an out-of-range element."
                    << " -- value of element " << std::distance(m_fromSetBeginsVec.begin(), it) << " was " << *it
                    << ". Max possible value should be " << m_toSetIndicesVec.size() << ".";
          }

          bValid = false;
        }
      }

      // Check that all elements of the toSetIndices vector point to valid set elements
      for(RelationVecConstIterator it = m_toSetIndicesVec.begin(), itEnd = m_toSetIndicesVec.end(); it != itEnd; ++it)
      {
        if( *it >= m_toSet->size() )
        {
          if(verboseOutput)
          {
            errSstr << "\n\t* toSetIndices had an out-of-range element."
                    << " -- value of element " << std::distance(m_toSetIndicesVec.begin(), it) << " was " << *it
                    << ". Max possible value should be " << m_toSet->size() << ".";
          }

          bValid = false;
        }
      }
    }


    if(verboseOutput)
    {
      std::stringstream sstr;

      sstr << "\n*** Detailed results of isValid on the relation.\n";
      if(bValid)
      {
        sstr << "Relation was valid." << std::endl;
      }
      else
      {
        sstr  << "Relation was NOT valid.\n"
              << errSstr.str()
              << std::endl;
      }

      const int fromSz = m_fromSetBeginsVec.size();
      sstr << "\n** fromSetBeginsVec vec w/ size " << fromSz << ": ";
      for(SetPosition i = 0; i< fromSz; ++i)
        sstr << m_fromSetBeginsVec[i] << " ";

      const int toSz = m_toSetIndicesVec.size();
      sstr << "\n** toSetIndices vec w/ size " << toSz << ": ";
      for(SetPosition i = 0; i< toSz; ++i)
        sstr << m_toSetIndicesVec[i] << " ";

      sstr  << "\n** relation data "
            << "(total size: " << m_toSetIndicesVec.size() << "):\n";

      bool hasData = m_fromSet != AXOM_NULLPTR
          && !m_fromSetBeginsVec.empty()
          && !m_toSetIndicesVec.empty();
      if(hasData)
      {
        for(int i = 0; i < m_fromSet->size(); ++i)
        {
          sstr << "\t[ ";
          RelationSet fSet = (*this)[i];
          for(int j = 0; j < fSet.size(); ++j)
          {
            sstr << fSet[j] << " ";
          }
          sstr << "]\n";
        }
      }
      SLIC_DEBUG(sstr.str());

    }

    return bValid;
  }

} // namespace slam
} // namespace axom
