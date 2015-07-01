/**
 * \file StaticVariableRelation.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */

#include "StaticVariableRelation.hpp"

#include <sstream>
#include <iostream>
#include <iterator>

namespace asctoolkit {
namespace meshapi {




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

      sstr << "\n** fromSetBeginsVec vec w/ size " << m_fromSetBeginsVec.size() << ": ";
      std:: copy( m_fromSetBeginsVec.begin(), m_fromSetBeginsVec.end(), std::ostream_iterator<SetPosition>(sstr, " "));

      sstr << "\n** toSetIndices vec w/ size " << m_toSetIndicesVec.size() << ": ";
      std:: copy( m_toSetIndicesVec.begin(),  m_toSetIndicesVec.end(),  std::ostream_iterator<SetPosition>(sstr, " "));

      std::cout << sstr.str() << std::endl;

    }

    return bValid;
  }

} // namespace meshapi
} // namespace asctoolkit
