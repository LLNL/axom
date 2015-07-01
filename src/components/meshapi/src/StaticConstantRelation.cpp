/**
 * \file StaticConstantRelation.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */

#include "StaticConstantRelation.hpp"

#include <sstream>
#include <iostream>
#include <iterator>

namespace asctoolkit {
namespace meshapi {

  StaticConstantRelation::StaticConstantRelation (Set* fromSet, Set* toSet)
      : m_stride(SetPosition()), m_fromSet(fromSet), m_toSet(toSet)
  {}

  void StaticConstantRelation::bindRelationData(const RelationVec & toOffsets, const SetPosition stride)
  {
    m_stride = stride;

    m_toSetIndicesVec.clear();
    m_toSetIndicesVec.reserve(toOffsets.size());
    std::copy(toOffsets.begin(), toOffsets.end(), std::back_inserter(m_toSetIndicesVec));
  }

  bool StaticConstantRelation::isValid(bool verboseOutput) const
  {
    bool bValid = true;

    std::stringstream errSstr;

    if( *m_fromSet == s_nullSet || *m_toSet == s_nullSet)
    {
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

      // Check that the toSetIndices vector has the right size
      if( static_cast<SetPosition>(m_toSetIndicesVec.size()) != (m_stride * m_fromSet->size()) )
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* toSetIndices has the wrong size."
                  << "\n\t-- from set size is: " << m_fromSet->size()
                  << "\n\t-- constant stride is: " << m_stride
                  << "\n\t-- expected relation size: " << (m_stride * m_fromSet->size())
                  << "\n\t-- actual size: " << m_toSetIndicesVec.size()
          ;
        }
        bValid = false;
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

      if(bValid)
      {
        sstr << "(static,constant) Relation with stride " << m_stride << " was valid." << std::endl;
      }
      else
      {
        sstr  << "Relation was NOT valid.\n"
              << errSstr.str()
              << std::endl;
      }

      sstr << "\n*** Detailed results of isValid on the relation.\n";
      if(m_fromSet)
        sstr << "\n** fromSet has size " << m_fromSet->size() << ": ";
      if(m_toSet)
        sstr << "\n** toSet has size " << m_toSet->size() << ": ";

      sstr << "\n** toSetIndices vec w/ size " << m_toSetIndicesVec.size() << ": ";
      std::copy(m_toSetIndicesVec.begin(), m_toSetIndicesVec.end(), std::ostream_iterator<SetPosition>(sstr, " "));

      std::cout << sstr.str() << std::endl;

    }

    return bValid;
  }

} // namespace meshapi
} // namespace asctoolkit
