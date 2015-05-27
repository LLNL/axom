/*
 * Set.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <stdexcept>
#include <sstream>

#include "meshapi/Set.hpp"

namespace asctoolkit {
namespace meshapi {

IndirectionSet::Index& Set::at( size_type idx )
{
  if(idx >= m_entities.size())
  {
    std::stringstream sstr;
    sstr<< "MeshAPI::Set -- requested out of range element at position "
        << idx << ", but set only has " << m_entities.size() << " elements.";
    throw std::out_of_range(sstr.str());
  }
  return m_entities[idx];
}

IndirectionSet::Index const& Set::at( size_type idx ) const
{
  if(idx >= m_entities.size())
  {
    std::stringstream sstr;
    sstr<< "MeshAPI::Set -- requested out of range element at position "
        << idx << ", but set only has " << m_entities.size() << " elements.";
    throw std::out_of_range(sstr.str());
  }
  return m_entities[idx];
}



} /* namespace meshapi */
} /* namespace asctoolkit */

