/*
 * OrderedSet.cpp
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include <stdexcept>
#include <sstream>
#include "OrderedSet.hpp"

namespace asctoolkit {
namespace meshapi {

OrderedSet::Index OrderedSet::at( Index idx )
{
  if(idx >= size())
  {
    std::stringstream sstr;
    sstr<< "MeshAPI::OrderedSet -- requested out of range element at position "
        << idx << ", but set only has " << size() << " elements.";
    throw std::out_of_range(sstr.str());
  }
  return idx;
}


} /* namespace meshapi */
} /* namespace asctoolkit */

