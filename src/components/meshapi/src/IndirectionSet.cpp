/**
 * \file IndirectionSet.cpp
 *
 *  \brief Implementation of the IndirectionSet class
 */

#include <stdexcept>
#include <sstream>

#include "IndirectionSet.hpp"

namespace asctoolkit {
namespace meshapi {

  const NullSet IndirectionSet::s_nullSet;

  bool IndirectionSet::isValid(bool verboseOutput) const
  {
    bool bValid = true;

    std::stringstream errStr;

    // Not much to check here since we are essentially wrapping around an array/vector

    if(verboseOutput)
    {
      if( !bValid)
        std::cout << " There was a problem: " << errStr.str() << std::endl;
    }

    return bValid;
  }

} /* namespace meshapi */
} /* namespace asctoolkit */
