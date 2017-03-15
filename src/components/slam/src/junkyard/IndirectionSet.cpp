/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file IndirectionSet.cpp
 *
 *  \brief Implementation of the IndirectionSet class
 */

#include <stdexcept>
#include <sstream>

#include "IndirectionSet.hpp"

namespace axom {
namespace slam {

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

} /* namespace slam */
} /* namespace axom */
