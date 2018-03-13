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
 * \file IndirectionSet.cpp
 *
 *  \brief Implementation of the IndirectionSet class
 */

#include <stdexcept>
#include <sstream>

#include "IndirectionSet.hpp"

namespace axom
{
namespace slam
{

const NullSet IndirectionSet::s_nullSet;

bool IndirectionSet::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream errStr;

  // Not much to check here since we are essentially wrapping around an
  // array/vector

  if(verboseOutput)
  {
    if( !bValid)
      std::cout << " There was a problem: " << errStr.str() << std::endl;
  }

  return bValid;
}

} /* namespace slam */
} /* namespace axom */
