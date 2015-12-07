/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/**
 *  \file SidreWrapperHelpers.hpp
 *
 *  \brief File used to contain helper functions for C API wrappers.
 *         User code should not include this file.
 *
 */

#ifndef SIDREWRAPPERHELPERS_HPP_
#define SIDREWRAPPERHELPERS_HPP_

#include "Bitwidth_Style_Types.h"

// SiDRe project headers
#include "sidre/SidreTypes.hpp"
#include "sidre/SidreTypes.h"

namespace asctoolkit
{
namespace sidre
{

// This version is used for the C.
inline TypeID getTypeID( const int typeID )
{

  return static_cast<TypeID>(typeID);
}

} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* SIDREWRAPPERHELPERS_HPP_ */
