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
 * \file OrderedSet.cpp
 */

#include "OrderedSet.hpp"

namespace asctoolkit {
namespace slam {

namespace policies {
  const NullSet NoSubset::s_nullSet;
  NullSet VirtualParentSubset::s_nullSet;

} // namespace policies
} // namespace slam
} // namespace asctoolkit
