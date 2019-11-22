// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __ZOO_CLIPPING_TABLES_H__
#define __ZOO_CLIPPING_TABLES_H__

/**
 * \file ZooClippingTables.hpp
 * 
 * \brief Contains the definitions for the clipping cases and enumerator
 *        for the shape types in the zoo.
 */

#include <vector>

namespace axom
{
namespace mir
{
  extern const int quadClipTable[16][19];
  extern const int triangleClipTable[8][10];
  
  extern const std::vector<std::vector<int> > triangleClipTableVec;
  extern const std::vector<std::vector<int> > quadClipTableVec;
}
}
#endif
