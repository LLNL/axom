// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __ZOO_CLIPPING_TABLES_H__
#define __ZOO_CLIPPING_TABLES_H__

/**
 * \file ZooClippingTables.hpp
 * 
 * \brief Contains the defintions for the clipping cases and enumerator
 *        for the shape types in the zoo.
 */

namespace axom
{
namespace mir
{

  enum Shape
  {
    Triangle,
    Quad,
    Tetrahedron,
    Triangular_Prism,
    Pyramid,
    Hexahedron
  };

  extern const int quadClipTable[16][19];
  extern const int triangleClipTable[8][10];
}
}
#endif