// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_PRIMAL_BVH_MORTON_H_
#define AXOM_PRIMAL_BVH_MORTON_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"    // for fixed bitwidth types

// C/C++ includes
#include <cmath>                  // for fmin()/fmax

namespace axom
{
namespace primal
{
namespace bvh
{

//Returns 30 bit morton code for coordinates for
// x, y, and z are expecting to be between [0,1]
AXOM_HOST_DEVICE
axom::int32 morton32_encode( axom::float32 x,
                             axom::float32 y,
                             axom::float32 z=0.0f );

//Returns 30 bit morton code for coordinates for
//coordinates in the unit cude
AXOM_HOST_DEVICE
axom::int64 morton64_encode( axom::float32 x,
                             axom::float32 y,
                             axom::float32 z=0.0f );

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */

#endif
