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

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{


//expands 10-bit unsigned int into 30 bits
inline AXOM_HOST_DEVICE
axom::int32 expand_bits32(axom::int32 x32)
{
  x32 = (x32 | (x32 << 16)) & 0x030000FF;
  x32 = (x32 | (x32 << 8)) & 0x0300F00F;
  x32 = (x32 | (x32 << 4)) & 0x030C30C3;
  x32 = (x32 | (x32 << 2)) & 0x09249249;
  return x32;
}

inline AXOM_HOST_DEVICE
axom::int64 expand_bits64(axom::int32 x)
{
  axom::int64 x64 = x & 0x1FFFFF;
  x64 = (x64 | x64 << 32) & 0x1F00000000FFFF;
  x64 = (x64 | x64 << 16) & 0x1F0000FF0000FF;
  x64 = (x64 | x64 << 8) & 0x100F00F00F00F00F;
  x64 = (x64 | x64 << 4) & 0x10c30c30c30c30c3;
  x64 = (x64 | x64 << 2) & 0x1249249249249249;

  return x64;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// MORTON FUNCTIONS
//------------------------------------------------------------------------------

//Returns 30 bit morton code for coordinates for
// x, y, and z are expecting to be between [0,1]
inline AXOM_HOST_DEVICE
axom::int32 morton32_encode( axom::float32 x,
                             axom::float32 y,
                             axom::float32 z=0.0f )
{
  //take the first 10 bits. Note, 2^10 = 1024
  x = fmin(fmax(x * 1024.0f, 0.0f), 1023.0f);
  y = fmin(fmax(y * 1024.0f, 0.0f), 1023.0f);
  z = fmin(fmax(z * 1024.0f, 0.0f), 1023.0f);

  //expand 10 bits to 30
  axom::int32 xx = expand_bits32((axom::int32)x);
  axom::int32 yy = expand_bits32((axom::int32)y);
  axom::int32 zz = expand_bits32((axom::int32)z);
  //interleave coordinates
  return (zz << 2 | yy << 1 | xx);
}

//Returns 30 bit morton code for coordinates for
//coordinates in the unit cude
inline AXOM_HOST_DEVICE
axom::int64 morton64_encode( axom::float32 x,
                             axom::float32 y,
                             axom::float32 z=0.0f )
{
  //take the first 21 bits. Note, 2^21= 2097152.0f
  x = fmin(fmax(x * 2097152.0f, 0.0f), 2097151.0f);
  y = fmin(fmax(y * 2097152.0f, 0.0f), 2097151.0f);
  z = fmin(fmax(z * 2097152.0f, 0.0f), 2097151.0f);

  //expand the 10 bits to 30
  axom::int64 xx = expand_bits64((axom::int32)x);
  axom::int64 yy = expand_bits64((axom::int32)y);
  axom::int64 zz = expand_bits64((axom::int32)z);

  //interleave coordinates
  return (zz << 2 | yy << 1 | xx);
}

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
