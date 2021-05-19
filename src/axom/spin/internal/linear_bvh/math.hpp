// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_MATH_H_
#define AXOM_SPIN_BVH_MATH_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"  // for fixed bitwidth types

// include math so we can use functions defined
// in both cuda and c
#include <math.h>

#ifndef __CUDACC__
  // make sure min / max resolve for both cuda and cpu
  #include <math.h>
  #include <string.h>  //resolve memcpy
using namespace std;
#endif

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
constexpr axom::int32 AXOM_INF_32 = 0x7f800000U;
constexpr axom::int32 AXOM_NG_INF_32 = 0xff800000U;

constexpr axom::int64 AXOM_INF_64 = 0x7ff0000000000000ULL;
constexpr axom::int64 AXOM_NG_INF_64 = 0xfff0000000000000ULL;

constexpr axom::int32 AXOM_NAN_32 = 0x7FC00000U;
constexpr axom::int64 AXOM_NAN_64 = 0x7FF8000000000000ULL;

constexpr float AXOM_EPSILON_32 = 1e-5f;
constexpr double AXOM_EPSILON_64 = 1e-9f;

union Bits32
{
  axom::float32 scalar;
  axom::int32 bits;
};

union Bits64
{
  axom::float64 scalar;
  axom::int64 bits;
};

template <typename T>
inline AXOM_HOST_DEVICE T epsilon()
{
  return 1;
}

template <>
inline AXOM_HOST_DEVICE axom::float32 epsilon<axom::float32>()
{
  return AXOM_EPSILON_32;
}

template <>
inline AXOM_HOST_DEVICE axom::float64 epsilon<axom::float64>()
{
  return AXOM_EPSILON_64;
}

inline AXOM_HOST_DEVICE axom::float32 nan32()
{
  Bits32 nan;
  nan.bits = AXOM_NAN_32;
  return nan.scalar;
}

inline AXOM_HOST_DEVICE axom::float32 infinity32()
{
  Bits32 inf;
  inf.bits = AXOM_INF_32;
  return inf.scalar;
}

inline AXOM_HOST_DEVICE axom::float32 neg_infinity32()
{
  Bits32 ninf;
  ninf.bits = AXOM_NG_INF_32;
  return ninf.scalar;
}

inline AXOM_HOST_DEVICE axom::float64 nan64()
{
  Bits64 nan;
  nan.bits = AXOM_NAN_64;
  return nan.scalar;
}

inline AXOM_HOST_DEVICE axom::float64 infinity64()
{
  Bits64 inf;
  inf.bits = AXOM_INF_64;
  return inf.scalar;
}

inline AXOM_HOST_DEVICE axom::float64 neg_infinity64()
{
  Bits64 ninf;
  ninf.bits = AXOM_NG_INF_64;
  return ninf.scalar;
}

template <typename T>
inline AXOM_HOST_DEVICE T infinity();

template <>
inline AXOM_HOST_DEVICE axom::float32 infinity<axom::float32>()
{
  return infinity32();
}

template <>
inline AXOM_HOST_DEVICE axom::float64 infinity<axom::float64>()
{
  return infinity64();
}

//
// count leading zeros
//
inline AXOM_HOST_DEVICE axom::int32 clz(axom::int32 x)
{
  axom::int32 y;
  axom::int32 n = 32;
  y = x >> 16;
  if(y != 0)
  {
    n = n - 16;
    x = y;
  }
  y = x >> 8;
  if(y != 0)
  {
    n = n - 8;
    x = y;
  }
  y = x >> 4;
  if(y != 0)
  {
    n = n - 4;
    x = y;
  }
  y = x >> 2;
  if(y != 0)
  {
    n = n - 2;
    x = y;
  }
  y = x >> 1;
  if(y != 0) return axom::int32(n - 2);
  return axom::int32(n - x);
}

inline AXOM_HOST_DEVICE axom::float64 pi()
{
  return 3.14159265358979323846264338327950288;
}

inline AXOM_HOST_DEVICE axom::float32 rcp(axom::float32 f) { return 1.0f / f; }

inline AXOM_HOST_DEVICE axom::float64 rcp(axom::float64 f) { return 1.0 / f; }

inline AXOM_HOST_DEVICE axom::float64 rcp_safe(axom::float64 f)
{
  return rcp((fabs(f) < 1e-8) ? 1e-8 : f);
}

inline AXOM_HOST_DEVICE axom::float32 rcp_safe(axom::float32 f)
{
  return rcp((fabs(f) < 1e-8f) ? 1e-8f : f);
}

template <typename T>
inline AXOM_HOST_DEVICE T clamp(const T &val, const T &min_val, const T &max_val)
{
  return min(max_val, max(min_val, val));
}

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */
#endif
