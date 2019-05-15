// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_PRIMAL_BVH_VEC_H_
#define AXOM_PRIMAL_BVH_VEC_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"    // for fixed bitwidth types

#include <assert.h>
#include <iostream>

#ifdef __CUDACC__
#include "math.h"
#else
#include <math.h>
#endif

namespace axom
{
namespace primal
{
namespace bvh
{

template<typename T, axom::int32 S>
class Vec
{

public:
  T m_data[S];

  //
  //  No contructors so this is a POD type
  //

  inline AXOM_HOST_DEVICE bool operator==(const Vec<T,S> &other) const
  {
    bool e = true;
    for(int i = 0; i < S; ++i)
    {
      if(m_data[i] != other[i]) e = false;
    }
    return e;
  }

  template<typename TT, axom::int32 SS>
  friend std::ostream& operator<<(std::ostream &os, const Vec<TT,SS> &vec);

  inline AXOM_HOST_DEVICE void operator=(const Vec<T,S> &other)
  {
    for(int i = 0; i < S; ++i)
    {
      m_data[i] = other.m_data[i];
    }
  }

  inline AXOM_HOST_DEVICE const T& operator[](const axom::int32 &i) const
  {
    assert(i > -1 && i < S);
    return m_data[i];
  }

  inline AXOM_HOST_DEVICE T& operator[](const axom::int32 &i)
  {
    assert(i > -1 && i < S);
    return m_data[i];
  }

  // scalar mult /  div
  inline AXOM_HOST_DEVICE Vec<T,S> operator*(const T &s) const
  {
    Vec<T,S> res;

    for(int i = 0; i < S; ++i)
    {
      res[i] = m_data[i] * s;
    }

    return res;
  }

  inline AXOM_HOST_DEVICE Vec<T,S> operator/(const T &s) const
  {
    Vec<T,S> res;

    for(int i = 0; i < S; ++i)
    {
      res[i] = m_data[i] / s;
    }

    return res;
  }

  inline AXOM_HOST_DEVICE void operator*=(const T &s)
  {
    for(int i = 0; i < S; ++i)
    {
      m_data[i] *= s;
    }
  }

  inline AXOM_HOST_DEVICE void operator/=(const T &s)
  {
    for(int i = 0; i < S; ++i)
    {
      m_data[i] /= s;
    }

  }

  // vector add / sub

  inline AXOM_HOST_DEVICE Vec<T,S> operator+(const Vec<T,S> &other) const
  {
    Vec<T,S> res;

    for(int i = 0; i < S; ++i)
    {
      res[i] = m_data[i] + other[i];
    }

    return res;
  }

  inline AXOM_HOST_DEVICE Vec<T,S> operator-(const Vec<T,S> &other) const
  {
    Vec<T,S> res;

    for(int i = 0; i < S; ++i)
    {
      res[i] = m_data[i] - other[i];
    }

    return res;
  }

  inline AXOM_HOST_DEVICE void operator+=(const Vec<T,S> &other)
  {

    for(int i = 0; i < S; ++i)
    {
      m_data[i] += other[i];
    }

  }

  inline AXOM_HOST_DEVICE void operator-=(const Vec<T,S> &other)
  {

    for(int i = 0; i < S; ++i)
    {
      m_data[i] -= other[i];
    }

  }

  inline AXOM_HOST_DEVICE Vec<T,S> operator-(void) const
  {
    Vec<T,S> res;

    for(int i = 0; i < S; ++i)
    {
      res[i] = -m_data[i];
    }

    return res;
  }

  inline AXOM_HOST_DEVICE T magnitude() const
  {
    T sum = T(0);

    for(int i = 0; i < S; ++i)
    {
      sum += m_data[i] * m_data[i];
    }

    return sqrtf(sum);
  }


  inline AXOM_HOST_DEVICE void normalize()
  {
    T mag = magnitude();
    *this /= mag;
  }

};

// vector utility functions
template<typename T, axom::int32 S>
inline AXOM_HOST_DEVICE T dot(const Vec<T,S> &a, const Vec<T,S> &b)
{
  T res = T(0);

  for(int i = 0; i < S; ++i)
  {
    res += a[i] * b[i];
  }

  return res;
}

template<typename T>
inline AXOM_HOST_DEVICE Vec<T,3> cross(const Vec<T,3> &a, const Vec<T,3> &b)
{
  Vec<T,3> res;
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
  return res;
}

template<typename TT, axom::int32 SS>
std::ostream& operator<<(std::ostream &os, const Vec<TT,SS> &vec)
{
  os<<"[";
  for(int i = 0; i < SS; ++i)
  {
    os<<vec[i];
    if(i != SS - 1) os<<", ";
  }
  os<<"]";
  return os;
}

// typedefs
using Vec2i  = Vec< axom::int32,2 >;
using Vec2li = Vec< axom::int64,2 >;
using Vec2f  = Vec< axom::float32,2 >;
using Vec2d  = Vec< axom::float64,2 >;

using Vec3i  = Vec< axom::int32,3 >;
using Vec3li = Vec< axom::int64,3 >;
using Vec3f  = Vec< axom::float32,3 >;
using Vec3d  = Vec< axom::float64,3 >;

using Vec4i   = Vec< axom::int32,4 >;
using Vec4li  = Vec< axom::int64,4 >;
using Vec4f   = Vec< axom::float32,4 >;
using Vec4d   = Vec< axom::float64,4 >;

inline AXOM_HOST_DEVICE
Vec2i make_vec2i( const axom::int32 &a,
                  const axom::int32 &b  )
{
  Vec2i res;
  res[0] = a;
  res[1] = b;
  return res;
}

inline AXOM_HOST_DEVICE
Vec2li make_vec2li( const axom::int64 &a,
                    const axom::int64 &b  )
{
  Vec2li res;
  res[0] = a;
  res[1] = b;
  return res;
}

inline AXOM_HOST_DEVICE
Vec2f make_vec2f( const axom::float32 &a,
                  const axom::float32 &b  )
{
  Vec2f res;
  res[0] = a;
  res[1] = b;
  return res;
}

inline AXOM_HOST_DEVICE
Vec2d make_vec2d( const axom::float64 &a,
                  const axom::float64 &b  )
{
  Vec2d res;
  res[0] = a;
  res[1] = b;
  return res;
}

inline AXOM_HOST_DEVICE
Vec3i make_vec3i( const axom::int32 &a,
                  const axom::int32 &b,
                  const axom::int32 &c  )
{
  Vec3i res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  return res;
}

inline AXOM_HOST_DEVICE
Vec3li make_vec3li( const axom::int64 &a,
                    const axom::int64 &b,
                    const axom::int64 &c  )
{
  Vec3li res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  return res;
}

inline AXOM_HOST_DEVICE
Vec3f make_vec3f( const axom::float32 &a,
                  const axom::float32 &b,
                  const axom::float32 &c  )
{
  Vec3f res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  return res;
}

inline AXOM_HOST_DEVICE
Vec3d make_vec3d( const axom::float64 &a,
                  const axom::float64 &b,
                  const axom::float64 &c  )
{
  Vec3d res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  return res;
}

inline AXOM_HOST_DEVICE
Vec4i make_vec4i( const axom::int32 &a,
                  const axom::int32 &b,
                  const axom::int32 &c,
                  const axom::int32 &d  )
{
  Vec4i res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  res[3] = d;
  return res;
}

inline AXOM_HOST_DEVICE
Vec4li make_vec4li( const axom::int64 &a,
                    const axom::int64 &b,
                    const axom::int64 &c,
                    const axom::int64 &d  )
{
  Vec4li res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  res[3] = d;
  return res;
}

inline AXOM_HOST_DEVICE
Vec4f make_vec4f( const axom::float32 &a,
                  const axom::float32 &b,
                  const axom::float32 &c,
                  const axom::float32 &d  )
{
  Vec4f res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  res[3] = d;
  return res;
}

inline AXOM_HOST_DEVICE
Vec4d make_vec4d( const axom::float64 &a,
                  const axom::float64 &b,
                  const axom::float64 &c,
                  const axom::float64 &d  )
{
  Vec4d res;
  res[0] = a;
  res[1] = b;
  res[2] = c;
  res[3] = d;
  return res;
}

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */

#endif
