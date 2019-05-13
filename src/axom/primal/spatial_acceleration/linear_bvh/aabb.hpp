// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_PRIMAL_BVH_AABB_H_
#define AXOM_PRIMAL_BVH_AABB_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"    // for fixed bitwidth types
#include "axom/primal/spatial_acceleration/linear_bvh/range.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/vec.hpp"

#include <iostream>

namespace axom
{
namespace primal
{
namespace bvh
{

class AABB
{

public:
  Range m_x;
  Range m_y;
  Range m_z;

  AXOM_HOST_DEVICE
  void include(const AABB &other)
  {
    m_x.include(other.m_x);
    m_y.include(other.m_y);
    m_z.include(other.m_z);
  }

  AXOM_HOST_DEVICE
  void include(const Vec3f &point)
  {
    m_x.include(point[0]);
    m_y.include(point[1]);
    m_z.include(point[2]);
  }

  AXOM_HOST_DEVICE
  void expand(const axom::float32 &epsilon)
  {
    assert(epsilon > 0.f);
    m_x.include(m_x.min() - epsilon);
    m_x.include(m_x.max() + epsilon);
    m_y.include(m_y.min() - epsilon);
    m_y.include(m_y.max() + epsilon);
    m_z.include(m_z.min() - epsilon);
    m_z.include(m_z.max() + epsilon);
  }

  AXOM_HOST_DEVICE
  void scale(const axom::float32 &scale)
  {
    assert(scale >= 1.f);
    m_x.scale(scale);
    m_y.scale(scale);
    m_z.scale(scale);
  }

  AXOM_HOST_DEVICE
  Vec3f center() const
  {
    return make_vec3f(m_x.center(),
                      m_y.center(),
                      m_z.center());
  }

  //AXOM_HOST_DEVICE
  //AABB identity() const
  //{
  //  return AABB();
  //}

  //AXOM_HOST_DEVICE
  //AABB operator+(const AABB &other) const
  //{
  //  AABB res = *this;
  //  res.include(other);
  //  return res;
  //}

  friend std::ostream& operator<<(std::ostream &os, const AABB &aabb);
};

inline std::ostream& operator<<(std::ostream &os, const AABB &aabb)
{
  os<<"[";
  os<<aabb.m_x.min()<<", ";
  os<<aabb.m_y.min()<<", ";
  os<<aabb.m_z.min();
  os<<"] - [";
  os<<aabb.m_x.max()<<", ";
  os<<aabb.m_y.max()<<", ";
  os<<aabb.m_z.max();
  os<<"]";
  return os;
}

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
