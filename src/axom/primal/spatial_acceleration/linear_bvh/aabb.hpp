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

// Forward declarations
template < typename FloatType, int NDIMS >
class AABB
{
  AXOM_STATIC_ASSERT_MSG( (NDIMS >= 2 && NDIMS <=3),
      "The AABB class is supported only in 2D and 3D." );
};

template < typename FloatType, int NDIMS >
std::ostream& operator<<( std::ostream &os,
                          const AABB< FloatType, NDIMS > &aabb);

//------------------------------------------------------------------------------
// 3D AABB Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class AABB< FloatType, 3 >
{

public:
  Range< double > m_x;
  Range< double > m_y;
  Range< double > m_z;

  AXOM_HOST_DEVICE
  void include(const AABB &other)
  {
    m_x.include(other.m_x);
    m_y.include(other.m_y);
    m_z.include(other.m_z);
  }

  AXOM_HOST_DEVICE
  void include(const Vec<FloatType,3>& point)
  {
    m_x.include(point[0]);
    m_y.include(point[1]);
    m_z.include(point[2]);
  }

  AXOM_HOST_DEVICE
  void expand(const FloatType& epsilon)
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
  void scale(const FloatType& scale)
  {
    assert(scale >= 1.f);
    m_x.scale(scale);
    m_y.scale(scale);
    m_z.scale(scale);
  }

  AXOM_HOST_DEVICE
  Vec< FloatType, 3 > center() const
  {
    return make_vec< FloatType, 3 > ( m_x.center(), m_y.center(),m_z.center() );
  }

  std::ostream& print( std::ostream& os ) const
  {
    os << "x-range: " << m_x << "; "
       << "y-range: " << m_y << "; "
       << "z-range: " << m_z;
    return os;
  }

};

//------------------------------------------------------------------------------
// 2D AABB Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class AABB< FloatType, 2 >
{

public:
  Range< double > m_x;
  Range< double > m_y;


  AXOM_HOST_DEVICE
  void include(const AABB &other)
  {
    m_x.include(other.m_x);
    m_y.include(other.m_y);
  }

  AXOM_HOST_DEVICE
  void include(const Vec<FloatType,2>& point)
  {
    m_x.include(point[0]);
    m_y.include(point[1]);
  }

  AXOM_HOST_DEVICE
  void expand(const FloatType& epsilon)
  {
    assert(epsilon > 0.f);
    m_x.include(m_x.min() - epsilon);
    m_x.include(m_x.max() + epsilon);
    m_y.include(m_y.min() - epsilon);
    m_y.include(m_y.max() + epsilon);
  }

  AXOM_HOST_DEVICE
  void scale(const FloatType& scale)
  {
    assert(scale >= 1.f);
    m_x.scale(scale);
    m_y.scale(scale);
  }

  AXOM_HOST_DEVICE
  Vec< FloatType, 3 > center() const
  {
    return make_vec< FloatType, 2 > ( m_x.center(), m_y.center() );
  }

  std::ostream& print( std::ostream& os ) const
  {
    os << "x-range: " << m_x << "; "
       << "y-range: " << m_y;
    return os;
  }
};

//------------------------------------------------------------------------------
// FREE FUNCTIONS
//------------------------------------------------------------------------------
template < typename FloatType, int NDIMS >
std::ostream& operator<<( std::ostream &os,
                          const AABB< FloatType, NDIMS > &aabb)
{
  AXOM_STATIC_ASSERT( NDIMS >=1 && NDIMS <= 3 );
  AXOM_STATIC_ASSERT( std::is_floating_point< FloatType >::value );

  return aabb.print( os );
}


} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
