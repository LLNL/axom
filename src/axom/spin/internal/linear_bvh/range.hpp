// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_SPIN_BVH_RANGE_H_
#define AXOM_SPIN_BVH_RANGE_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"    // for fixed bitwidth types
#include "axom/spin/internal/linear_bvh/math.hpp"

#include <iostream>

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{

// Forward declarations
template < typename FloatType >
class Range;

template < typename FloatType >
std::ostream& operator<<( std::ostream &os,
                          const Range< FloatType >& range);

template < typename FloatType >
class Range
{
  AXOM_STATIC_ASSERT( std::is_floating_point< FloatType >::value );

public:


  AXOM_HOST_DEVICE
  FloatType min() const
  {
    return m_min;
  }

  AXOM_HOST_DEVICE
  FloatType max() const
  {
    return m_max;
  }

  AXOM_HOST_DEVICE
  bool is_empty() const
  {
    return m_min > m_max;
  }

  AXOM_HOST_DEVICE
  void include(const FloatType& val)
  {
    m_min = fmin( m_min, FloatType(val) );
    m_max = fmax( m_max, FloatType(val) );
  }

  AXOM_HOST_DEVICE
  void include(const Range &other)
  {
    if(!other.is_empty())
    {
      include(other.min());
      include(other.max());
    }
  }

  AXOM_HOST_DEVICE
  Range identity() const
  {
    return Range();
  }

  AXOM_HOST_DEVICE
  FloatType center() const
  {
    if(is_empty())
    {
      return nan32();
    }
    else return 0.5f * (m_min + m_max);
  }

  AXOM_HOST_DEVICE
  FloatType length() const
  {
    if(is_empty())
    {
      return nan32();
    }
    else return m_max - m_min;
  }

  AXOM_HOST_DEVICE
  void scale(FloatType scale)
  {
    if(is_empty())
    {
      return;
    }

    FloatType c = center();
    FloatType delta = scale * 0.5f * length();
    include(c - delta);
    include(c + delta);
  }


  AXOM_HOST_DEVICE
  Range operator+(const Range &other) const
  {
    Range res;
    res.include(*this);
    res.include(other);
    return res;
  }

private:
  FloatType m_min = infinity32();
  FloatType m_max = neg_infinity32();

};

template < typename FloatType >
std::ostream& operator<<( std::ostream &os,
                          const Range< FloatType >& range)
{
  os<<"[";
  os<<range.min()<<", ";
  os<<range.max();
  os<<"]";
  return os;
}


} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif
