/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef AXOM_PRIMAL_BVH_RANGE_H_
#define AXOM_PRIMAL_BVH_RANGE_H_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"    // for fixed bitwidth types
#include "axom/primal/spatial_acceleration/linear_bvh/math.hpp"

#include <iostream>

namespace axom
{
namespace primal
{
namespace bvh
{

class Range
{
protected:
  axom::common::float32 m_min = infinity32();
  axom::common::float32 m_max = neg_infinity32();
public:


  AXOM_HOST_DEVICE
  axom::common::float32 min() const
  {
    return m_min;
  }

  AXOM_HOST_DEVICE
  axom::common::float32 max() const
  {
    return m_max;
  }

  AXOM_HOST_DEVICE
  bool is_empty() const
  {
    return m_min > m_max;
  }

  template<typename T>
  AXOM_HOST_DEVICE
  void include(const T &val)
  {
    m_min = fmin(m_min, axom::common::float32(val));
    m_max = fmax(m_max, axom::common::float32(val));
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
  axom::common::float32 center() const
  {
    if(is_empty())
    {
      return nan32();
    }
    else return 0.5f * (m_min + m_max);
  }

  AXOM_HOST_DEVICE
  axom::common::float32 length() const
  {
    if(is_empty())
    {
      return nan32();
    }
    else return m_max - m_min;
  }

  AXOM_HOST_DEVICE
  void scale(axom::common::float32 scale)
  {
    if(is_empty())
    {
      return;
    }

    axom::common::float32 c = center();
    axom::common::float32 delta = scale * 0.5f * length();
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


  friend std::ostream& operator<<(std::ostream &os, const Range &range);

};

inline std::ostream& operator<<(std::ostream &os, const Range &range)
{
  os<<"[";
  os<<range.min()<<", ";
  os<<range.max();
  os<<"]";
  return os;
}

} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
