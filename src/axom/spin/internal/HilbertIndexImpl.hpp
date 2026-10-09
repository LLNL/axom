// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! \file HilbertIndexImpl.hpp
 *  \brief Scalar Skilling axes/transpose transform and bit-plane packing.
 */

#include "axom/core/Macros.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <cstdint>
#include <type_traits>

namespace axom::spin::internal::hilbert
{
/// Skilling transform for inputs that satisfy Hilbertizer's preconditions.
template <typename CoordType, typename IndexType, int DIM>
struct HilbertTransform
{
  using PointType = primal::Point<CoordType, DIM>;
  // Use unsigned words of at least 32 bits to avoid integer promotion in shifts.
  using WorkType =
    typename std::conditional<(sizeof(IndexType) <= 4), std::uint32_t, std::uint64_t>::type;

  AXOM_HOST_DEVICE static IndexType toIndex(const PointType& point, int bits)
  {
    WorkType axes[DIM];
    for(int axis = 0; axis < DIM; ++axis)
    {
      axes[axis] = static_cast<WorkType>(point[axis]);
    }

    // Undo orientation changes, from the largest coordinate bit downwards.
    for(int bit = bits - 1; bit > 0; --bit)
    {
      const WorkType q = WorkType {1} << bit;
      const WorkType p = q - 1;
      for(int axis = 0; axis < DIM; ++axis)
      {
        if(axes[axis] & q)
        {
          axes[0] ^= p;
        }
        else
        {
          const WorkType exchange = (axes[0] ^ axes[axis]) & p;
          axes[0] ^= exchange;
          axes[axis] ^= exchange;
        }
      }
    }

    // Gray-encode the axes, then apply the orientation correction.
    for(int axis = 1; axis < DIM; ++axis)
    {
      axes[axis] ^= axes[axis - 1];
    }
    WorkType correction = 0;
    for(int bit = bits - 1; bit > 0; --bit)
    {
      const WorkType q = WorkType {1} << bit;
      if(axes[DIM - 1] & q)
      {
        correction ^= q - 1;
      }
    }
    for(int axis = 0; axis < DIM; ++axis)
    {
      axes[axis] ^= correction;
    }

    // Pack each plane with axis zero in its most significant bit.
    // Shift one bit at a time to support keys that fill the unsigned word.
    WorkType index = 0;
    for(int bit = bits - 1; bit >= 0; --bit)
    {
      for(int axis = 0; axis < DIM; ++axis)
      {
        index = (index << 1) | ((axes[axis] >> bit) & WorkType {1});
      }
    }
    return static_cast<IndexType>(index);
  }

  AXOM_HOST_DEVICE static PointType toPoint(IndexType index, int bits)
  {
    WorkType axes[DIM] = {};
    WorkType remaining = static_cast<WorkType>(index);
    for(int bit = 0; bit < bits; ++bit)
    {
      for(int axis = DIM - 1; axis >= 0; --axis)
      {
        axes[axis] |= (remaining & WorkType {1}) << bit;
        remaining >>= 1;
      }
    }

    // Decode the Gray code. Stop at axis 1 to avoid reading axes[-1].
    const WorkType correction = axes[DIM - 1] >> 1;
    for(int axis = DIM - 1; axis > 0; --axis)
    {
      axes[axis] ^= axes[axis - 1];
    }
    axes[0] ^= correction;

    // Restore orientation in the reverse order of the forward transform.
    for(int bit = 1; bit < bits; ++bit)
    {
      const WorkType q = WorkType {1} << bit;
      const WorkType p = q - 1;
      for(int axis = DIM - 1; axis >= 0; --axis)
      {
        if(axes[axis] & q)
        {
          axes[0] ^= p;
        }
        else
        {
          const WorkType exchange = (axes[0] ^ axes[axis]) & p;
          axes[0] ^= exchange;
          axes[axis] ^= exchange;
        }
      }
    }

    PointType point;
    for(int axis = 0; axis < DIM; ++axis)
    {
      point[axis] = static_cast<CoordType>(axes[axis]);
    }
    return point;
  }
};
}  // namespace axom::spin::internal::hilbert
