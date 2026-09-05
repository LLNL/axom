// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file SizeChecks.hpp
 * \brief Internal checks for SLAM size and component-count arithmetic.
 */

#include "axom/core/Types.hpp"
#include "axom/core/Macros.hpp"
#include "axom/slic.hpp"

#include <limits>
#include <utility>

namespace axom::slam::detail
{
/// Whether a nonnegative product fits, without evaluating an overflowing product.
template <typename Int>
AXOM_HOST_DEVICE constexpr bool nonnegativeProductFits(Int a, Int b)
{
  return a >= 0 && b >= 0 && (a == 0 || b <= std::numeric_limits<Int>::max() / a);
}

/// Validate the count in both the map's indexing type and ArrayView's size type.
template <typename Position, typename Stride>
bool mapStorageSize(Position size, Stride stride, Position& count)
{
  if(stride <= 0 || !std::in_range<Position>(stride) ||
     !nonnegativeProductFits(size, static_cast<Position>(stride)))
  {
    return false;
  }
  count = size * static_cast<Position>(stride);
  return std::in_range<axom::IndexType>(count);
}

template <typename Position, typename Stride>
Position checkedMapStorageSize(Position size, Stride stride)
{
  Position count {};
  SLIC_ERROR_IF(!mapStorageSize(size, stride, count),
                "SLAM map requires a nonnegative size, positive component count, and "
                "representable storage size; received size "
                  << size << " and components " << stride);
  return count;
}

template <typename Position, typename Stride>
Position checkedMapStride(Stride stride)
{
  SLIC_ERROR_IF(stride <= 0 || !std::in_range<Position>(stride),
                "SLAM map component count must be positive and representable; received " << stride);
  return static_cast<Position>(stride);
}
}  // namespace axom::slam::detail
