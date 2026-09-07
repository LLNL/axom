// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file Ranges.hpp
 *
 * \brief Host-only integration between SLAM types and the C++20 ranges library.
 *
 * This header is intentionally excluded from axom/slam.hpp so device-facing translation units
 * do not acquire a dependency on the host standard library's ranges implementation.
 */

#include "axom/slam/RangeSet.hpp"

#include <ranges>

namespace std::ranges
{
/**
 * \brief Range iterators can outlive sets with no indirection or parent binding.
 *
 * The iterator stores the range's offset, stride and size by value. This covers
 * both RangeSet and PositionSet. Indirection-backed sets and parent-bound subsets
 * have additional lifetime requirements and are not marked as borrowed here.
 */
template <typename PositionType, typename ElementType, typename OffsetPolicy, typename StridingPolicy, typename InterfacePolicy>
inline constexpr bool enable_borrowed_range<
  axom::slam::GenericRangeSet<PositionType,
                              ElementType,
                              OffsetPolicy,
                              StridingPolicy,
                              axom::slam::policies::NoIndirection<PositionType, ElementType>,
                              axom::slam::policies::NoSubset,
                              InterfacePolicy>> = true;
}  // namespace std::ranges
