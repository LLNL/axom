// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UTILITIES_HPP_
#define AXOM_MIR_UTILITIES_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/memory_management.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#include <RAJA/RAJA.hpp>

/// This header is for device utility functions for MIR.

// Q: does this belong in core with a better name?

namespace axom
{
namespace mir
{
namespace utilities
{

/**
 * \brief This function makes a unique array of values from an input list of keys.
 *
 * \tparam ExecSpace The execution space.
 * \tparam KeyType   The data type for the keys.
 *
 * \param keys_orig_view The input view that contains the input keys to be made unique.
 * \param[out] skeys     A sorted unique array of keys produced from keys_orig_view.
 * \param[out] sindices  An array of indices that indicate where in the original view the keys came from.
 *
 * \note This code is adapted from Ascent/DevilRay.
 *
 */
template <typename ExecSpace, typename KeyType>
void unique(const axom::ArrayView<KeyType> &keys_orig_view, axom::Array<KeyType> &skeys, axom::Array<axom::IndexType> &sindices)
{
  using for_policy = typename axom::execution_space<ExecSpace>::for_policy;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Make a copy of the keys and make original indices.
  const auto n = keys_orig_view.size();
  axom::Array<KeyType> keys(n, n, allocatorID);
  axom::Array<axom::IndexType> indices(n, n, allocatorID);
  auto keys_view = keys.view();
  auto indices_view = indices.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    keys_view[i] = keys_orig_view[i];
    indices_view[i] = i;
  });

  // Sort the keys, indices in place.
  RAJA::sort_pairs<for_policy>(RAJA::make_span(keys_view, n),
                               RAJA::make_span(indices_view, n));

  // Make a mask array for where differences occur.
  axom::Array<axom::IndexType> mask(n, n, allocatorID);
  auto mask_view = mask.view();
  RAJA::ReduceSum<reduce_policy, axom::IndexType> mask_sum(0);
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    const axom::IndexType different = (keys_view[i] != keys_view[i - 1]) ? 1 : 0;
    const axom::IndexType m = (i >= 1) ? different : 1;
    mask_view[i] = m;
    mask_sum += m;
  });

  // Do a scan on the mask array to build an offset array.
  axom::Array<axom::IndexType> offsets(n, n, allocatorID);
  auto offsets_view = offsets.view();
  RAJA::exclusive_scan<for_policy>(RAJA::make_span(mask_view, n),
                                   RAJA::make_span(offsets_view, n),
                                   RAJA::operators::plus<axom::IndexType>{});

  // Allocate the output arrays.
  const axom::IndexType newsize = mask_sum.get();
  skeys = axom::Array<KeyType>(newsize, newsize, allocatorID);
  sindices = axom::Array<axom::IndexType>(newsize, newsize, allocatorID);

  // Iterate over the mask/offsets to store values at the right
  // offset in the new array.
  auto skeys_view = skeys.view();
  auto sindices_view = sindices.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    if(mask_view[i])
    {
      skeys_view[offsets_view[i]] = keys_view[i];
      sindices_view[offsets_view[i]] = indices_view[i];
    }
  });
}

} // end namespace utilities
} // end namespace mir
} // end namespace axom

#endif
