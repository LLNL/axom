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

//------------------------------------------------------------------------------
/**
 * \brief This class and its specializations provide a type trait that lets us
 *        determine the type that should be used to accumulate values when we
 *        do floating point math.
 *
 * \note this belongs in algorithm utilities, maybe core.
 */
template <typename T>
struct accumulation_traits { using type = float; };

template <>
struct accumulation_traits<double> { using type = double; };

template <>
struct accumulation_traits<long> { using type = double; };

template <>
struct accumulation_traits<unsigned long> { using type = double; };

//-------------------------------------------------------------------------
template <typename T>
AXOM_HOST_DEVICE
int32 bsearch(T value, const axom::ArrayView<T> &view)
{
  int32 index = -1;
  int32 left = 0;
  int32 right = view.size() - 1;
  while(left <= right)
  {
    int32 m = (left + right) / 2;
    if(view[m] < value)
      left = m + 1;
    else if(view[m] > value)
      right = m - 1;
    else
    {
      index = m;
      break;
    }
  }

  return index;
}

//------------------------------------------------------------------------------
/// Based on a Jenkins hash, modified to include length and hash forwards and
/// backwards to make int64 rather than int32.
AXOM_HOST_DEVICE
uint64 hash_bytes(const uint8 *data, uint32 length)
{
  uint32 hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const uint8 *>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  uint32 hashr = hash;
  for(uint32 i = 0; i < length; i++)
  {
    hash += data[i];
    hash += hash << 10;
    hash ^= hash >> 6;

    hashr += data[length - 1 - i];
    hashr += hashr << 10;
    hashr ^= hashr >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  hashr += hashr << 3;
  hashr ^= hashr >> 11;
  hashr += hashr << 15;

  return (static_cast<uint64>(hash) << 32) | hashr;
}

//------------------------------------------------------------------------------
template <typename ValueType, typename IndexType>
AXOM_HOST_DEVICE
uint32
sort_values(ValueType *v, IndexType n)
{
  for(IndexType i = 0; i < n-1; i++)
  {
    const IndexType m = n - i - 1;
    for(IndexType j = 0; j < m; j++)
    {
      if(v[j] > v[j+1])
      {
        axom::utilities::swap(v[j], v[j+1]);
      }
    }
  }
}

//------------------------------------------------------------------------------
template <typename ValueType>
uint64
AXOM_HOST_DEVICE
make_name_1(ValueType id)
{
  return hash_bytes(reinterpret_cast<uint8*>(&id), sizeof(ValueType));
};

//------------------------------------------------------------------------------
template <typename ValueType>
AXOM_HOST_DEVICE
uint64
make_name_2(ValueType id0, ValueType id1)
{
  ValueType data[2] = {id0, id1};
  if(id1 < id0)
  {
    data[0] = id1;
    data[1] = id0;
  }
  return hash_bytes(reinterpret_cast<uint8*>(data), 2 * sizeof(ValueType));
};

//-------------------------------------------------------------------------
template <typename ViewType>
uint64
AXOM_HOST_DEVICE
make_name_n(const ViewType &view, int32 start, int32 n)
{
  using value_type = typename ViewType::value_type;

  if(n == 2)
    return make_name_2(view[start], view[start + 1]);

  value_type v[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // pick largest number of blends.
  for(int32 i = 0; i < n; i++)
  {
    v[i] = view[start + i];
  }
  sort_values(v, n);

  return hash_bytes(reinterpret_cast<uint8*>(v), n * sizeof(value_type));
}

//------------------------------------------------------------------------------
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
