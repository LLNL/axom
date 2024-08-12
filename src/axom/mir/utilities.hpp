// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UTILITIES_HPP_
#define AXOM_MIR_UTILITIES_HPP_

#include "axom/core.hpp"
//#include "axom/core/Array.hpp"
//#include "axom/core/ArrayView.hpp"
//#include "axom/core/memory_management.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#include <RAJA/RAJA.hpp>

#include <cstdint>

// NOTE: Longer term, we should hide more RAJA functionality behind axom wrappers
//       so we can write serial versions for when RAJA is not enabled.
namespace axom
{
template <typename ExecSpace, typename ContiguousMemoryContainer>
struct scans
{
  inline void exclusive_scan(const ContiguousMemoryContainer &input,
                             ContiguousMemoryContainer &output)
  {
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    assert(input.size() == output.size());
    RAJA::exclusive_scan<loop_policy>(
      RAJA::make_span(input.data(), input.size()),
      RAJA::make_span(output.data(), output.size()));
  }
};

template <typename ExecSpace, typename ContiguousMemoryContainer>
inline void exclusive_scan(const ContiguousMemoryContainer &input,
                           ContiguousMemoryContainer &output)
{
  scans<ExecSpace, ContiguousMemoryContainer> s;
  s.exclusive_scan(input, output);
}

}  // end namespace axom

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
struct accumulation_traits
{
  using value_type = float;
};

template <>
struct accumulation_traits<double>
{
  using value_type = double;
};

template <>
struct accumulation_traits<long>
{
  using value_type = double;
};

template <>
struct accumulation_traits<unsigned long>
{
  using value_type = double;
};

//------------------------------------------------------------------------------
/**
 * \brief Use binary search to find the index of the \a value in the supplied
 *        sorted view.
 *
 * \param[in] value The search value.
 * \param[in] view The view that contains the sorted search data values.
 *
 * \return The index where value was located in view or -1 if not found.
 */
template <typename T>
AXOM_HOST_DEVICE std::int32_t bsearch(T value, const axom::ArrayView<T> &view)
{
  std::int32_t index = -1;
  std::int32_t left = 0;
  std::int32_t right = view.size() - 1;
  while(left <= right)
  {
    std::int32_t m = (left + right) / 2;
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
/**
 * \brief Hash a stream of bytes into a uint64_t hash value.
 *
 * \param[in] data The bytes to be hashed.
 * \param[in] length The number of bytes to be hashed.
 *
 * \return A uint64_t hash for the byte stream.
 *
 * \note The byte stream is hashed using a Jenkins-hash algorithm forwards and
 *       backwards and the two results are merged into a uint64_t. The length is
 *       also part of the hash to guard against a lot of repeated values in the
 *       byte stream hashing to the same thing.
 *
 * \note We make this function inline since it is not a template and we want to
 *       use it in both host and device code.
 */
AXOM_HOST_DEVICE
inline std::uint64_t hash_bytes(const std::uint8_t *data, std::uint32_t length)
{
  std::uint32_t hash = 0;

  // Build the length into the hash.
  const auto ldata = reinterpret_cast<const std::uint8_t *>(&length);
  for(int e = 0; e < 4; e++)
  {
    hash += ldata[e];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  std::uint32_t hashr = hash;
  for(std::uint32_t i = 0; i < length; i++)
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

  return (static_cast<std::uint64_t>(hash) << 32) | hashr;
}

//------------------------------------------------------------------------------
/**
 * \brief A very basic sort for small arrays.
 * \param[inout] v The values to be sorted in place.
 * \param[in] n The number of values to be sorted.
 */
template <typename ValueType, typename IndexType>
AXOM_HOST_DEVICE void sort_values(ValueType *v, IndexType n)
{
  for(IndexType i = 0; i < n - 1; i++)
  {
    const IndexType m = n - i - 1;
    for(IndexType j = 0; j < m; j++)
    {
      if(v[j] > v[j + 1])
      {
        axom::utilities::swap(v[j], v[j + 1]);
      }
    }
  }
}

//------------------------------------------------------------------------------
/**
 * \brief Make a hashed "name" for one id.
 *
 * \param[in] id The id we're hashing.
 * \return A hashed name for the id.
 */
template <typename ValueType>
AXOM_HOST_DEVICE std::uint64_t make_name_1(ValueType id)
{
  return hash_bytes(reinterpret_cast<std::uint8_t *>(&id), sizeof(ValueType));
};

//------------------------------------------------------------------------------
/**
 * \brief Make a hashed "name" for two ids.
 *
 * \param[in] id0 The first id we're hashing.
 * \param[in] id1 The second id we're hashing.
 * \return A hashed name for the ids.
 */
template <typename ValueType>
AXOM_HOST_DEVICE std::uint64_t make_name_2(ValueType id0, ValueType id1)
{
  ValueType data[2] = {id0, id1};
  if(id1 < id0)
  {
    data[0] = id1;
    data[1] = id0;
  }
  return hash_bytes(reinterpret_cast<std::uint8_t *>(data),
                    2 * sizeof(ValueType));
};

//------------------------------------------------------------------------------
/**
 * \brief Make a hashed "name" for multiple ids.
 *
 * \tparam MaxValues The largest expected number of values.
 *
 * \param[in] values The ids that are being hashed.
 * \param[in] start The starting index for ids.
 * \param[in] n The number if ids being hashed.
 *
 * \return A hashed name for the ids.
 */
template <typename ValueType, std::uint32_t MaxValues = 14>
AXOM_HOST_DEVICE std::uint64_t make_name_n(const ValueType *values,
                                           std::uint32_t n)
{
  assert(n <= MaxValues);
  if(n == 2) return make_name_2(values[0], values[1]);

  // Make sure the values are sorted before hashing.
  ValueType sorted[MaxValues];
  for(std::uint32_t i = 0; i < n; i++)
  {
    sorted[i] = values[i];
  }
  sort_values(sorted, n);

  return hash_bytes(reinterpret_cast<std::uint8_t *>(sorted),
                    n * sizeof(ValueType));
}

//------------------------------------------------------------------------------
/**
 * \brief This class implements a naming policy that uses some hashing functions
 *        to produce a "name" for an array of ids.
 */
class HashNaming
{
public:
  using KeyType = std::uint64_t;

  /**
   * \brief Make a name from an array of ids.
   *
   * \param ids The array of ids.
   * \param numIds The number of ids in the array.
   *
   * \return The name that describes the array of ids.
   *
   * \note Different make_name_* functions are used because we can skip most
   *       sorting for 1,2 element arrays.
   */
  AXOM_HOST_DEVICE
  inline static KeyType makeName(const IndexType *ids, IndexType numIds)
  {
    KeyType name {};
    if(numIds == 1)
    {
      name = axom::mir::utilities::make_name_1(ids[0]);
    }
    else if(numIds == 2)
    {
      name = axom::mir::utilities::make_name_2(ids[0], ids[1]);
    }
    else
    {
      name = axom::mir::utilities::make_name_n(ids, numIds);
    }
    return name;
  }
};

//------------------------------------------------------------------------------
/**
 * \brief This function makes a unique array of values from an input list of keys.
 *
 * \tparam ExecSpace The execution space.
 * \tparam KeyType   The data type for the keys.
 *
 * \param[in] keys_orig_view The input view that contains the input keys to be made unique.
 * \param[out] skeys     A sorted unique array of keys produced from keys_orig_view.
 * \param[out] sindices  An array of indices that indicate where in the original view the keys came from.
 *
 */
template <typename ExecSpace, typename KeyType>
void unique(const axom::ArrayView<KeyType> &keys_orig_view,
            axom::Array<KeyType> &skeys,
            axom::Array<axom::IndexType> &sindices)
{
  using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Make a copy of the keys and make original indices.
  const auto n = keys_orig_view.size();
  axom::Array<KeyType> keys(n, n, allocatorID);
  axom::Array<axom::IndexType> indices(n, n, allocatorID);
  auto keys_view = keys.view();
  auto indices_view = indices.view();
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) {
      keys_view[i] = keys_orig_view[i];
      indices_view[i] = i;
    });

  // Sort the keys, indices in place.
  RAJA::sort_pairs<loop_policy>(RAJA::make_span(keys_view.data(), n),
                                RAJA::make_span(indices_view.data(), n));

  // Make a mask array for where differences occur.
  axom::Array<axom::IndexType> mask(n, n, allocatorID);
  auto mask_view = mask.view();
  RAJA::ReduceSum<reduce_policy, axom::IndexType> mask_sum(0);
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) {
      const axom::IndexType m =
        (i >= 1) ? ((keys_view[i] != keys_view[i - 1]) ? 1 : 0) : 1;
      mask_view[i] = m;
      mask_sum += m;
    });

  // Do a scan on the mask array to build an offset array.
  axom::Array<axom::IndexType> offsets(n, n, allocatorID);
  auto offsets_view = offsets.view();
  RAJA::exclusive_scan<loop_policy>(RAJA::make_span(mask_view.data(), n),
                                    RAJA::make_span(offsets_view.data(), n),
                                    RAJA::operators::plus<axom::IndexType> {});

  // Allocate the output arrays.
  const axom::IndexType newsize = mask_sum.get();
  skeys = axom::Array<KeyType>(newsize, newsize, allocatorID);
  sindices = axom::Array<axom::IndexType>(newsize, newsize, allocatorID);

  // Iterate over the mask/offsets to store values at the right
  // offset in the new array.
  auto skeys_view = skeys.view();
  auto sindices_view = sindices.view();
  axom::for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(axom::IndexType i) {
      if(mask_view[i])
      {
        skeys_view[offsets_view[i]] = keys_view[i];
        sindices_view[offsets_view[i]] = indices_view[i];
      }
    });
}

}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
