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
#if 0
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
#else
template <typename IndexT, int MAXIDS = 14>
class HashNaming
{
public:
  using KeyType = std::uint64_t;
  using IndexType = IndexT;

  // The top 2 bits are reserved for the key type.
  constexpr static KeyType KeyIDSingle = 0;
  constexpr static KeyType KeyIDPair = KeyType(1) << 62;
  constexpr static KeyType KeyIDPack = KeyType(2) << 62;
  constexpr static KeyType KeyIDHash = KeyType(3) << 62;

  // The rest of the bits can be divided in various ways.
  constexpr static KeyType KeyMask = KeyType(3) << 62;
  constexpr static KeyType PayloadMask = ~KeyMask;
  constexpr static KeyType Max15Bit = (KeyType(1) << 15) - 1;
  constexpr static KeyType Max16Bit = (KeyType(1) << 16) - 1;
  constexpr static KeyType Max20Bit = (KeyType(1) << 20) - 1;
  constexpr static KeyType Max31Bit = (KeyType(1) << 31) - 1;
  constexpr static KeyType Max32Bit = (KeyType(1) << 32) - 1;

  class View
  {
  public:
    using KeyType = HashNaming::KeyType;

    AXOM_HOST_DEVICE
    KeyType makeName(const IndexType* p, int n) const
    {
      KeyType retval{};
      if (n == 1)
        retval = make_name_1(p[0]);
      else if (n == 2)
        retval = make_name_2(p[0], p[1]);
      else
        retval = make_name_n(p, n);
      return retval;
    }

    AXOM_HOST_DEVICE
    void setMaxId(IndexType m)
    {
      m_maxId = static_cast<KeyType>(m);
    }

  private:
    AXOM_HOST_DEVICE
    inline KeyType make_name_1(IndexType p0) const
    {
      assert(p0 < PayloadMask);
      // Store p0 in the key as a 62-bit integer
      KeyType k0 = (static_cast<KeyType>(p0) & PayloadMask);
      return KeyIDSingle | k0;
    }

    AXOM_HOST_DEVICE
    inline KeyType make_name_2(IndexType p0, IndexType p1) const
    {
      assert(p0 <= Max31Bit && p1 <= Max31Bit);
      // Store p0 and p1 both in the 64-bit key as 31-bit integers
      KeyType k0 = (static_cast<KeyType>(std::min(p0, p1)) & Max31Bit);
      KeyType k1 = (static_cast<KeyType>(std::max(p0, p1)) & Max31Bit);
      return KeyIDPair | (k0 << 31) | k1;  
    }
 
    AXOM_HOST_DEVICE
    KeyType make_name_n(const IndexType *p, int n) const
    {
      KeyType retval {};
      if(n == 3 && m_maxId <= Max20Bit)
      {
        // We can pack 3 values into the id lossless
        IndexType sorted[3];
        sorted[0] = p[0];
        sorted[1] = p[1];
        sorted[2] = p[2];
        axom::utilities::Sorting<IndexType, 3>::sort(sorted);

        KeyType k0 = static_cast<KeyType>(sorted[0]) & Max20Bit;
        KeyType k1 = static_cast<KeyType>(sorted[1]) & Max20Bit;
        KeyType k2 = static_cast<KeyType>(sorted[2]) & Max20Bit;
        constexpr KeyType len = KeyType(3 - 1) << 60;
        retval = KeyIDPack | len | (k0 << 40) | (k1 << 20) | k2;
      }
      else if (n == 4 && m_maxId <= Max15Bit)
      {
          // We can pack 4 values into the id lossless
          IndexType sorted[4];
          sorted[0] = p[0];
          sorted[1] = p[1];
          sorted[2] = p[2];
          sorted[3] = p[3];
          axom::utilities::Sorting<IndexType, 4>::sort(sorted);

          KeyType k0 = static_cast<KeyType>(sorted[0]) & Max15Bit;
          KeyType k1 = static_cast<KeyType>(sorted[1]) & Max15Bit;
          KeyType k2 = static_cast<KeyType>(sorted[2]) & Max15Bit;
          KeyType k3 = static_cast<KeyType>(sorted[3]) & Max15Bit;
          constexpr KeyType len = KeyType(4 - 1) << 60;
          retval = KeyIDPack | len | (k0 << 45) | (k1 << 30) | (k2 << 15) | k3;
      }
      else if(m_maxId < Max16Bit)
      {
         // Narrow to 16-bit, sort
         std::uint16_t sorted[MAXIDS];
         for(int i = 0; i < n; i++)
           sorted[i] = static_cast<std::uint16_t>(p[i]);
         axom::utilities::Sorting<std::uint16_t, MAXIDS>::sort(sorted, n);

         // Make a hash from the narrowed ids
         void *ptr = static_cast<void *>(sorted);
         KeyType k0 = axom::mir::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n << 1);
         retval = KeyIDHash | (k0 & PayloadMask);
      }
      else if (m_maxId < Max32Bit)
      {
          // Narrow to 32-bit, sort
          std::uint32_t sorted[MAXIDS];
          for (int i = 0; i < n; i++)
            sorted[i] = static_cast<std::uint32_t>(p[i]);
          axom::utilities::Sorting<std::uint32_t, MAXIDS>::sort(sorted, n);

          // Make a hash from the narrowed ids
          void *ptr = static_cast<void *>(sorted);
          KeyType k0 = axom::mir::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n << 2);
          retval = KeyIDHash | (k0 & PayloadMask);
      }
      else
      {
          IndexType sorted[MAXIDS];
          for (int i = 0; i < n; i++)
            sorted[i] = p[i];
          axom::utilities::Sorting<IndexType, MAXIDS>::sort(sorted, n);

          // Make a hash from the ids
          void *ptr = static_cast<void *>(sorted);
          KeyType k0 = axom::mir::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n * sizeof(IndexType));
          retval = KeyIDHash | (k0 & PayloadMask);
      }
      return retval;
    }

    KeyType m_maxId {};
  };

  // Host-callable methods

  KeyType makeName(const IndexType* p, int n) const
  {
    return m_view.makeName(p, n);
  }

  void setMaxId(IndexType n)
  {
    m_view.setMaxId(n);
  }

  View view()
  {
    return m_view;
  }

  static std::string toString(KeyType key)
  {
      std::stringstream ss;
      auto kt = key & KeyMask;
      if(kt == KeyIDSingle)
      {
        auto id = key & PayloadMask;
        ss << "single("<< std::hex << id << ")";
      }
      else if(kt == KeyIDPair)
      {
        auto payload = key & PayloadMask;
        auto p0 = (payload >> 31) & Max31Bit;
        auto p1 = payload & Max31Bit;
        ss << "pair(" << std::hex << p0 << ", " << p1 << ")";
      }
      else if(kt == KeyIDHash)
      {
        ss << "hash(" << std::hex << key << ")";
      }
      else if(kt == KeyIDPack)
      {
        auto npts = ((key >> 60) & 3) + 1;
        if(npts == 3)
        {
          auto p0 = (key >> 40) & Max20Bit;
          auto p1 = (key >> 20) & Max20Bit;         
          auto p2 = (key) & Max20Bit;
          ss << "pack(" << std::hex << p0 << ", " << p1 << ", " << p2 << ")";
        }
        else if(npts == 4)
        {
          auto p0 = (key >> 45) & Max15Bit;
          auto p1 = (key >> 30) & Max15Bit;
          auto p2 = (key >> 15) & Max15Bit;         
          auto p3 = (key) & Max15Bit;
          ss << "pack(" << std::hex << p0 << ", " << p1 << ", " << p2 << ", " << p3 << ")";
        }
      }
      return ss.str();
  }

  View m_view {};
};

#endif

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
