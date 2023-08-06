// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_Detail_FlatTable_Hpp
#define Axom_Core_Detail_FlatTable_Hpp

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/BitUtilities.hpp"

namespace axom
{
namespace detail
{
namespace flat_map
{
struct QuadraticProbing
{
  /*!
     * \brief Returns the next offset to jump given an iteration count.
     *
     *  Suppose we start with some probe point H + i^2, then the next probe
     *  point is H + (i*1)^2 = H + i^2 + 2*i + 1.
     *
     *  We return the offset from H(i) to H(i+1), which is 2*i+1.
     */
  int getNext(int iter) const { return 2 * iter + 1; }
};

// Boost::unordered_flat_map uses a 128-bit chunk of metadata for each
// group of 15 buckets.
// This is split up into an "overflow bit", and 15 bytes representing the
// state of each bucket.
//
// Each bucket byte has the value:
// - 0 if the bucket is "deleted"
// - 1 to signal the end of the bucket array.
// - if the bucket has an element, a reduced hash in the range [2, 255].
//
// The overflow bit acts as a Bloom filter to determine if probing should
// terminate.
struct GroupBucket
{
  constexpr static std::uint8_t Empty = 0;
  constexpr static std::uint8_t Sentinel = 1;
  constexpr static int InvalidSlot = -1;

  GroupBucket() : data {0ULL, 0ULL} { }

  int getEmptyBucket() const
  {
    for(int i = 0; i < 15; i++)
    {
      if(metadata.buckets[i] == GroupBucket::Empty)
      {
        return i;
      }
    }
    // Bucket not found.
    return InvalidSlot;
  }

  int getHashBucket(std::uint8_t hash) const
  {
    std::uint8_t reducedHash = reduceHash(hash);
    for(int i = 0; i < 15; i++)
    {
      if(metadata.buckets[i] == reducedHash)
      {
        return i;
      }
    }
    return InvalidSlot;
  }

  void setBucket(int index, std::uint8_t hash)
  {
    metadata.buckets[index] = reduceHash(hash);
  }

  void setOverflow(std::uint8_t hash)
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    metadata.ofw |= hashOfwBit;
  }

  bool getMaybeOverflowed(std::uint8_t hash) const
  {
    std::uint8_t hashOfwBit = 1 << (hash % 8);
    return (metadata.ofw & hashOfwBit);
  }

  bool hasSentinel() const { return metadata.buckets[14] == Sentinel; }

  void setSentinel() { metadata.buckets[14] = Sentinel; }

  // We need to map hashes in the range [0, 255] to [2, 255], since 0 and 1
  // are taken by the "empty" and "sentinel" values respectively.
  static std::uint8_t reduceHash(std::uint8_t hash)
  {
    return (hash < 2) ? (hash + 8) : hash;
  }

  union alignas(16)
  {
    struct
    {
      std::uint8_t ofw;
      std::uint8_t buckets[15];
    } metadata;
    std::uint64_t data[2];
    static_assert(
      sizeof(metadata) == sizeof(data),
      "FlatMap::SwissTable::Bucket: sizeof(data_bytes) != sizeof(data)");
  };
};

static_assert(sizeof(GroupBucket) == 16,
              "FlatMap::SwissTable::Bucket: size != 16 bytes");
static_assert(std::alignment_of<GroupBucket>::value == 16,
              "FlatMap::SwissTable::Bucket: alignment != 16 bytes");
static_assert(std::is_standard_layout<GroupBucket>::value,
              "FlatMap::SwissTable::Bucket: not standard layout");

struct FlatMapIndex
{
  int group_index;
  int bucket_index;
};

template <typename HashType, typename ProbePolicy = QuadraticProbing>
struct SequentialLookupPolicy
{
  constexpr static int NO_MATCH = -1;

  /*!
     * \brief Inserts a hash into the first empty bucket in an array of groups
     *  for an open-addressing hash map.
     *
     * \param [in] ngroups_pow_2 the number of groups, expressed as a power of 2
     * \param [in] groups the array of metadata for the groups in the hash map
     * \param [in] hash the hash to insert
     *
     * \return a FlatMapIndex with group and bucket ID, or NO_MATCH if no empty
     *  buckets were found.
     */
  FlatMapIndex insert(int ngroups_pow_2,
                      ArrayView<GroupBucket> groups,
                      HashType hash) const
  {
    // We use the k MSBs of the hash as the initial group probe point,
    // where ngroups = 2^k.
    int group_divisor = 1 << ((CHAR_BIT * sizeof(HashType)) - ngroups_pow_2);
    int curr_group = hash / group_divisor;

    std::uint8_t hash_8 {hash};
    int iteration = 0;
    while(curr_group != NO_MATCH)
    {
      int empty_bucket = groups[curr_group].getEmptyBucket();
      if(empty_bucket != GroupBucket::InvalidSlot)
      {
        groups[curr_group].setBucket(empty_bucket, hash_8);
        return FlatMapIndex {curr_group, empty_bucket};
      }
      else if(!groups[curr_group].hasSentinel())
      {
        // Set the overflow bit and continue probing.
        groups[curr_group].setOverflow(hash_8);
        curr_group =
          (curr_group + ProbePolicy {}.getNext(iteration)) % groups.size();
        iteration++;
      }
      else
      {
        // Stop probing.
        curr_group = NO_MATCH;
      }
    }
    return FlatMapIndex {NO_MATCH, NO_MATCH};
  }

  /*!
     * \brief Finds the next potential bucket index for a given hash in a group
     *  array for an open-addressing hash map.
     *
     * \param [in] ngroups_pow_2 the number of groups, expressed as a power of 2
     * \param [in] groups the array of metadata for the groups in the hash map
     * \param [in] hash the hash to insert
     *
     * \return a FlatMapIndex with group and bucket ID, or NO_MATCH if no empty
     *  buckets were found.
     */
  FlatMapIndex find(int ngroups_pow_2,
                    ArrayView<const GroupBucket> groups,
                    HashType hash,
                    int last = NO_MATCH) const
  {
    int curr_group = last;
    int iteration = 0;
    if(curr_group == NO_MATCH)
    {
      // We use the k MSBs of the hash as the initial group probe point,
      // where ngroups = 2^k.
      int group_divisor = 1 << ((CHAR_BIT * sizeof(HashType)) - ngroups_pow_2);
      curr_group = hash / group_divisor;
    }
    std::uint8_t hash_8 {hash};
    while(curr_group != NO_MATCH)
    {
      int bucket = groups[curr_group].getHashBucket(hash_8);
      if(bucket != GroupBucket::InvalidSlot)
      {
        return FlatMapIndex {curr_group, bucket};
      }
      if(!groups[curr_group].getMaybeOverflowed(hash_8) ||
         groups[curr_group].hasSentinel())
      {
        // Stop probing if the "overflow" bit is not set or if the sentinel
        // is set for a bucket.
        curr_group = NO_MATCH;
      }
      else
      {
        curr_group =
          (curr_group + ProbePolicy {}.getNext(iteration)) % groups.size();
        iteration++;
      }
    }
    // Unable to find a matching entry for a hash.
    return FlatMapIndex {NO_MATCH, NO_MATCH};
  }
};

}  // namespace flat_map
}  // namespace detail
}  // namespace axom

#endif  // Axom_Core_Detail_FlatTable_Hpp
