// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_HPP
#define Axom_Core_FlatMap_HPP

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/detail/FlatTable.hpp"

namespace axom
{
template <typename KeyType, typename ValueType, typename Hash = std::hash<KeyType>>
class FlatMap
  : detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>
{
private:
  using LookupPolicy =
    detail::flat_map::SequentialLookupPolicy<typename Hash::result_type>;
  using LookupPolicy::NO_MATCH;

  constexpr static int BucketsPerGroup = detail::flat_map::GroupBucket::Size;

  template <bool Const>
  class IteratorImpl;

  using KeyValuePair = std::pair<const KeyType, ValueType>;

  template <bool Const>
  friend class IteratorImpl;

  using MixedHash = detail::flat_map::HashMixer64<KeyType, Hash>;

public:
  using size_type = IndexType;
  using value_type = KeyValuePair;
  using iterator = IteratorImpl<false>;
  using const_iterator = IteratorImpl<true>;

  // Constructors
  FlatMap() : FlatMap(MIN_NUM_BUCKETS) { }

  explicit FlatMap(IndexType bucket_count);

  template <typename InputIt>
  FlatMap(InputIt first, InputIt last, IndexType bucket_count = -1)
    : FlatMap(std::distance(first, last), first, last, bucket_count)
  { }

  explicit FlatMap(std::initializer_list<value_type> init,
                   IndexType bucket_count = -1)
    : FlatMap(init.begin(), init.end(), bucket_count)
  { }

  void swap(FlatMap& other)
  {
    axom::utilities::swap(m_numGroups2, other.m_numGroups2);
    axom::utilities::swap(m_size, other.m_size);
    axom::utilities::swap(m_metadata, other.m_metadata);
    axom::utilities::swap(m_buckets, other.m_buckets);
    axom::utilities::swap(m_loadCount, other.m_loadCount);
  }

  // Iterators
  iterator begin()
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return iterator(this, firstBucketIndex);
  }
  const_iterator begin() const
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return const_iterator(this, firstBucketIndex);
  }
  const_iterator cbegin() const
  {
    IndexType firstBucketIndex = this->nextValidIndex(m_metadata, NO_MATCH);
    return const_iterator(this, firstBucketIndex);
  }

  iterator end() { return iterator(this, m_buckets.size()); }
  const_iterator end() const { return const_iterator(this, m_buckets.size()); }
  const_iterator cend() const { return const_iterator(this, m_buckets.size()); }

  // Capacity
  bool empty() const { return m_size == 0; }
  IndexType size() const { return m_size; }

  // Lookup
  iterator find(const KeyType& key);
  const_iterator find(const KeyType& key) const;

  ValueType& at(const KeyType& key) { return this->find(key)->second; }
  const ValueType& at(const KeyType& key) const
  {
    return this->find(key)->second;
  }

  ValueType& operator[](const KeyType& key) { return at(key); }
  const ValueType& operator[](const KeyType& key) const { return at(key); }

  IndexType count(const KeyType& key) const { return (find(key) != end()); }
  bool contains(const KeyType& key) const { return (find(key) != end()); }

  // Modifiers
  void clear();
  std::pair<iterator, bool> insert(const value_type& value)
  {
    return emplaceImpl(false, value.first, value.second);
  }
  std::pair<iterator, bool> insert(value_type&& value)
  {
    return emplaceImpl(false, std::move(value.first), std::move(value.second));
  }
  template <typename InputPair>
  std::pair<iterator, bool> insert(InputPair&& pair)
  {
    return emplace(std::forward<InputPair>(pair));
  }
  template <typename... InputArgs>
  std::pair<iterator, bool> emplace(InputArgs&&... pair)
  {
    return emplaceImpl(false, std::forward<InputArgs>(pair)...);
  }

  template <typename InputIt>
  void insert(InputIt first, InputIt last);

  template <typename... Args>
  std::pair<iterator, bool> insert_or_assign(const KeyType& key, Args&&... args)
  {
    return emplaceImpl(true, KeyType {key}, std::forward<Args>(args)...);
  }
  template <typename... Args>
  std::pair<iterator, bool> insert_or_assign(KeyType&& key, Args&&... args)
  {
    return emplaceImpl(true, std::move(key), std::forward<Args>(args)...);
  }

  template <typename... Args>
  std::pair<iterator, bool> try_emplace(const KeyType& key, Args&&... args)
  {
    return emplaceImpl(false, KeyType {key}, std::forward<Args>(args)...);
  }
  template <typename... Args>
  std::pair<iterator, bool> try_emplace(KeyType&& key, Args&&... args)
  {
    return emplaceImpl(false, std::move(key), std::forward<Args>(args)...);
  }

  iterator erase(iterator pos) { erase(const_iterator {pos}); }
  iterator erase(const_iterator pos);
  IndexType erase(const KeyType& key)
  {
    const_iterator it = find(key);
    if(it != end())
    {
      erase(it);
      return 1;
    }
    return 0;
  }

  // Hashing
  IndexType bucket_count() const { return m_buckets.size(); }
  double load_factor() const
  {
    return ((double)m_loadCount) / m_buckets.size();
  }
  double max_load_factor() const { return MAX_LOAD_FACTOR; }
  void rehash(IndexType count)
  {
    FlatMap rehashed(m_size, cbegin(), cend(), count);
    this->swap(rehashed);
  }
  void reserve(IndexType count) { rehash(std::ceil(count / MAX_LOAD_FACTOR)); }

private:
  template <typename InputIt>
  FlatMap(IndexType num_elems, InputIt first, InputIt last, IndexType bucket_count);

  template <typename UKeyType, typename... Args>
  std::pair<iterator, bool> emplaceImpl(bool assign_on_existence,
                                        UKeyType&& key,
                                        Args&&... args);

  constexpr static IndexType MIN_NUM_BUCKETS {29};

  IndexType m_numGroups2;  // Number of groups of 15 buckets, expressed as a power of 2
  IndexType m_size;
  axom::Array<detail::flat_map::GroupBucket> m_metadata;
  axom::Array<KeyValuePair> m_buckets;

  // Boost flat_unordered_map uses a fixed load factor.
  constexpr static double MAX_LOAD_FACTOR = 0.875;
  std::uint64_t m_loadCount;
};

template <typename KeyType, typename ValueType, typename Hash>
template <bool Const>
class FlatMap<KeyType, ValueType, Hash>::IteratorImpl
{
private:
  using MapType = FlatMap<KeyType, ValueType, Hash>;

  template <bool OtherConst>
  friend class IteratorImpl;

  friend class FlatMap<KeyType, ValueType, Hash>;

public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = MapType::value_type;
  using difference_type = IndexType;

  using DataType = std::conditional_t<Const, const value_type, value_type>;
  using MapConstType = std::conditional_t<Const, const MapType, MapType>;
  using pointer = DataType*;
  using reference = DataType&;

public:
  IteratorImpl(MapConstType* map, IndexType internalIdx)
    : m_map(map)
    , m_internalIdx(internalIdx)
  {
    assert(m_internalIdx >= 0 && m_internalIdx <= m_map->m_buckets.size());
  }

  template <bool UConst = Const, typename Enable = std::enable_if_t<UConst>>
  IteratorImpl(IteratorImpl<!Const> non_const_iter)
    : m_map(non_const_iter.m_map)
    , m_internalIdx(non_const_iter.m_internalIdx)
  { }

  friend bool operator==(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return lhs.m_internalIdx == rhs.m_internalIdx;
  }

  friend bool operator!=(const IteratorImpl& lhs, const IteratorImpl& rhs)
  {
    return lhs.m_internalIdx != rhs.m_internalIdx;
  }

  IteratorImpl& operator++()
  {
    m_internalIdx = m_map->nextValidIndex(m_map->m_metadata, m_internalIdx);
    return *this;
  }

  IteratorImpl operator++(int)
  {
    IteratorImpl next = *this;
    ++(*this);
    return next;
  }

  reference operator*() const { return m_map->m_buckets[m_internalIdx]; }

  pointer operator->() const { return &(m_map->m_buckets[m_internalIdx]); }

private:
  MapConstType* m_map;
  IndexType m_internalIdx;
};

template <typename KeyType, typename ValueType, typename Hash>
FlatMap<KeyType, ValueType, Hash>::FlatMap(IndexType bucket_count)
  : m_size(0)
  , m_loadCount(0)
{
  int minBuckets = MIN_NUM_BUCKETS;
  bucket_count = std::max(minBuckets, bucket_count);
  // Get the smallest power-of-two number of groups satisfying:
  // N * GroupSize - 1 >= minBuckets
  // TODO: we should add a leadingZeros overload for 64-bit integers
  {
    std::int32_t numGroups =
      std::ceil((bucket_count + 1) / (double)BucketsPerGroup);
    m_numGroups2 = 31 - (axom::utilities::leadingZeros(numGroups));
  }

  IndexType numGroupsRounded = 1 << m_numGroups2;
  IndexType numBuckets = numGroupsRounded * BucketsPerGroup - 1;
  m_metadata.resize(numGroupsRounded);
  m_metadata[numGroupsRounded - 1].setSentinel();
  m_buckets.resize(numBuckets);
}

template <typename KeyType, typename ValueType, typename Hash>
template <typename InputIt>
FlatMap<KeyType, ValueType, Hash>::FlatMap(IndexType num_elems,
                                           InputIt first,
                                           InputIt last,
                                           IndexType bucket_count)
  : FlatMap(std::max(num_elems, bucket_count))
{
  insert(first, last);
}

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::find(const KeyType& key) -> iterator
{
  auto hash = MixedHash {}(key);
  iterator found_iter = end();
  this->probeIndex(m_numGroups2,
                   m_metadata,
                   hash,
                   [&](IndexType bucket_index) -> bool {
                     if(this->m_buckets[bucket_index].first == key)
                     {
                       found_iter = iterator(this, bucket_index);
                       // Stop tracking.
                       return false;
                     }
                     return true;
                   });
  return found_iter;
}

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::find(const KeyType& key) const
  -> const_iterator
{
  auto hash = MixedHash {}(key);
  const_iterator found_iter = end();
  this->probeIndex(m_numGroups2,
                   m_metadata,
                   hash,
                   [&](IndexType bucket_index) -> bool {
                     if(this->m_buckets[bucket_index].first == key)
                     {
                       found_iter = const_iterator(this, bucket_index);
                       // Stop tracking.
                       return false;
                     }
                     return true;
                   });
  return found_iter;
}

template <typename KeyType, typename ValueType, typename Hash>
template <typename InputIt>
void FlatMap<KeyType, ValueType, Hash>::insert(InputIt first, InputIt last)
{
  while(first != last)
  {
    insert(*first);
    first++;
  }
}

template <typename KeyType, typename ValueType, typename Hash>
template <typename UKeyType, typename... InputArgs>
auto FlatMap<KeyType, ValueType, Hash>::emplaceImpl(bool assign_on_existence,
                                                    UKeyType&& key,
                                                    InputArgs&&... args)
  -> std::pair<iterator, bool>
{
  static_assert(std::is_convertible<UKeyType, KeyType>::value,
                "UKeyType -> KeyType not convertible");
  auto hash = MixedHash {}(key);
  // Resize to double the number of bucket groups if insertion would put us
  // above the maximum load factor.
  if(((m_loadCount + 1) / (double)m_buckets.size()) >= MAX_LOAD_FACTOR)
  {
    IndexType newNumGroups = m_metadata.size() * 2;
    rehash(newNumGroups * BucketsPerGroup - 1);
  }

  bool keyExistsAlready = false;
  IndexType foundBucketIndex = NO_MATCH;
  auto FindExistingElem = [&, this](IndexType bucket_index) -> bool {
    if(this->m_buckets[bucket_index].first == key)
    {
      keyExistsAlready = true;
      foundBucketIndex = bucket_index;
      // Exit out of probing, we can't insert if the key already exists.
      return true;
    }
    return false;
  };

  IndexType newBucket =
    this->probeEmptyIndex(m_numGroups2, m_metadata, hash, FindExistingElem);
  if(!keyExistsAlready)
  {
    foundBucketIndex = newBucket;
    // Add a hash to the corresponding bucket slot.
    this->setBucketHash(m_metadata, newBucket, hash);
    m_size++;
    m_loadCount++;
  }
  iterator keyIterator = iterator(this, foundBucketIndex);
  if(!keyExistsAlready)
  {
    new(&m_buckets[foundBucketIndex])
      KeyValuePair(std::move(key), std::forward<InputArgs>(args)...);
  }
  else if(keyExistsAlready && assign_on_existence)
  {
    m_buckets[foundBucketIndex].second =
      ValueType(std::forward<InputArgs>(args)...);
  }
  return {keyIterator, !keyExistsAlready};
}

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::erase(const_iterator pos) -> iterator
{
  assert(pos != end());
  auto hash = MixedHash {}(pos->first);

  bool midSequence = this->clearBucket(m_metadata, pos.m_internalIdx, hash);
  pos->~KeyValuePair();
  if(!midSequence)
  {
    m_loadCount--;
  }
  return ++iterator(this, pos.m_internalIdx);
}

}  // namespace axom

#endif  // Axom_Core_FlatMap_HPP
