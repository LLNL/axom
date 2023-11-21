// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_HPP
#define Axom_Core_FlatMap_HPP

#include <tuple>
#include <utility>
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

  /*!
   * \brief Constructs a FlatMap with no elements.
   */
  FlatMap() : FlatMap(MIN_NUM_BUCKETS) { }

  /*!
   * \brief Constructs a FlatMap with at least a given number of buckets.
   *
   * \param [in] bucket_count the minimum number of buckets to allocate
   */
  explicit FlatMap(IndexType bucket_count);

  /*!
   * \brief Constructs a FlatMap with a range of elements.
   *
   * \param [in] first iterator pointing to the beginning of the range
   * \param [in] last iterator pointing to the end of the range
   * \param [in] bucket_count minimum number of buckets to allocate (optional)
   */
  template <typename InputIt>
  FlatMap(InputIt first, InputIt last, IndexType bucket_count = -1)
    : FlatMap(std::distance(first, last), first, last, bucket_count)
  { }

  /*!
   * \brief Constructs a FlatMap with a range of elements.
   *
   * \param [in] init a list of pairs to initialize the map with
   * \param [in] bucket_count minimum number of buckets to allocate (optional)
   */
  explicit FlatMap(std::initializer_list<value_type> init,
                   IndexType bucket_count = -1)
    : FlatMap(init.begin(), init.end(), bucket_count)
  { }

  /*!
   * \brief Move constructor for a FlatMap instance.
   *
   * \param other the FlatMap to move data from
   */
  FlatMap(FlatMap&& other) : FlatMap() { swap(other); }

  /*!
   * \brief Move assignment operator for a FlatMap instance.
   *
   * \param other the FlatMap to move data from
   */
  FlatMap& operator=(FlatMap&& other) { swap(other); }

  /*!
   * \brief Copy constructor for a FlatMap instance.
   *
   * \param other the FlatMap to copy data from
   */
  FlatMap(const FlatMap& other)
    : m_numGroups2(other.m_numGroups2)
    , m_size(other.m_size)
    , m_metadata(other.m_metadata)
    , m_buckets(other.m_buckets.size())
    , m_loadCount(other.m_loadCount)
  {
    // Copy all elements.
    IndexType index = this->nextValidIndex(m_metadata, NO_MATCH);
    while(index < bucket_count())
    {
      new(&m_buckets[index].data) KeyValuePair(other.m_buckets[index].get());
      index = this->nextValidIndex(m_metadata, index);
    }
  }

  /*!
   * \brief Copy assignment operator for a FlatMap instance.
   *
   * \param other the FlatMap to copy data from
   */
  FlatMap& operator=(const FlatMap& other)
  {
    if(*this != other)
    {
      FlatMap new_map(other);
      swap(new_map);
    }
  }

  /// \brief Destructor for a FlatMap instance.
  ~FlatMap()
  {
    // Destroy all elements.
    IndexType index = this->nextValidIndex(m_metadata, NO_MATCH);
    while(index < bucket_count())
    {
      m_buckets[index].get().~KeyValuePair();
      index = this->nextValidIndex(m_metadata, index);
    }

    // Unlike in clear() we don't need to reset metadata here.
  }

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

  iterator end() { return iterator(this, bucket_count()); }
  const_iterator end() const { return const_iterator(this, bucket_count()); }
  const_iterator cend() const { return const_iterator(this, bucket_count()); }

  // Capacity
  bool empty() const { return m_size == 0; }
  IndexType size() const { return m_size; }

  // Lookup
  iterator find(const KeyType& key);
  const_iterator find(const KeyType& key) const;

  ValueType& at(const KeyType& key)
  {
    auto it = this->find(key);
    if(it == end())
    {
      throw std::out_of_range {"axom::FlatMap::find(): key not found"};
    }
    return it->second;
  }
  const ValueType& at(const KeyType& key) const
  {
    auto it = this->find(key);
    if(it == end())
    {
      throw std::out_of_range {"axom::FlatMap::find(): key not found"};
    }
    return it->second;
  }

  ValueType& operator[](const KeyType& key)
  {
    return this->try_emplace(key).first->second;
  }
  const ValueType& operator[](const KeyType& key) const
  {
    return this->try_emplace(key).first->second;
  }

  IndexType count(const KeyType& key) const { return (find(key) != end()); }
  bool contains(const KeyType& key) const { return (find(key) != end()); }

  // Modifiers
  void clear()
  {
    // Destroy all elements.
    IndexType index = this->nextValidIndex(m_metadata, NO_MATCH);
    while(index < bucket_count())
    {
      m_buckets[index].get().~KeyValuePair();
      index = this->nextValidIndex(m_metadata, index);
    }

    // Also reset metadata.
    for(int group_index = 0; group_index < m_metadata.size(); group_index++)
    {
      m_metadata[group_index] = detail::flat_map::GroupBucket {};
    }
    m_metadata[m_metadata.size() - 1].setSentinel();
    m_size = 0;
    m_loadCount = 0;
  }

  std::pair<iterator, bool> insert(const value_type& value)
  {
    auto emplace_pos = getEmplacePos(value.first);
    emplaceImpl(emplace_pos, false, value);
    return emplace_pos;
  }
  std::pair<iterator, bool> insert(value_type&& value)
  {
    auto emplace_pos = getEmplacePos(value.first);
    emplaceImpl(emplace_pos, false, std::move(value));
    return emplace_pos;
  }
  template <typename InputPair>
  std::pair<iterator, bool> insert(InputPair&& pair)
  {
    return emplace(std::forward<InputPair>(pair));
  }
  template <typename... InputArgs>
  std::pair<iterator, bool> emplace(InputArgs&&... pair)
  {
    KeyValuePair kv {std::forward<InputArgs>(pair)...};
    auto emplace_pos = getEmplacePos(kv.first);
    emplaceImpl(emplace_pos, false, std::move(kv));
    return emplace_pos;
  }

  template <typename InputIt>
  void insert(InputIt first, InputIt last);

  template <typename... Args>
  std::pair<iterator, bool> insert_or_assign(const KeyType& key, Args&&... args)
  {
    auto emplace_pos = getEmplacePos(key);
    emplaceImpl(emplace_pos,
                true,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
    return emplace_pos;
  }
  template <typename... Args>
  std::pair<iterator, bool> insert_or_assign(KeyType&& key, Args&&... args)
  {
    auto emplace_pos = getEmplacePos(key);
    emplaceImpl(emplace_pos,
                true,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
    return emplace_pos;
  }

  template <typename... Args>
  std::pair<iterator, bool> try_emplace(const KeyType& key, Args&&... args)
  {
    auto emplace_pos = getEmplacePos(key);
    emplaceImpl(emplace_pos,
                false,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
    return emplace_pos;
  }
  template <typename... Args>
  std::pair<iterator, bool> try_emplace(KeyType&& key, Args&&... args)
  {
    auto emplace_pos = getEmplacePos(key);
    emplaceImpl(emplace_pos,
                false,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
    return emplace_pos;
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
  double load_factor() const { return ((double)m_loadCount) / bucket_count(); }
  double max_load_factor() const { return MAX_LOAD_FACTOR; }
  void rehash(IndexType count)
  {
    FlatMap rehashed(m_size,
                     std::make_move_iterator(begin()),
                     std::make_move_iterator(end()),
                     count);
    this->swap(rehashed);
  }
  void reserve(IndexType count) { rehash(std::ceil(count / MAX_LOAD_FACTOR)); }

private:
  template <typename InputIt>
  FlatMap(IndexType num_elems, InputIt first, InputIt last, IndexType bucket_count);

  std::pair<iterator, bool> getEmplacePos(const KeyType& key);

  template <typename... Args>
  void emplaceImpl(const std::pair<iterator, bool>& pos,
                   bool assign_on_existence,
                   Args&&... args);

  constexpr static IndexType MIN_NUM_BUCKETS {29};

  IndexType m_numGroups2;  // Number of groups of 15 buckets, expressed as a power of 2
  IndexType m_size;
  axom::Array<detail::flat_map::GroupBucket> m_metadata;

  // Storage details:
  struct alignas(KeyValuePair) PairStorage
  {
    unsigned char data[sizeof(KeyValuePair)];

    const KeyValuePair& get() const
    {
      return *(reinterpret_cast<const KeyValuePair*>(&data));
    }

    KeyValuePair& get() { return *(reinterpret_cast<KeyValuePair*>(&data)); }
  };

  axom::Array<PairStorage> m_buckets;

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
    assert(m_internalIdx >= 0 && m_internalIdx <= m_map->bucket_count());
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

  reference operator*() const { return m_map->m_buckets[m_internalIdx].get(); }

  pointer operator->() const
  {
    return &(m_map->m_buckets[m_internalIdx].get());
  }

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
                     if(this->m_buckets[bucket_index].get().first == key)
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
                     if(this->m_buckets[bucket_index].get().first == key)
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
auto FlatMap<KeyType, ValueType, Hash>::getEmplacePos(const KeyType& key)
  -> std::pair<iterator, bool>
{
  auto hash = MixedHash {}(key);
  // Resize to double the number of bucket groups if insertion would put us
  // above the maximum load factor.
  if(((m_loadCount + 1) / (double)bucket_count()) >= MAX_LOAD_FACTOR)
  {
    IndexType newNumGroups = m_metadata.size() * 2;
    rehash(newNumGroups * BucketsPerGroup - 1);
  }

  bool keyExistsAlready = false;
  IndexType foundBucketIndex = NO_MATCH;
  auto FindExistingElem = [&, this](IndexType bucket_index) -> bool {
    if(this->m_buckets[bucket_index].get().first == key)
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
  return {keyIterator, !keyExistsAlready};
}

template <typename KeyType, typename ValueType, typename Hash>
template <typename... Args>
void FlatMap<KeyType, ValueType, Hash>::emplaceImpl(
  const std::pair<iterator, bool>& pos,
  bool assign_on_existence,
  Args&&... args)
{
  IndexType bucketIndex = pos.first.m_internalIdx;
  bool keyExistsAlready = !pos.second;

  if(!keyExistsAlready)
  {
    new(&m_buckets[bucketIndex].data) KeyValuePair(std::forward<Args>(args)...);
  }
  else if(keyExistsAlready && assign_on_existence)
  {
    KeyValuePair kv(std::forward<Args>(args)...);
    m_buckets[bucketIndex].get().second = std::move(kv.second);
  }
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
