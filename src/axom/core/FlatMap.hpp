// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_FlatMap_HPP
#define Axom_Core_FlatMap_HPP

#include <tuple>
#include <type_traits>
#include <utility>
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/detail/FlatTable.hpp"

namespace axom
{
/*!
 * \class FlatMap
 *
 * \brief Provides a generic associative key-value container.
 *
 *  The FlatMap class is a container which maps unique keys to a single value.
 *  It supports insertion, removal, and lookup of key-value pairs in amortized
 *  constant time.
 *
 * \note FlatMap is designed to be a largely drop-in replacement for
 *  std::unordered_map. However, FlatMap is internally represented as an open-
 *  addressing, quadratic probing hash map; thus, the API differs from
 *  unordered_map in certain regards:
 *
 *   - Operations which insert a new element may invalidate references and
 *     iterators to existing elements in the container, if a resize of the
 *     underlying array is triggered.
 *   - Methods which only make sense for a closed-addressing hash map, such as
 *     begin/end(bucket_index), or bucket(key), are not implemented.
 *
 * \tparam KeyType the type of the keys to hold
 * \tparam ValueType the type of the values to hold
 * \tparam Hash the hash to use with the key type
 *
 * \pre KeyType must be EqualityComparable
 * \pre Hash is invocable with an instance of KeyType, and returns an integer
 *  value (32- or 64-bit)
 */
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
  using key_type = KeyType;
  using mapped_type = ValueType;
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
   *
   * \pre KeyType must be copy-constructible
   * \pre ValueType must be copy-constructible
   */
  FlatMap(const FlatMap& other)
    : m_numGroups2(other.m_numGroups2)
    , m_size(other.m_size)
    , m_metadata(other.m_metadata)
    , m_buckets(other.m_buckets.size())
    , m_loadCount(other.m_loadCount)
  {
    static_assert(std::is_copy_constructible<KeyType>::value,
                  "Cannot copy an axom::FlatMap when key type is not "
                  "copy-constructible.");
    static_assert(std::is_copy_constructible<ValueType>::value,
                  "Cannot copy an axom::FlatMap when value type is not "
                  "copy-constructible.");
    // Copy all elements.
    const auto metadata = m_metadata.view();
    IndexType index = this->nextValidIndex(metadata, NO_MATCH);
    while(index < bucket_count())
    {
      new(&m_buckets[index].data) KeyValuePair(other.m_buckets[index].get());
      index = this->nextValidIndex(metadata, index);
    }
  }

  /*!
   * \brief Copy assignment operator for a FlatMap instance.
   *
   * \param other the FlatMap to copy data from
   *
   * \pre KeyType must be copy-constructible
   * \pre ValueType must be copy-constructible
   */
  FlatMap& operator=(const FlatMap& other)
  {
    static_assert(std::is_copy_constructible<KeyType>::value,
                  "Cannot copy an axom::FlatMap when key type is not "
                  "copy-constructible.");
    static_assert(std::is_copy_constructible<ValueType>::value,
                  "Cannot copy an axom::FlatMap when value type is not "
                  "copy-constructible.");
    if(*this != other)
    {
      FlatMap new_map(other);
      swap(new_map);
    }
    return *this;
  }

  /// \brief Destructor for a FlatMap instance.
  ~FlatMap()
  {
    // Destroy all elements.
    const auto metadata = m_metadata.view();
    IndexType index = this->nextValidIndex(metadata, NO_MATCH);
    while(index < bucket_count())
    {
      m_buckets[index].get().~KeyValuePair();
      index = this->nextValidIndex(metadata, index);
    }

    // Unlike in clear() we don't need to reset metadata here.
  }

  /*!
   * \brief Swaps the contents of one FlatMap with another.
   */
  void swap(FlatMap& other)
  {
    axom::utilities::swap(m_numGroups2, other.m_numGroups2);
    axom::utilities::swap(m_size, other.m_size);
    m_metadata.swap(other.m_metadata);
    m_buckets.swap(other.m_buckets);
    axom::utilities::swap(m_loadCount, other.m_loadCount);
  }

  /*!
   * \brief Returns an iterator to the first valid object in the bucket array.
   */
  /// @{
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
  /// @}

  /*!
   * \brief Returns an iterator to "one past" the last valid object in the
   *  bucket array.
   */
  /// @{
  iterator end() { return iterator(this, bucket_count()); }
  const_iterator end() const { return const_iterator(this, bucket_count()); }
  const_iterator cend() const { return const_iterator(this, bucket_count()); }
  /// @}

  /*!
   * \brief Returns true if there are no entries in the FlatMap, false
   *  otherwise.
   */
  bool empty() const { return m_size == 0; }

  /*!
   * \brief Returns the number of entries stored in the FlatMap.
   */
  IndexType size() const { return m_size; }

  /*!
   * \brief Try to find an entry with a given key.
   *
   * \param [in] key the key to search for
   *
   * \return An iterator pointing to the corresponding key-value pair, or end()
   *  if the key wasn't found.
   */
  /// @{
  iterator find(const KeyType& key);
  const_iterator find(const KeyType& key) const;
  /// @}

  /*!
   * \brief Try to find an entry with a given key.
   *
   * \param [in] key the key to search for
   *
   * \return A reference to the corresponding value.
   * \throw std::out_of_range if the key is not found.
   */
  /// @{
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
  /// @}

  /*!
   * \brief Find an entry with a given key.
   *
   *  If a corresponding value does not exist, a default value for the value
   *  type will be inserted for the given key.
   *
   * \param [in] key the key to search for
   *
   * \return A reference to the corresponding value.
   *
   * \pre ValueType is default-constructible
   */
  /// @{
  ValueType& operator[](const KeyType& key)
  {
    static_assert(std::is_default_constructible<ValueType>::value,
                  "Cannot use axom::FlatMap::operator[] when value type is not "
                  "default-constructible.");
    return this->try_emplace(key).first->second;
  }
  const ValueType& operator[](const KeyType& key) const
  {
    static_assert(std::is_default_constructible<ValueType>::value,
                  "Cannot use axom::FlatMap::operator[] when value type is not "
                  "default-constructible.");
    return this->try_emplace(key).first->second;
  }
  /// @}

  /*!
   * \brief Return the number of entries matching a given key.
   *
   *  This method will always return 0 or 1.
   *
   * \param [in] key the key to search for
   */
  IndexType count(const KeyType& key) const
  {
    return contains(key) ? IndexType {1} : IndexType {0};
  }

  /*!
   * \brief Return true if the FlatMap contains a key, false otherwise.
   *
   * \param [in] key the key to search for
   */
  bool contains(const KeyType& key) const { return (find(key) != end()); }

  /*!
   * \brief Erases all elements from the FlatMap.
   */
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

  /*!
   * \brief Inserts a key-value pair into the FlatMap.
   *
   *  If the key already exists in the FlatMap, insertion is skipped.
   *  Otherwise, the key-value mapping is inserted into the FlatMap.
   *
   * \param [in] value the key-value pair to insert
   *
   * \return A pair consisting of:
   *   - an iterator pointing to either the existing key-value pair, or the
   *     newly-inserted pair
   *   - true if a new pair was inserted, false otherwise
   */
  /// @{
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
  /// @}

  /*!
   * \brief Inserts a range of key-value pairs into the FlatMap.
   *
   *  If the key already exists in the FlatMap, insertion is skipped.
   *  Otherwise, the key-value mapping is inserted into the FlatMap.
   *
   * \param [in] first the beginning of the range of pairs
   * \param [in] last the end of the range of pairs
   */
  template <typename InputIt>
  void insert(InputIt first, InputIt last);

  /*!
   * \brief Inserts a key-value pair into the FlatMap.
   *
   *  If the key already exists, assigns the value to the existing key in the
   *  FlatMap.
   *
   * \param [in] key the key to insert or assign
   * \param [in] args arguments to construct the value with
   *
   * \return A pair consisting of:
   *   - an iterator pointing to either the existing key-value pair, or the
   *     newly-inserted pair
   *   - true if a new pair was inserted, false otherwise
   */
  /// {@
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
  /// @}

  /*!
   * \brief Inserts a key-value pair into the FlatMap.
   *
   *  If the key already exists in the FlatMap, insertion is skipped.
   *  Otherwise, the key-value mapping is inserted into the FlatMap.
   *
   *  Compared to emplace(), this method only moves-from the value arguments
   *  if the key does not exist; otherwise, the input arguments are left as-is.
   *
   * \param [in] key the key to insert or assign
   * \param [in] args arguments to construct the value with.
   *
   * \return A pair consisting of:
   *   - an iterator pointing to either the existing key-value pair, or the
   *     newly-inserted pair
   *   - true if a new pair was inserted, false otherwise
   */
  /// {@
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
  /// @}

  /*!
   * \brief Remove a key-value pair from the FlatMap, specified by iterator.
   *
   * \param pos the iterator pointing to the key-value pair to remove
   *
   * \return an iterator to the next valid entry after the removed entry
   */
  /// {@
  iterator erase(iterator pos) { return erase(const_iterator {pos}); }
  iterator erase(const_iterator pos);
  /// @}

  /*!
   * \brief Remove a key-value pair from the FlatMap, specified by key.
   *
   *  If the key doesn't exist in the FlatMap, does nothing.
   *
   * \param key the key to remove
   *
   * \return 1 if an entry was removed, 0 otherwise
   */
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

  /*!
   * \brief Returns the number of buckets allocated in the FlatMap.
   *
   *  The maximum number of elements that can be stored in the FlatMap without
   *  resizing and rehashing is bucket_count() * max_load_factor().
   */
  IndexType bucket_count() const { return m_buckets.size(); }

  /*!
   * \brief Returns the current load factor of the FlatMap.
   */
  double load_factor() const { return ((double)m_loadCount) / bucket_count(); }

  /*!
   * \brief Returns the maximum load factor of the FlatMap.
   */
  double max_load_factor() const { return MAX_LOAD_FACTOR; }

  /*!
   * \brief Explicitly rehash the FlatMap with a given number of buckets.
   *
   * \param count the minimum number of buckets to allocate for the rehash
   */
  void rehash(IndexType count)
  {
    FlatMap rehashed(m_size,
                     std::make_move_iterator(begin()),
                     std::make_move_iterator(end()),
                     count);
    this->swap(rehashed);
  }

  /*!
   * \brief Reallocate and rehash the FlatMap, such that up to the specified
   *  number of elements may be inserted without a rehash.
   *
   * \param count the number of elements to fit without a rehash
   */
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
  using PairStorage = detail::flat_map::TypeErasedStorage<KeyValuePair>;
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
  IndexType minBuckets = MIN_NUM_BUCKETS;
  bucket_count = axom::utilities::max(minBuckets, bucket_count);
  // Get the smallest power-of-two number of groups satisfying:
  // N * GroupSize - 1 >= minBuckets
  // TODO: we should add a countl_zero overload for 64-bit integers
  {
    std::int32_t numGroups =
      std::ceil((bucket_count + 1) / (double)BucketsPerGroup);
    m_numGroups2 = 31 - (axom::utilities::countl_zero(numGroups));
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

  // If the key already exists, return the existing iterator.
  iterator existing_elem = this->find(key);
  if(existing_elem != this->end())
  {
    return {existing_elem, false};
  }
  // Resize to double the number of bucket groups if insertion would put us
  // above the maximum load factor.
  if(((m_loadCount + 1) / (double)bucket_count()) >= MAX_LOAD_FACTOR)
  {
    IndexType newNumGroups = m_metadata.size() * 2;
    rehash(newNumGroups * BucketsPerGroup - 1);
  }

  // Get an empty index to place the element into.
  IndexType newBucket = this->probeEmptyIndex(m_numGroups2, m_metadata, hash);

  // Add a hash to the corresponding bucket slot.
  this->setBucketHash(m_metadata, newBucket, hash);
  m_size++;
  m_loadCount++;

  // Return an iterator pointing to the just-found empty bucket.
  iterator keyIterator = iterator(this, newBucket);
  return {keyIterator, true};
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
