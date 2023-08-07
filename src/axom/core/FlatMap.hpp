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

  struct KeyValuePair
  {
    const KeyType first;
    ValueType second;
  };

  template <bool Const>
  friend class IteratorImpl;

public:
  using size_type = IndexType;
  using value_type = KeyValuePair;
  using iterator = IteratorImpl<false>;
  using const_iterator = IteratorImpl<true>;

  // Constructors
  FlatMap();

  explicit FlatMap(IndexType bucket_count);

  template <typename InputIt>
  FlatMap(InputIt first, InputIt last, IndexType bucket_count = -1);

  explicit FlatMap(std::initializer_list<value_type> init,
                   IndexType bucket_count = -1);

  // Iterators
  iterator begin();
  const_iterator begin() const;
  const_iterator cbegin() const;

  iterator end();
  const_iterator end() const;
  const_iterator cend() const;

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
    return emplace(value);
  }
  std::pair<iterator, bool> insert(value_type&& value)
  {
    return emplace(std::move(value));
  }
  template <typename InputPair>
  std::pair<iterator, bool> insert(InputPair&& pair)
  {
    return emplace(std::forward<InputPair>(pair));
  }
  template <typename... InputArgs>
  std::pair<iterator, bool> emplace(InputArgs&&... pair)
  {
    KeyValuePair kv {pair...};
    emplaceImpl(std::move(kv.first), std::move(kv.second));
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

  // Hashing
  IndexType bucket_count() const;
  double load_factor() const;
  double max_load_factor() const;
  void rehash(IndexType count);
  void reserve(IndexType count);

private:
  template <typename... Args>
  std::pair<iterator, bool> emplaceImpl(bool assign_on_existence,
                                        KeyType&& key,
                                        Args&&... args);

  IndexType m_numGroups2;  // Number of groups of 15 buckets, expressed as a power of 2
  IndexType m_size;
  axom::Array<detail::flat_map::GroupBucket> m_metadata;
  axom::Array<KeyValuePair> m_buckets;
};

template <typename KeyType, typename ValueType, typename Hash>
template <bool Const>
class FlatMap<KeyType, ValueType, Hash>::IteratorImpl
{
private:
  using MapType = FlatMap<KeyType, ValueType, Hash>;

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
    // TODO
  }

  IteratorImpl operator++(int)
  {
    // TODO
  }

  reference operator*() const { return m_map->m_buckets[m_internalIdx]; }

  pointer operator->() const { return &(m_map->m_buckets[m_internalIdx]); }

private:
  MapConstType* m_map;
  IndexType m_internalIdx;
};

template <typename KeyType, typename ValueType, typename Hash>
auto FlatMap<KeyType, ValueType, Hash>::find(const KeyType& key) -> iterator
{
  auto hash = Hash {}(key);
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
  auto hash = Hash {}(key);
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
template <typename... InputArgs>
auto FlatMap<KeyType, ValueType, Hash>::emplaceImpl(bool assign_on_existence,
                                                    KeyType&& key,
                                                    InputArgs&&... args)
  -> std::pair<iterator, bool>
{
  auto hash = Hash {}(key);
  // TODO: ensure that we have enough space to insert with given load factor

  bool keyExistsAlready = false;
  IndexType foundBucketIndex = NO_MATCH;
  auto FindExistingElem = [&, this](IndexType bucket_index) -> bool {
    if(this->m_buckets[bucket_index].first == key)
    {
      keyExistsAlready = true;
      foundBucketIndex = bucket_index;
      // Exit out of probing, we can't insert if the key already exists.
      return false;
    }
    return true;
  };

  IndexType newBucket =
    this->probeEmptyIndex(m_numGroups2, m_metadata, hash, FindExistingElem);
  if(!keyExistsAlready)
  {
    foundBucketIndex = newBucket;
    // Add a hash to the corresponding bucket slot.
    this->setBucketHash(m_metadata, newBucket, hash);
  }
  iterator keyIterator = iterator(this, foundBucketIndex);
  if(!keyExistsAlready)
  {
    new(&m_buckets[foundBucketIndex])
      KeyValuePair(std::move(key), std::forward<InputArgs>(args)...);
  }
  else if(keyExistsAlready && assign_on_existence)
  {
    m_buckets[foundBucketIndex] =
      KeyValuePair(std::move(key), std::forward<InputArgs>(args)...);
  }
  return {keyIterator, !keyExistsAlready};
}

}  // namespace axom

#endif  // Axom_Core_FlatMap_HPP
