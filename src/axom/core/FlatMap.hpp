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
  using size_type = IndexType;
  using value_type = std::pair<const KeyType, ValueType>;
  using iterator = void;
  using const_iterator = void;

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
  bool empty() const;
  IndexType size() const;

  // Lookup
  iterator find(const KeyType& key);
  const_iterator find(const KeyType& key) const;

  ValueType& at(const KeyType& key);
  const ValueType& at(const KeyType& key) const;

  ValueType& operator[](const KeyType& key);
  const ValueType& operator[](const KeyType& key) const;

  IndexType count(const KeyType& key) const;
  bool contains(const KeyType& key) const;

  // Modifiers
  void clear();
  std::pair<iterator, bool> insert(const value_type& value);
  std::pair<iterator, bool> insert(value_type&& value);
  template <typename InputPair>
  std::pair<iterator, bool> insert(InputPair&& pair);
  template <typename InputPair>
  std::pair<iterator, bool> emplace(InputPair&& pair);
  template <typename InputIt>
  void insert(InputIt first, InputIt last);

  // Hashing
  IndexType bucket_count() const;
  double load_factor() const;
  double max_load_factor() const;
  void rehash(IndexType count);
  void reserve(IndexType count);
};

}  // namespace axom

#endif  // Axom_Core_FlatMap_HPP
