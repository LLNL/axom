// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/FlatMap.hpp"

// gtest includes
#include "gtest/gtest.h"

inline void flatmap_get_value(double key, std::string& out)
{
  out = std::to_string(key);
}

inline void flatmap_get_value(int key, std::string& out)
{
  out = std::to_string(key);
}

template <typename T, typename U>
inline void flatmap_get_value(T key, U& out)
{
  out = key;
}

template <typename FlatMapType>
class core_flatmap : public ::testing::Test
{
public:
  using MapType = FlatMapType;
  using KeyType = typename FlatMapType::key_type;
  using ValueType = typename FlatMapType::mapped_type;

  template <typename T>
  KeyType getKey(T input)
  {
    KeyType key;
    flatmap_get_value(input, key);
    return key;
  }

  template <typename T>
  ValueType getValue(T input)
  {
    ValueType val;
    flatmap_get_value(input, val);
    return val;
  }

  ValueType getDefaultValue() { return ValueType(); }
};

using MyTypes = ::testing::Types<axom::FlatMap<int, double>,
                                 axom::FlatMap<int, std::string>,
                                 axom::FlatMap<std::string, double>,
                                 axom::FlatMap<std::string, std::string>>;

TYPED_TEST_SUITE(core_flatmap, MyTypes);

AXOM_TYPED_TEST(core_flatmap, default_init)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;
  EXPECT_EQ(0, int_to_dbl.size());
  EXPECT_EQ(true, int_to_dbl.empty());
}

AXOM_TYPED_TEST(core_flatmap, insert_only)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  const int NUM_ELEMS = 100;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);
    // Initial insertion of a given key should succeed.
    auto initial_insert = int_to_dbl.insert({key, value});
    EXPECT_EQ(int_to_dbl.size(), i + 1);
    EXPECT_EQ(initial_insert.first, int_to_dbl.find(key));
    EXPECT_EQ(value, int_to_dbl.at(key));
    EXPECT_TRUE(initial_insert.second);

    int current_bucket_capacity = int_to_dbl.bucket_count();

    // Inserting a duplicate key should not change the value.
    auto value_dup = this->getValue(i * 10.0 + 5.0);
    auto duplicate_insert = int_to_dbl.insert({key, value_dup});
    EXPECT_EQ(int_to_dbl.size(), i + 1);
    EXPECT_EQ(duplicate_insert.first, int_to_dbl.find(key));
    EXPECT_EQ(value, int_to_dbl.at(key));
    EXPECT_FALSE(duplicate_insert.second);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    auto value_indexed = int_to_dbl[key];
    EXPECT_EQ(value_indexed, value);
    EXPECT_EQ(int_to_dbl.size(), i + 1);

    // Check that a rehash didn't occur on the second insertion.
    EXPECT_EQ(duplicate_insert.first, initial_insert.first);
    EXPECT_EQ(current_bucket_capacity, int_to_dbl.bucket_count());
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_or_assign)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  const int NUM_ELEMS = 100;

  // Test insert behavior of FlatMap::insert_or_assign.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    auto result = int_to_dbl.insert_or_assign(key, value);
    EXPECT_EQ(i + 1, int_to_dbl.size());
    EXPECT_EQ(value, int_to_dbl.at(key));
    EXPECT_EQ(result.first, int_to_dbl.find(key));
    EXPECT_TRUE(result.second);
  }

  // Test assign behavior of FlatMap::insert_or_assign.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 7.0);

    auto result = int_to_dbl.insert_or_assign(key, value);
    EXPECT_EQ(value, int_to_dbl.at(key));
    EXPECT_EQ(result.first, int_to_dbl.find(key));
    EXPECT_FALSE(result.second);
  }

  // Assignments should not change size of FlatMap.
  EXPECT_EQ(NUM_ELEMS, int_to_dbl.size());
}

TEST(core_flatmap_moveonly, try_emplace)
{
  axom::FlatMap<int, std::unique_ptr<double>> int_to_dbl;

  const int NUM_ELEMS = 40;

  // Test behavior when key does not exist.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    std::unique_ptr<double> value {new double {i + 10.0}};
    auto result = int_to_dbl.try_emplace(i, std::move(value));
    EXPECT_EQ(*(int_to_dbl[i]), i + 10.0);
    EXPECT_EQ(result.first, int_to_dbl.find(i));
    EXPECT_TRUE(result.second);
    // Value should have been moved.
    EXPECT_EQ(value.get(), nullptr);
  }

  // Test behavior when key already exists.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    std::unique_ptr<double> value {new double {i + 20.0}};
    auto result = int_to_dbl.try_emplace(i, std::move(value));
    EXPECT_EQ(*(int_to_dbl[i]), i + 10.0);
    EXPECT_EQ(result.first, int_to_dbl.find(i));
    EXPECT_FALSE(result.second);
    // Since key already exists, value should NOT be moved.
    EXPECT_NE(value.get(), nullptr);
    EXPECT_EQ(*value, i + 20.0);
  }
}

AXOM_TYPED_TEST(core_flatmap, initializer_list)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl {{this->getKey(0), this->getValue(10.0)},
                      {this->getKey(1), this->getValue(20.0)},
                      {this->getKey(2), this->getValue(30.0)}};

  EXPECT_EQ(3, int_to_dbl.size());

  // Check consistency of added values.
  const double expected_str[3] {10.0, 20.0, 30.0};
  for(int i = 0; i < 3; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(expected_str[i]);

    auto iterator = int_to_dbl.find(key);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, key);
    EXPECT_EQ(iterator->second, value);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    auto indexed_value = int_to_dbl[key];
    EXPECT_EQ(indexed_value, value);
    EXPECT_EQ(int_to_dbl.size(), 3);
  }
}

AXOM_TYPED_TEST(core_flatmap, index_operator_default)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  const int NUM_ELEMS = 100;

  auto expected_default_value = this->getDefaultValue();
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto default_value = int_to_dbl[key];

    EXPECT_EQ(default_value, expected_default_value);

    auto new_value = this->getValue(i * 10.0 + 5.0);
    int_to_dbl[key] = new_value;
  }

  EXPECT_EQ(NUM_ELEMS, int_to_dbl.size());

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    auto iterator = int_to_dbl.find(key);
    EXPECT_EQ(iterator->second, value);
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_clear)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  // Insert enough elements to trigger a resize of the buckets.
  // This allows us to test that a clear() doesn't reset the allocated buckets.
  int NUM_ELEMS_RESIZE = 100;
  EXPECT_GT(NUM_ELEMS_RESIZE, int_to_dbl.bucket_count());

  for(int i = 0; i < NUM_ELEMS_RESIZE; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    int_to_dbl[key] = value;
  }

  EXPECT_EQ(NUM_ELEMS_RESIZE, int_to_dbl.size());

  int buckets_before_clear = int_to_dbl.bucket_count();

  int_to_dbl.clear();

  EXPECT_EQ(int_to_dbl.size(), 0);
  EXPECT_EQ(int_to_dbl.load_factor(), 0.0);
  EXPECT_EQ(int_to_dbl.bucket_count(), buckets_before_clear);
  for(int i = 0; i < NUM_ELEMS_RESIZE; i++)
  {
    auto key = this->getKey(i);
    auto iterator = int_to_dbl.find(key);
    EXPECT_EQ(iterator, int_to_dbl.end());
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_move)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    int_to_dbl[key] = value;
  }

  MapType moved_to_map = std::move(int_to_dbl);

  EXPECT_EQ(int_to_dbl.size(), 0);
  EXPECT_EQ(int_to_dbl.load_factor(), 0);
  EXPECT_EQ(moved_to_map.size(), NUM_ELEMS);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);
    EXPECT_EQ(moved_to_map[key], value);

    auto old_it = int_to_dbl.find(key);
    EXPECT_EQ(old_it, int_to_dbl.end());
  }
}

TEST(core_flatmap_moveonly, init_and_move_moveonly)
{
  axom::FlatMap<int, std::unique_ptr<double>> int_to_dbl;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    int_to_dbl.emplace(i, new double {i + 10.0});
  }

  axom::FlatMap<int, std::unique_ptr<double>> int_to_dbl_move =
    std::move(int_to_dbl);

  EXPECT_EQ(int_to_dbl.size(), 0);
  EXPECT_EQ(int_to_dbl.load_factor(), 0);
  EXPECT_EQ(int_to_dbl_move.size(), NUM_ELEMS);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(*(int_to_dbl_move[i]), i + 10.0);

    auto old_it = int_to_dbl.find(i);
    EXPECT_EQ(old_it, int_to_dbl.end());
  }
}

AXOM_TYPED_TEST(core_flatmap, init_and_copy)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;
  int NUM_ELEMS = 40;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    int_to_dbl[key] = value;
  }

  int expected_buckets = int_to_dbl.bucket_count();

  MapType int_to_dbl_copy = int_to_dbl;

  EXPECT_EQ(int_to_dbl.size(), NUM_ELEMS);
  EXPECT_EQ(int_to_dbl.bucket_count(), expected_buckets);
  EXPECT_EQ(int_to_dbl_copy.size(), NUM_ELEMS);
  EXPECT_EQ(int_to_dbl_copy.bucket_count(), expected_buckets);
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i + 10.0);

    EXPECT_EQ(int_to_dbl[key], value);
    EXPECT_EQ(int_to_dbl_copy[key], value);
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_until_rehash)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  const int INIT_CAPACITY = int_to_dbl.bucket_count();
  const double LOAD_FACTOR = int_to_dbl.max_load_factor();
  const int SIZE_NO_REHASH = LOAD_FACTOR * INIT_CAPACITY;

  for(int i = 0; i < SIZE_NO_REHASH; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    int_to_dbl.insert({key, value});
  }
  EXPECT_EQ(int_to_dbl.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(int_to_dbl.size(), SIZE_NO_REHASH);

  // Next insert should trigger a rehash.
  {
    auto key_rehash = this->getKey(SIZE_NO_REHASH);
    auto value_rehash = this->getValue(2. * SIZE_NO_REHASH + 1);

    int_to_dbl.insert({key_rehash, value_rehash});
  }
  EXPECT_GT(int_to_dbl.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(int_to_dbl.size(), SIZE_NO_REHASH + 1);

  // Check consistency of values.
  for(int i = 0; i < SIZE_NO_REHASH + 1; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    auto iterator = int_to_dbl.find(key);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, key);
    EXPECT_EQ(iterator->second, value);
  }
}

AXOM_TYPED_TEST(core_flatmap, insert_then_delete)
{
  using MapType = typename TestFixture::MapType;
  MapType int_to_dbl;

  const int INIT_CAPACITY = int_to_dbl.bucket_count();
  const double LOAD_FACTOR = int_to_dbl.max_load_factor();
  const int NUM_INSERTS = LOAD_FACTOR * INIT_CAPACITY * 4;

  for(int i = 0; i < NUM_INSERTS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    int_to_dbl.insert({key, value});
  }
  EXPECT_EQ(int_to_dbl.size(), NUM_INSERTS);
  EXPECT_GE(int_to_dbl.bucket_count(), NUM_INSERTS);

  for(int i = 0; i < NUM_INSERTS; i += 3)
  {
    auto key = this->getKey(i);

    auto iterator_to_remove = int_to_dbl.find(key);
    auto one_after_elem = iterator_to_remove;
    one_after_elem++;
    // Delete every third entry starting from 0, inclusive.
    // (i.e. keys 0, 3, 6, ...)
    auto deleted_iterator = int_to_dbl.erase(iterator_to_remove);

    EXPECT_EQ(deleted_iterator, one_after_elem);
  }

  // Check consistency of values.
  for(int i = 0; i < NUM_INSERTS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(2. * i + 1);

    auto iterator = int_to_dbl.find(key);
    if(i % 3 == 0)
    {
      EXPECT_EQ(iterator, int_to_dbl.end());
      EXPECT_EQ(0, int_to_dbl.count(key));
      EXPECT_EQ(false, int_to_dbl.contains(key));
    }
    else
    {
      EXPECT_NE(iterator, int_to_dbl.end());
      EXPECT_EQ(iterator->first, key);
      EXPECT_EQ(iterator->second, value);
      EXPECT_EQ(1, int_to_dbl.count(key));
      EXPECT_EQ(true, int_to_dbl.contains(key));
    }
  }
}
