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

TEST(core_flatmap, default_init)
{
  axom::FlatMap<int, std::string> int_to_dbl;
  EXPECT_EQ(0, int_to_dbl.size());
  EXPECT_EQ(true, int_to_dbl.empty());
}

TEST(core_flatmap, insert_only)
{
  axom::FlatMap<int, double> int_to_dbl;

  int_to_dbl.insert({0, 10.0});
  EXPECT_EQ(1, int_to_dbl.size());

  int_to_dbl.insert({1, 20.0});
  EXPECT_EQ(2, int_to_dbl.size());

  int_to_dbl.insert({2, 30.0});
  EXPECT_EQ(3, int_to_dbl.size());

  // Check consistency of added values.
  const double expected_str[3] {10.0, 20.0, 30.0};
  for(int i = 0; i < 3; i++)
  {
    auto iterator = int_to_dbl.find(i);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, i);
    EXPECT_EQ(iterator->second, expected_str[i]);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    double value = int_to_dbl[i];
    EXPECT_EQ(value, expected_str[i]);
    EXPECT_EQ(int_to_dbl.size(), 3);
  }
}

TEST(core_flatmap_str, insert_only)
{
  axom::FlatMap<std::string, std::string> int_to_dbl;

  int_to_dbl.insert({std::to_string(0), std::to_string(10.0)});
  EXPECT_EQ(1, int_to_dbl.size());

  int_to_dbl.insert({std::to_string(1), std::to_string(20.0)});
  EXPECT_EQ(2, int_to_dbl.size());

  int_to_dbl.insert({std::to_string(2), std::to_string(30.0)});
  EXPECT_EQ(3, int_to_dbl.size());

  // Check consistency of added values.
  const double expected_str[3] {10.0, 20.0, 30.0};
  for(int i = 0; i < 3; i++)
  {
    std::string key = std::to_string(i);
    std::string expected_value = std::to_string(expected_str[i]);
    auto iterator = int_to_dbl.find(key);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, key);
    EXPECT_EQ(iterator->second, expected_value);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    std::string value = int_to_dbl[key];
    EXPECT_EQ(value, expected_value);
    EXPECT_EQ(int_to_dbl.size(), 3);
  }
}

TEST(core_flatmap, initializer_list)
{
  axom::FlatMap<int, double> int_to_dbl {{0, 10.0}, {1, 20.0}, {2, 30.0}};

  EXPECT_EQ(3, int_to_dbl.size());

  // Check consistency of added values.
  const double expected_str[3] {10.0, 20.0, 30.0};
  for(int i = 0; i < 3; i++)
  {
    auto iterator = int_to_dbl.find(i);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, i);
    EXPECT_EQ(iterator->second, expected_str[i]);

    // Using operator[] with an already-existing key should return the
    // existing value and not add a value.
    double value = int_to_dbl[i];
    EXPECT_EQ(value, expected_str[i]);
    EXPECT_EQ(int_to_dbl.size(), 3);
  }
}

TEST(core_flatmap, index_operator_default)
{
  axom::FlatMap<int, double> int_to_dbl;

  int NUM_ELEMS = 10;

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    double default_value = int_to_dbl[i];
    EXPECT_EQ(default_value, 0);
    int_to_dbl[i] = i + 10.0;
  }

  EXPECT_EQ(NUM_ELEMS, int_to_dbl.size());

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto iterator = int_to_dbl.find(i);
    EXPECT_EQ(iterator->second, i + 10.0);
  }
}

TEST(core_flatmap, init_and_clear)
{
  axom::FlatMap<int, double> int_to_dbl;

  // Insert enough elements to trigger a resize of the buckets.
  // This allows us to test that a clear() doesn't reset the allocated buckets.
  int NUM_ELEMS_RESIZE = 40;
  EXPECT_GT(NUM_ELEMS_RESIZE, int_to_dbl.bucket_count());

  for(int i = 0; i < NUM_ELEMS_RESIZE; i++)
  {
    int_to_dbl[i] = i + 10.0;
  }

  EXPECT_EQ(NUM_ELEMS_RESIZE, int_to_dbl.size());

  int buckets_before_clear = int_to_dbl.bucket_count();

  int_to_dbl.clear();

  EXPECT_EQ(int_to_dbl.size(), 0);
  EXPECT_EQ(int_to_dbl.load_factor(), 0.0);
  EXPECT_EQ(int_to_dbl.bucket_count(), buckets_before_clear);
  for(int i = 0; i < 3; i++)
  {
    auto iterator = int_to_dbl.find(i);
    EXPECT_EQ(iterator, int_to_dbl.end());
  }
}

TEST(core_flatmap, insert_until_rehash)
{
  axom::FlatMap<int, double> int_to_dbl;

  const int INIT_CAPACITY = int_to_dbl.bucket_count();
  const double LOAD_FACTOR = int_to_dbl.max_load_factor();
  const int SIZE_NO_REHASH = LOAD_FACTOR * INIT_CAPACITY;

  for(int i = 0; i < SIZE_NO_REHASH; i++)
  {
    int_to_dbl.insert({i, 2. * i + 1});
  }
  EXPECT_EQ(int_to_dbl.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(int_to_dbl.size(), SIZE_NO_REHASH);

  // Next insert should trigger a rehash.
  int_to_dbl.insert({SIZE_NO_REHASH, 2. * SIZE_NO_REHASH + 1});
  EXPECT_GT(int_to_dbl.bucket_count(), INIT_CAPACITY);
  EXPECT_EQ(int_to_dbl.size(), SIZE_NO_REHASH + 1);

  // Check consistency of values.
  for(int i = 0; i < SIZE_NO_REHASH + 1; i++)
  {
    auto iterator = int_to_dbl.find(i);
    EXPECT_NE(iterator, int_to_dbl.end());
    EXPECT_EQ(iterator->first, i);
    EXPECT_EQ(iterator->second, 2. * i + 1);
  }
}

TEST(core_flatmap, insert_then_delete)
{
  axom::FlatMap<int, double> int_to_dbl;

  const int INIT_CAPACITY = int_to_dbl.bucket_count();
  const double LOAD_FACTOR = int_to_dbl.max_load_factor();
  const int NUM_INSERTS = LOAD_FACTOR * INIT_CAPACITY * 4;

  for(int i = 0; i < NUM_INSERTS; i++)
  {
    int_to_dbl.insert({i, 2. * i + 1});
  }
  EXPECT_EQ(int_to_dbl.size(), NUM_INSERTS);
  EXPECT_GE(int_to_dbl.bucket_count(), NUM_INSERTS);

  for(int i = 0; i < NUM_INSERTS; i += 3)
  {
    // Delete every third entry starting from 0, inclusive.
    // (i.e. keys 0, 3, 6, ...)
    int_to_dbl.erase(i);
  }

  // Check consistency of values.
  for(int i = 0; i < NUM_INSERTS; i++)
  {
    auto iterator = int_to_dbl.find(i);
    if(i % 3 == 0)
    {
      EXPECT_EQ(iterator, int_to_dbl.end());
      EXPECT_EQ(0, int_to_dbl.count(i));
      EXPECT_EQ(false, int_to_dbl.contains(i));
    }
    else
    {
      EXPECT_NE(iterator, int_to_dbl.end());
      EXPECT_EQ(iterator->first, i);
      EXPECT_EQ(iterator->second, 2. * i + 1);
      EXPECT_EQ(1, int_to_dbl.count(i));
      EXPECT_EQ(true, int_to_dbl.contains(i));
    }
  }
}
