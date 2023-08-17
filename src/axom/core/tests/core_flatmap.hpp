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
