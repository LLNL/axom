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
