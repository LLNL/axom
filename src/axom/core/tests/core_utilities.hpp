// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <algorithm>
#include <random>
#include <numeric>

TEST(core_utilities, insertion_sort_int)
{
  constexpr int NUM_ITERS = 100;
  constexpr int NUM_INTS = 32;
  // Create some RNG state used by std::shuffle
  std::random_device rd;
  std::mt19937 mersenne_twister(rd());

  // Run this for a few iterations to test that sorting works on different
  // random shuffles.
  for(int iter = 0; iter < NUM_ITERS; iter++)
  {
    // Create an array of 32 integers
    int data[NUM_INTS];

    // Fill it with uniform numbers from 0...31
    std::iota(data, data + NUM_INTS, 1);

    // Randomly shuffle the data.
    std::shuffle(data, data + NUM_INTS, mersenne_twister);

    // Sort the data.
    axom::utilities::insertionSort(data, NUM_INTS);

    // The result should be our ordered range 0...31
    for(int i = 0; i < NUM_INTS; i++)
    {
      ASSERT_EQ(data[i], i + 1);
    }
  }
}

TEST(core_utilities, insertion_sort_doubles)
{
  constexpr int NUM_ITERS = 100;
  constexpr int NUM_DBLS = 64;
  // Create some RNG state used by std::shuffle
  std::random_device rd;
  std::mt19937 mersenne_twister(rd());

  // Run this for a few iterations to test that sorting works on different
  // random shuffles.
  for(int iter = 0; iter < NUM_ITERS; iter++)
  {
    // Create an array of 64 doubles
    int data[NUM_DBLS];

    // Fill it with random doubles in the range of [0, 1024)
    std::generate_n(data, NUM_DBLS, []() -> double {
      return axom::utilities::random_real(0.0, 1024.0);
    });

    // Randomly shuffle the data.
    std::shuffle(data, data + NUM_DBLS, mersenne_twister);

    // Sort the data.
    axom::utilities::insertionSort(data, NUM_DBLS);

    // Check that the range is now sorted.
    for(int i = 0; i < NUM_DBLS - 1; i++)
    {
      ASSERT_LE(data[i], data[i + 1]);
    }
  }
}
