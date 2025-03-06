// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/Sorting.hpp"

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

TEST(core_utilities, qsort_sort_double)
{
  constexpr int NUM_ITERS = 1000;
  constexpr int MAX_SIZE = 100;

  // Run this for a few iterations to test that sorting works on different random shuffles.
  for(int n = 1; n < MAX_SIZE; ++n)
  {
    std::vector<double> vec;
    vec.resize(n);

    for(int iter = 0; iter < NUM_ITERS; iter++)
    {
      // Fill it with random doubles in the range of [0, 1024)
      std::generate_n(vec.data(), n, []() -> double {
        return axom::utilities::random_real(0., 1024.);
      });
      
      // sort the first two iterations in descending/ascending order
      if(iter == 0)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a > b; });
      }
      else if(iter == 1)
      {
        std::sort(vec.begin(), vec.end(), [](int a, int b) { return a < b; });
      }
      
      axom::utilities::Sorting<double, MAX_SIZE>::qsort(vec.data(), n);
    }

    // Check that the results are sorted
    for(int i = 1; i < n; i++)
    {
      EXPECT_TRUE(vec[i - 1] <= vec[i]) << "data not sorted!";
    }
  }
}
