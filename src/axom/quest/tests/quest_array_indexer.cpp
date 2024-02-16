// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Array.hpp"
#include "axom/quest/ArrayIndexer.hpp"
#include "axom/fmt.hpp"

// Google test include
#include "gtest/gtest.h"

// Test strides and permutations.
TEST(quest_array_indexer, quest_strides_and_permutatations)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {2};

    axom::ArrayIndexer<axom::IndexType, DIM> colIndexer(lengths, 'c');
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {0};
    axom::StackArray<axom::IndexType, DIM> colStrides {1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }

    axom::ArrayIndexer<axom::IndexType, DIM> rowIndexer(lengths, 'r');
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {0};
    axom::StackArray<axom::IndexType, DIM> rowStrides {1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};

    axom::ArrayIndexer<axom::IndexType, DIM> colIndexer(lengths, 'c');
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {0, 1};
    axom::StackArray<axom::IndexType, DIM> colStrides {2, 1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }

    axom::ArrayIndexer<axom::IndexType, DIM> rowIndexer(lengths, 'r');
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {1, 0};
    axom::StackArray<axom::IndexType, DIM> rowStrides {1, 3};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
  }
  {
    // 3D
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};

    axom::ArrayIndexer<axom::IndexType, DIM> colIndexer(lengths, 'c');
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {0, 1, 2};
    axom::StackArray<axom::IndexType, DIM> colStrides {6, 2, 1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }

    axom::ArrayIndexer<axom::IndexType, DIM> rowIndexer(lengths, 'r');
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {2, 1, 0};
    axom::StackArray<axom::IndexType, DIM> rowStrides {1, 5, 15};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
  }
}

// Test column-major offsets.
TEST(quest_array_indexer, quest_col_major_offset)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {3};
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'c');
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      axom::StackArray<axom::IndexType, DIM> indices {i};
      EXPECT_EQ(ai.toFlatIndex(indices), offset);
      EXPECT_EQ(ai.toMultiIndex(offset), indices);
      ++offset;
    }
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'c');
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(ai.toFlatIndex(indices), offset);
        EXPECT_EQ(ai.toMultiIndex(offset), indices);
        ++offset;
      }
    }
  }
  {
    // 3D
    axom::StackArray<axom::IndexType, 3> lengths {5, 3, 2};
    axom::ArrayIndexer<axom::IndexType, 3> ai(lengths, 'c');
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        for(int k = 0; k < lengths[2]; ++k)
        {
          axom::StackArray<axom::IndexType, 3> indices {i, j, k};
          EXPECT_EQ(ai.toFlatIndex(indices), offset);
          EXPECT_EQ(ai.toMultiIndex(offset), indices);
          ++offset;
        }
      }
    }
  }
}

// Test row-major offsets.
TEST(quest_array_indexer, quest_row_major_offset)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {3};
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'r');
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      axom::StackArray<axom::IndexType, DIM> indices {i};
      EXPECT_EQ(ai.toFlatIndex(indices), offset);
      EXPECT_EQ(ai.toMultiIndex(offset), indices);
      ++offset;
    }
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'r');
    axom::IndexType offset = 0;
    for(int j = 0; j < lengths[1]; ++j)
    {
      for(int i = 0; i < lengths[0]; ++i)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(ai.toFlatIndex(indices), offset);
        EXPECT_EQ(ai.toMultiIndex(offset), indices);
        ++offset;
      }
    }
  }
  {
    // 3D
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'r');
    axom::IndexType offset = 0;
    for(int k = 0; k < lengths[2]; ++k)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        for(int i = 0; i < lengths[0]; ++i)
        {
          axom::StackArray<axom::IndexType, DIM> indices {i, j, k};
          EXPECT_EQ(ai.toFlatIndex(indices), offset);
          EXPECT_EQ(ai.toMultiIndex(offset), indices);
          ++offset;
        }
      }
    }
  }
}

template <int DIM>
typename std::enable_if<DIM == 1>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::ArrayIndexer<axom::IndexType, DIM>& ai)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> i;
  auto& n = i[fastestDirs[0]];
  axom::IndexType offset = 0;
  for(n = 0; n < lengths[fastestDirs[0]]; ++n)
  {
    SLIC_INFO(axom::fmt::format("offset {}   i {}", offset, i));
    EXPECT_EQ(ai.toMultiIndex(offset), i);
    EXPECT_EQ(ai.toFlatIndex(i), offset);
    ++offset;
  }
}

template <int DIM>
typename std::enable_if<DIM == 2>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::ArrayIndexer<axom::IndexType, DIM>& ai)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> ij;
  auto& m = ij[fastestDirs[1]];
  auto& n = ij[fastestDirs[0]];
  axom::IndexType offset = 0;
  for(m = 0; m < lengths[fastestDirs[1]]; ++m)
  {
    for(n = 0; n < lengths[fastestDirs[0]]; ++n)
    {
      SLIC_INFO(axom::fmt::format("offset {}   ij {}", offset, ij));
      EXPECT_EQ(ai.toMultiIndex(offset), ij);
      EXPECT_EQ(ai.toFlatIndex(ij), offset);
      ++offset;
    }
  }
}

template <int DIM>
typename std::enable_if<DIM == 3>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::ArrayIndexer<axom::IndexType, DIM>& ai)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> ijk;
  auto& l = ijk[fastestDirs[2]];  // Slowest dir
  auto& m = ijk[fastestDirs[1]];
  auto& n = ijk[fastestDirs[0]];  // Fastest dir
  axom::IndexType offset = 0;
  for(l = 0; l < lengths[fastestDirs[2]]; ++l)
  {
    for(m = 0; m < lengths[fastestDirs[1]]; ++m)
    {
      for(n = 0; n < lengths[fastestDirs[0]]; ++n)
      {
        SLIC_INFO(axom::fmt::format("offset {}   ijk {}", offset, ijk));
        EXPECT_EQ(ai.toMultiIndex(offset), ijk);
        EXPECT_EQ(ai.toFlatIndex(ijk), offset);
        ++offset;
      }
    }
  }
}

template <int DIM>
void check_arbitrary_strides(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs)
{
  // fastestDirs should be a permutation.
  SLIC_INFO(axom::fmt::format("Testing lengths {} with fastestDirs {}",
                              lengths,
                              fastestDirs));

  axom::StackArray<axom::IndexType, DIM> strides;
  axom::IndexType currentStride = 1;
  for(int d = 0; d < DIM; ++d)
  {
    axom::IndexType currentDir = fastestDirs[d];
    strides[currentDir] = currentStride;
    currentStride *= lengths[currentDir];
  }

  axom::ArrayIndexer<axom::IndexType, DIM> ai(strides);

  const auto slowestDirs = ai.slowestDirs();
  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_EQ(slowestDirs[d], fastestDirs[DIM - d - 1]);
  }

  check_arbitrary_strides_nested_loops(lengths, fastestDirs, ai);
}

// Test arbitrary strides.
TEST(quest_array_indexer, quest_arbitrary_strides)
{
  {
    constexpr int DIM = 1;
    // ArrayIndexer is overkill for 1D, but think of it as a smoke test.
    axom::StackArray<axom::IndexType, DIM> lengths {3};

    axom::StackArray<axom::IndexType, DIM> fastestDirs {0};
    check_arbitrary_strides(lengths, fastestDirs);
  }
  {
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};

    // Row-major ordering.
    axom::StackArray<axom::IndexType, DIM> rowFastestDirs {0, 1};
    check_arbitrary_strides(lengths, rowFastestDirs);

    // Col-major ordering.
    axom::StackArray<axom::IndexType, DIM> colFastestDirs {1, 0};
    check_arbitrary_strides(lengths, colFastestDirs);

    // A few random orderings.
    axom::StackArray<axom::IndexType, DIM> rand1FastestDirs {0, 1};
    check_arbitrary_strides(lengths, rand1FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand2FastestDirs {1, 0};
    check_arbitrary_strides(lengths, rand2FastestDirs);
  }
  {
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};

    // Row-major ordering.
    axom::StackArray<axom::IndexType, DIM> rowFastestDirs {0, 1, 2};
    check_arbitrary_strides(lengths, rowFastestDirs);

    // Col-major ordering.
    axom::StackArray<axom::IndexType, DIM> colFastestDirs {2, 1, 0};
    check_arbitrary_strides(lengths, colFastestDirs);

    // A few random orderings.
    axom::StackArray<axom::IndexType, DIM> rand1FastestDirs {0, 2, 1};
    check_arbitrary_strides(lengths, rand1FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand2FastestDirs {2, 0, 1};
    check_arbitrary_strides(lengths, rand2FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand3FastestDirs {1, 0, 2};
    check_arbitrary_strides(lengths, rand3FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand4FastestDirs {1, 2, 0};
    check_arbitrary_strides(lengths, rand4FastestDirs);
  }
}

// Test column-major element offsets with Array's.
TEST(quest_array_indexer, quest_array_match)
{
  // No test for 1D.  Array provides no non-trivial interface to test.

  {
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::Array<int, DIM> a(lengths);
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'c');
    for(axom::IndexType i = 0; i < lengths[0]; ++i)
    {
      for(axom::IndexType j = 0; j < lengths[1]; ++j)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(&a(i, j), &a.flatIndex(ai.toFlatIndex(indices)));
      }
    }
  }
  {
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};
    axom::Array<int, DIM> a(lengths);
    axom::ArrayIndexer<axom::IndexType, DIM> ai(lengths, 'c');
    for(axom::IndexType i = 0; i < lengths[0]; ++i)
    {
      for(axom::IndexType j = 0; j < lengths[1]; ++j)
      {
        for(axom::IndexType k = 0; j < lengths[2]; ++j)
        {
          axom::StackArray<axom::IndexType, DIM> indices {i, j, k};
          EXPECT_EQ(&a(i, j, k), &a.flatIndex(ai.toFlatIndex(indices)));
        }
      }
    }
  }
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  return result;
}
