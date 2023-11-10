// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
  axom::StackArray<axom::IndexType, 3> lengths{5,3,2};

  axom::ArrayIndexer<axom::IndexType, 3> colIndexer(lengths, 'c');
  axom::StackArray<std::uint16_t, 3> colSlowestDirs{0,1,2};
  axom::StackArray<axom::IndexType, 3> colStrides{6,2,1};
  for(int d=0; d<3; ++d)
  {
    EXPECT_EQ( colIndexer.slowestDirs()[d], colSlowestDirs[d] );
    EXPECT_EQ( colIndexer.strides()[d], colStrides[d] );
  }

  axom::ArrayIndexer<axom::IndexType, 3> rowIndexer(lengths, 'r');
  axom::StackArray<std::uint16_t, 3> rowSlowestDirs{2,1,0};
  axom::StackArray<axom::IndexType, 3> rowStrides{1,5,15};
  for(int d=0; d<3; ++d)
  {
    EXPECT_EQ( rowIndexer.slowestDirs()[d], rowSlowestDirs[d] );
    EXPECT_EQ( rowIndexer.strides()[d], rowStrides[d] );
  }
}

// Test column-major offsets.
TEST(quest_array_indexer, quest_col_major_offset)
{
  axom::StackArray<axom::IndexType, 3> lengths{5,3,2};
  axom::ArrayIndexer<axom::IndexType, 3> ai(lengths, 'c');
  axom::IndexType offset = 0;
  for(int i=0; i<lengths[0]; ++i)
  {
    for(int j=0; j<lengths[1]; ++j)
    {
      for(int k=0; k<lengths[2]; ++k)
      {
        axom::StackArray<axom::IndexType, 3> indices{i,j,k};
        EXPECT_EQ( ai.toFlatIndex(indices), offset );
        EXPECT_EQ( ai.toMultiIndex(offset), indices );
        ++offset;
      }
    }
  }
}

// Test row-major offsets.
TEST(quest_array_indexer, quest_row_major_offset)
{
  axom::StackArray<axom::IndexType, 3> lengths{5,3,2};
  axom::ArrayIndexer<axom::IndexType, 3> ai(lengths, 'r');
  axom::IndexType offset = 0;
  for(int k=0; k<lengths[2]; ++k)
  {
    for(int j=0; j<lengths[1]; ++j)
    {
      for(int i=0; i<lengths[0]; ++i)
      {
        axom::StackArray<axom::IndexType, 3> indices{i,j,k};
        EXPECT_EQ( ai.toFlatIndex(indices), offset );
        EXPECT_EQ( ai.toMultiIndex(offset), indices );
        ++offset;
      }
    }
  }
}

void check_arbitrary_strides(
  const axom::StackArray<axom::IndexType, 3>& lengths,
  const axom::StackArray<axom::IndexType, 3>& fastestDirs)
{
  // fastestDirs should be a permutation.
  SLIC_INFO(axom::fmt::format("Testing lengths {} with fastestDirs {}", lengths, fastestDirs));

  axom::StackArray<axom::IndexType, 3> strides;
  axom::IndexType currentStride = 1;
  for(int d=0; d<3; ++d)
  {
    axom::IndexType currentDir = fastestDirs[d];
    strides[currentDir] = currentStride;
    currentStride *= lengths[currentDir];
  }

  axom::ArrayIndexer<axom::IndexType, 3> ai(strides);

  const auto slowestDirs = ai.slowestDirs();
  for(int d=0; d<3; ++d)
  {
     EXPECT_EQ( slowestDirs[d], fastestDirs[3-d-1] );
  }

  axom::IndexType offset = 0;
  axom::StackArray<axom::IndexType, 3> ijk; // Natural indices.
  auto& l = ijk[fastestDirs[2]]; // Slowest dir
  auto& m = ijk[fastestDirs[1]];
  auto& n = ijk[fastestDirs[0]]; // Fastest dir
  for(l=0; l<lengths[fastestDirs[2]]; ++l)
  {
    for(m=0; m<lengths[fastestDirs[1]]; ++m)
    {
      for(n=0; n<lengths[fastestDirs[0]]; ++n)
      {
        SLIC_INFO(axom::fmt::format("offset {}   ijk {}", offset, ijk));
        EXPECT_EQ( ai.toMultiIndex(offset), ijk );
        EXPECT_EQ( ai.toFlatIndex(ijk), offset );
        ++offset;
      }
    }
  }
}

// Test arbitrary strides.
TEST(quest_array_indexer, quest_arbitrary_strides)
{
  axom::StackArray<axom::IndexType, 3> lengths{5,3,2};

  // Row-major ordering.
  axom::StackArray<axom::IndexType, 3> rowFastestDirs{0,1,2};
  check_arbitrary_strides(lengths, rowFastestDirs);

  // Col-major ordering.
  axom::StackArray<axom::IndexType, 3> colFastestDirs{2,1,0};
  check_arbitrary_strides(lengths, colFastestDirs);

  // A few random orderings.
  axom::StackArray<axom::IndexType, 3> rand1FastestDirs{0,2,1};
  check_arbitrary_strides(lengths, rand1FastestDirs);
  axom::StackArray<axom::IndexType, 3> rand2FastestDirs{2,0,1};
  check_arbitrary_strides(lengths, rand2FastestDirs);
  axom::StackArray<axom::IndexType, 3> rand3FastestDirs{1,0,2};
  check_arbitrary_strides(lengths, rand3FastestDirs);
  axom::StackArray<axom::IndexType, 3> rand4FastestDirs{1,2,0};
  check_arbitrary_strides(lengths, rand4FastestDirs);
}

// Test column-major element offsets with Array's.
TEST(quest_array_indexer, quest_array_match)
{
  axom::StackArray<axom::IndexType, 3> lengths{5,3,2};
  axom::Array<int, 3> a(lengths);
  axom::ArrayIndexer<axom::IndexType, 3> ai(lengths, 'c');
  for(axom::IndexType i=0; i<lengths[0]; ++i)
  {
    for(axom::IndexType j=0; j<lengths[1]; ++j)
    {
      for(axom::IndexType k=0; j<lengths[2]; ++j)
      {
        axom::StackArray<axom::IndexType, 3> indices{i,j,k};
        EXPECT_EQ( &a(i,j,k), &a.flatIndex(ai.toFlatIndex(indices)) );
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
