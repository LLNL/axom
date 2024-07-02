// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"
#include "axom/mir/views/StridedStructuredIndexing.hpp"

//----------------------------------------------------------------------

TEST(mir_views_indexing, strided_structured_indexing_2d)
{
  using Indexing2D = axom::mir::views::StridedStructuredIndexing<int, 2>;
  using LogicalIndex = typename axom::mir::views::StridedStructuredIndexing<int, 2>::LogicalIndex;
  LogicalIndex dims{6, 5};
  LogicalIndex origin{2, 2};
  LogicalIndex stride{1, 6};
  Indexing2D indexing(dims, origin, stride);

  EXPECT_EQ(indexing.dimensions(), 2);
  EXPECT_EQ(indexing.size(), dims[0] * dims[1]);

  const LogicalIndex logical0_0{0, 0};
  const auto index0_0 = indexing.LogicalIndexToIndex(logical0_0);
  EXPECT_EQ(index0_0, 14);

  const LogicalIndex logical2_2{2, 2};
  const auto index2_2 = indexing.LogicalIndexToIndex(logical2_2);
  EXPECT_EQ(index2_2, 28);

  LogicalIndex logical = indexing.IndexToLogicalIndex(index2_2);
  EXPECT_TRUE(logical == logical2_2);

  for(int j = 0; j < 3; j++)
  {
     for(int i = 0; i < 4; i++)
     {
        LogicalIndex logical{i, j};
        const auto flat = indexing.LogicalIndexToIndex(logical);
        const auto logical2 = indexing.IndexToLogicalIndex(flat);
        EXPECT_EQ(logical, logical2);
     }
  } 
}

//----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
