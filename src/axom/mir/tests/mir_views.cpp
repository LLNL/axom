// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"


TEST(mir_views, shape2conduitName)
{
  EXPECT_EQ(axom::mir::views::LineShape<int>::name, "line");
  EXPECT_EQ(axom::mir::views::LineShape<long>::name, "line");

  EXPECT_EQ(axom::mir::views::TriShape<int>::name, "tri");
  EXPECT_EQ(axom::mir::views::TriShape<long>::name, "tri");

  EXPECT_EQ(axom::mir::views::QuadShape<int>::name, "quad");
  EXPECT_EQ(axom::mir::views::QuadShape<long>::name, "quad");

  EXPECT_EQ(axom::mir::views::TetShape<int>::name, "tet");
  EXPECT_EQ(axom::mir::views::TetShape<long>::name, "tet");

  EXPECT_EQ(axom::mir::views::PyramidShape<int>::name, "pyramid");
  EXPECT_EQ(axom::mir::views::PyramidShape<long>::name, "pyramid");

  EXPECT_EQ(axom::mir::views::WedgeShape<int>::name, "wedge");
  EXPECT_EQ(axom::mir::views::WedgeShape<long>::name, "wedge");

  EXPECT_EQ(axom::mir::views::HexShape<int>::name, "hex");
  EXPECT_EQ(axom::mir::views::HexShape<long>::name, "hex");
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
