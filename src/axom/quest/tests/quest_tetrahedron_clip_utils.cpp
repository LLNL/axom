// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/detail/clipping/TetrahedronClipUtils.hpp"

#include "gtest/gtest.h"

#include <algorithm>

namespace
{

using axom::quest::experimental::detail::clipTetByVertexValues;
using axom::quest::experimental::detail::subdivideTetByMidpoints;
using axom::quest::experimental::detail::TetVertexValues;

TEST(quest_tetrahedron_clip_utils, affine_volume_cases)
{
  constexpr double volume = 2.5;

  const TetVertexValues outside {{-1.0, -2.0, -3.0, -4.0}};
  EXPECT_DOUBLE_EQ(0.0, clipTetByVertexValues(outside, volume));

  const TetVertexValues inside {{1.0, 2.0, 3.0, 4.0}};
  EXPECT_DOUBLE_EQ(volume, clipTetByVertexValues(inside, volume));

  const TetVertexValues oneInside {{1.0, -1.0, -1.0, -1.0}};
  EXPECT_DOUBLE_EQ(volume / 8.0, clipTetByVertexValues(oneInside, volume));

  const TetVertexValues threeInside {{-1.0, 1.0, 1.0, 1.0}};
  EXPECT_DOUBLE_EQ(7.0 * volume / 8.0, clipTetByVertexValues(threeInside, volume));

  const TetVertexValues twoInside {{1.0, 1.0, -1.0, -1.0}};
  EXPECT_DOUBLE_EQ(volume / 2.0, clipTetByVertexValues(twoInside, volume));
}

TEST(quest_tetrahedron_clip_utils, two_inside_is_permutation_invariant)
{
  constexpr double volume = 3.0;
  constexpr double expectedFraction = 37.0 / 105.0;
  TetVertexValues values {{-4.0, -2.0, 1.0, 3.0}};
  do
  {
    EXPECT_NEAR(volume * expectedFraction, clipTetByVertexValues(values, volume), 1e-14);
  } while(std::next_permutation(values.begin(), values.end()));
}

TEST(quest_tetrahedron_clip_utils, midpoint_subdivision_supports_aliased_input)
{
  using PointType = axom::quest::experimental::MeshClipperStrategy::Point3DType;
  using TetType = axom::quest::experimental::MeshClipperStrategy::TetrahedronType;

  const TetType parent(PointType {0.0, 0.0, 0.0},
                       PointType {1.0, 0.0, 0.0},
                       PointType {0.0, 1.0, 0.0},
                       PointType {0.0, 0.0, 1.0});
  axom::StackArray<TetType, 8> expected;
  subdivideTetByMidpoints(parent, expected.data());

  axom::StackArray<TetType, 8> aliased;
  aliased[0] = parent;
  subdivideTetByMidpoints(aliased[0], aliased.data());

  double childVolumeSum = 0.0;
  for(int child = 0; child < 8; ++child)
  {
    EXPECT_GT(aliased[child].signedVolume(), 0.0);
    childVolumeSum += aliased[child].volume();
    for(int vertex = 0; vertex < 4; ++vertex)
    {
      for(int dimension = 0; dimension < 3; ++dimension)
      {
        EXPECT_DOUBLE_EQ(expected[child][vertex][dimension], aliased[child][vertex][dimension]);
      }
    }
  }
  EXPECT_DOUBLE_EQ(parent.volume(), childVolumeSum);
}

}  // namespace
