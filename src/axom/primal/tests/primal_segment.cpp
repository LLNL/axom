// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"

#include "fmt/fmt.hpp"

#include <cmath>

namespace primal = axom::primal;

namespace
{
const double EPS = 1e-12;
}

TEST(primal_segment, basic_2d)
{
  const int DIM = 2;
  using CoordType = double;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using PointType = typename SegmentType::PointType;

  PointType a {1.1, 2.2};
  PointType b {3.3, 4.4};
  SegmentType s(a, b);

  SLIC_INFO("Testing segment " << s);

  EXPECT_EQ(a, s.source());
  EXPECT_EQ(b, s.target());
  EXPECT_NE(s.source(), s.target());
}

TEST(primal_segment, basic_3d)
{
  const int DIM = 3;
  using CoordType = double;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using PointType = typename SegmentType::PointType;

  PointType a {1.1, 2.1, 3.1};
  PointType b {1.8, 2.8, 3.8};
  SegmentType s(a, b);

  SLIC_INFO("Testing segment " << s);

  EXPECT_EQ(a, s.source());
  EXPECT_EQ(b, s.target());
  EXPECT_NE(s.source(), s.target());
}

TEST(primal_segment, length_2d)
{
  const int DIM = 2;
  using CoordType = double;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using PointType = typename SegmentType::PointType;

  SegmentType s0010(PointType {0, 0}, PointType {1, 0});
  EXPECT_NEAR(1., s0010.length(), EPS);

  // compare lengths of several basic segments
  {
    SegmentType s0001(PointType {0, 0}, PointType {0, 1});
    SegmentType s1000(PointType {1, 0}, PointType {0, 0});
    SegmentType s0100(PointType {0, 1}, PointType {0, 0});
    EXPECT_NEAR(s0010.length(), s0001.length(), EPS);
    EXPECT_NEAR(s0010.length(), s1000.length(), EPS);
    EXPECT_NEAR(s0010.length(), s0100.length(), EPS);
  }

  // Tests segment of length 0
  {
    SegmentType s1111(PointType {1, 1}, PointType {1, 1});
    EXPECT_NEAR(0., s1111.length(), EPS);
  }
}

TEST(primal_segment, interpolation_using_at_2d)
{
  const int DIM = 2;
  using CoordType = double;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using PointType = typename SegmentType::PointType;

  // Tests interpolation on zero-length segments
  {
    PointType a {1, 2};
    PointType b {2, 4};
    SegmentType s(a, b);

    EXPECT_EQ(a, s.at(0));
    EXPECT_EQ(b, s.at(1));

    PointType at_half {1.5, 3};
    EXPECT_NEAR(at_half[0], s.at(0.5)[0], EPS);
    EXPECT_NEAR(at_half[1], s.at(0.5)[1], EPS);

    PointType at_two {3, 6};
    EXPECT_NEAR(at_two[0], s.at(2.0)[0], EPS);
    EXPECT_NEAR(at_two[1], s.at(2.0)[1], EPS);
  }

  // Tests interpolation on zero-length segments
  {
    PointType a {1, 1};
    SegmentType s(a, a);
    for(int i = 0; i < 2; ++i)
    {
      EXPECT_NEAR(a[i], s.at(0.0)[i], EPS);
      EXPECT_NEAR(a[i], s.at(0.5)[i], EPS);
      EXPECT_NEAR(a[i], s.at(1.0)[i], EPS);
    }
  }
}

TEST(primal_segment, normal_2d)
{
  const int DIM = 2;
  using CoordType = double;
  using SegmentType = primal::Segment<CoordType, DIM>;
  using PointType = typename SegmentType::PointType;
  using VectorType = typename SegmentType::VectorType;

  {
    SCOPED_TRACE("typical normal");

    PointType a {1, 1};
    PointType b {4, 3};
    SegmentType s(a, b);

    VectorType normal {-2, 3};
    EXPECT_NEAR(normal[0], s.normal<2>()[0], EPS);
    EXPECT_NEAR(normal[1], s.normal<2>()[1], EPS);

    SegmentType reverse(b, a);

    VectorType reverse_normal {2, -3};
    EXPECT_NEAR(reverse_normal[0], reverse.normal<2>()[0], EPS);
    EXPECT_NEAR(reverse_normal[1], reverse.normal<2>()[1], EPS);
  }

  {
    SCOPED_TRACE("axis-aligned: +x");

    PointType a {0, 0};
    PointType b {1, 0};
    SegmentType s(a, b);

    VectorType normal {0, 1};
    EXPECT_NEAR(normal[0], s.normal<2>()[0], EPS);
    EXPECT_NEAR(normal[1], s.normal<2>()[1], EPS);
  }

  {
    SCOPED_TRACE("axis-aligned: -x");

    PointType a {0, 0};
    PointType b {-1, 0};
    SegmentType s(a, b);

    VectorType normal {0, -1};
    EXPECT_NEAR(normal[0], s.normal<2>()[0], EPS);
    EXPECT_NEAR(normal[1], s.normal<2>()[1], EPS);
  }

  {
    SCOPED_TRACE("axis-aligned: +y");

    PointType a {0, 0};
    PointType b {0, 1};
    SegmentType s(a, b);

    VectorType normal {1, 0};
    EXPECT_NEAR(normal[0], s.normal<2>()[0], EPS);
    EXPECT_NEAR(normal[1], s.normal<2>()[1], EPS);
  }

  {
    SCOPED_TRACE("axis-aligned: -y");

    PointType a {0, 0};
    PointType b {0, -1};
    SegmentType s(a, b);

    VectorType normal {-1, 0};
    EXPECT_NEAR(normal[0], s.normal<2>()[0], EPS);
    EXPECT_NEAR(normal[1], s.normal<2>()[1], EPS);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  int result = RUN_ALL_TESTS();
  return result;
}
