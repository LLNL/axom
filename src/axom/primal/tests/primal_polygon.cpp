// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/primal/geometry/Polygon.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic/interface/slic.hpp"

using namespace axom;

typedef primal::Polygon<double, 3> PolygonType;

//------------------------------------------------------------------------------
TEST(primal_polygon, polygon_empty)
{
  PolygonType poly;
  EXPECT_FALSE(poly.isValid());
}

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
