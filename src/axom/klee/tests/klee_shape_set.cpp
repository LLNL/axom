// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include <stdexcept>

#include "gtest/gtest.h"

namespace axom
{
namespace klee
{
namespace
{
TEST(ShapeSetTest, resolvePath_noPathSet)
{
  ShapeSet shapeSet;
  EXPECT_THROW(shapeSet.resolvePath("anyPath"), std::logic_error);
}

TEST(ShapeSetTest, resolvePath_startWithSimpleFileName)
{
  ShapeSet shapeSet;
  shapeSet.setPath("file.yaml");
  EXPECT_EQ("newPath.txt", shapeSet.resolvePath("newPath.txt"));
  EXPECT_EQ("d1/d2/newPath.txt", shapeSet.resolvePath("d1/d2/newPath.txt"));
  EXPECT_EQ("/abs/path/newPath.txt",
            shapeSet.resolvePath("/abs/path/newPath.txt"));
}

TEST(ShapeSetTest, resolvePath_startWithRelativeFileName)
{
  ShapeSet shapeSet;
  shapeSet.setPath("path/to/file.yaml");
  EXPECT_EQ("path/to/newPath.txt", shapeSet.resolvePath("newPath.txt"));
  EXPECT_EQ("path/to/d1/d2/newPath.txt",
            shapeSet.resolvePath("d1/d2/newPath.txt"));
  EXPECT_EQ("/abs/path/newPath.txt",
            shapeSet.resolvePath("/abs/path/newPath.txt"));
}

TEST(ShapeSetTest, resolvePath_startWithAbsoluteFileName)
{
  ShapeSet shapeSet;
  shapeSet.setPath("/path/to/file.yaml");
  EXPECT_EQ("/path/to/newPath.txt", shapeSet.resolvePath("newPath.txt"));
  EXPECT_EQ("/path/to/d1/d2/newPath.txt",
            shapeSet.resolvePath("d1/d2/newPath.txt"));
  EXPECT_EQ("/other/abs/path/newPath.txt",
            shapeSet.resolvePath("/other/abs/path/newPath.txt"));
}

TEST(ShapeSetTest, dimensions_getAndSet)
{
  ShapeSet shapeSet;

  // Can't query the ShapeSet dimensions until calling setDimensions
  EXPECT_THROW(shapeSet.getDimensions(), std::logic_error);

  {
    shapeSet.setDimensions(Dimensions::Two);

    EXPECT_EQ(Dimensions::Two, shapeSet.getDimensions());
    EXPECT_NE(Dimensions::Three, shapeSet.getDimensions());
  }

  {
    shapeSet.setDimensions(Dimensions::Three);

    EXPECT_NE(Dimensions::Two, shapeSet.getDimensions());
    EXPECT_EQ(Dimensions::Three, shapeSet.getDimensions());
  }
}

}  // namespace
}  // namespace klee
}  // namespace axom
