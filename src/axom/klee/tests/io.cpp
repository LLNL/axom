// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/io.hpp"

#include <sstream>

#include "gtest/gtest.h"

namespace axom { namespace klee { namespace {

ShapeSet readShapeSetFromString(std::string const &input) {
    std::istringstream istream(input);
    return readShapeSet(istream);
}

TEST(IOTest, readShapeSet_noShapes) {
    auto shapeSet = readShapeSetFromString("shapes: []");
    EXPECT_TRUE(shapeSet.getShapes().empty());
}

TEST(IOTest, readShapeSet_shapeWithNoReplacementLists) {
    auto shapeSet = readShapeSetFromString(R"(
shapes:
  - name: wheel
    material: steel
    geometry:
      format: test_format
      path: path/to/file.format
)");

    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    EXPECT_TRUE(shape.replaces("mat1"));
    EXPECT_TRUE(shape.replaces("mat2"));
    EXPECT_EQ("wheel", shape.getName());
    EXPECT_EQ("steel", shape.getMaterial());
    auto &geometry = shape.getGeometry();
    EXPECT_EQ("test_format", geometry.getFormat());
    EXPECT_EQ("path/to/file.format", geometry.getPath());
}

TEST(IOTest, readShapeSet_shapeWithReplacesList) {
    auto shapeSet = readShapeSetFromString(R"(
shapes:
  - name: wheel
    material: steel
    replaces: [mat1, mat2]
    geometry:
      format: test_format
      path: path/to/file.format
)");

    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    EXPECT_TRUE(shape.replaces("mat1"));
    EXPECT_TRUE(shape.replaces("mat2"));
    EXPECT_FALSE(shape.replaces("material_not_in_list"));
}

TEST(IOTest, readShapeSet_shapeWithDoesNotReplaceList) {
    auto shapeSet = readShapeSetFromString(R"(
shapes:
  - name: wheel
    material: steel
    does_not_replace: [mat1, mat2]
    geometry:
      format: test_format
      path: path/to/file.format
)");

    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    EXPECT_FALSE(shape.replaces("mat1"));
    EXPECT_FALSE(shape.replaces("mat2"));
    EXPECT_TRUE(shape.replaces("material_not_in_list"));
}

}}}
