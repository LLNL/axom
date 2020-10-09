// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IO.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "gtest/gtest.h"

#include "axom/klee/GeometryOperators.hpp"

namespace axom { namespace klee { namespace {

ShapeSet readShapeSetFromString(const std::string &input) {
    std::istringstream istream(input);
    return readShapeSet(istream);
}

TEST(IOTest, readShapeSet_noShapes) {
    auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
        shapes: [])");
    EXPECT_TRUE(shapeSet.getShapes().empty());
}

TEST(IOTest, readShapeSet_invalidDimensions) {
    EXPECT_THROW(readShapeSetFromString(R"(
        dimensions: 5
        shapes: [])"), std::invalid_argument);
}

TEST(IOTest, readShapeSet_shapeWithNoReplacementLists) {
    auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
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
    EXPECT_FALSE(geometry.getGeometryOperator());
}

TEST(IOTest, readShapeSet_shapeWithReplacesList) {
    auto shapeSet = readShapeSetFromString(R"(
        dimensions: 2
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
        dimensions: 2
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

TEST(IOTest, readShapeSet_file) {
    std::string fileName = "testFile.yaml";

    std::string fileContents = R"(
    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: test_format
          path: path/to/file.format)";
    std::ofstream fout{fileName};
    fout << fileContents;
    fout.close();

    auto shapeSet = readShapeSet(fileName);
    EXPECT_EQ(1u, shapeSet.getShapes().size());
    EXPECT_EQ("testFile.yaml", shapeSet.getPath());
}

TEST(IOTest, readShapeSet_shapeWithReplacesAndDoesNotReplaceLists) {
    EXPECT_THROW(readShapeSetFromString(R"(
      dimensions: 2
      shapes:
        - name: wheel
          material: steel
          replaces: [mat1, mat2]
          does_not_replace: [mat1, mat2]
          geometry:
            format: test_format
            path: path/to/file.format
    )"), std::invalid_argument);
}

TEST(IOTest, readShapeSet_geometryOperators) {
    auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            operators:
              - rotate: 90
              - translate: [10, 20]
    )");
    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    auto &geometryOperator = shape.getGeometry().getGeometryOperator();
    ASSERT_TRUE(geometryOperator);
    auto composite = std::dynamic_pointer_cast<const CompositeOperator>(
            geometryOperator);
    ASSERT_TRUE(composite);
    EXPECT_EQ(2u, composite->getOperators().size());
}

TEST(IOTest, readShapeSet_differentDimensions) {
    auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            initial_dimensions: 3
            operators:
              - slice:
                 x: 10
    )");
    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    auto &geometry = shape.getGeometry();
    EXPECT_EQ(Dimensions::Three, geometry.getInitialDimensions());
    EXPECT_EQ(Dimensions::Two, geometry.getDimensions());
    auto &geometryOperator = geometry.getGeometryOperator();
    ASSERT_TRUE(geometryOperator);
    auto composite = std::dynamic_pointer_cast<const CompositeOperator>(
            geometryOperator);
    ASSERT_TRUE(composite);
    EXPECT_EQ(1u, composite->getOperators().size());
    auto slice = std::dynamic_pointer_cast<const SliceOperator>(
            composite->getOperators()[0]);
    EXPECT_TRUE(slice);
}

TEST(IOTest, readShapeSet_namedGeometryOperators) {
    auto shapeSet = readShapeSetFromString(R"(
      dimensions: 2

      shapes:
        - name: wheel
          material: steel
          geometry:
            format: test_format
            path: path/to/file.format
            operators:
              - ref: my_operation

      named_operators:
        - name: my_operation
          value:
            - rotate: 90
            - translate: [10, 20]
    )");
    auto &shapes = shapeSet.getShapes();
    ASSERT_EQ(1u, shapes.size());
    auto &shape = shapes[0];
    auto &geometryOperator = shape.getGeometry().getGeometryOperator();
    ASSERT_TRUE(geometryOperator);
    auto composite = std::dynamic_pointer_cast<const CompositeOperator>(
            geometryOperator);
    ASSERT_TRUE(composite);
    EXPECT_EQ(1u, composite->getOperators().size());
}

}}}
