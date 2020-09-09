// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include <stdexcept>

#include "gtest/gtest.h"

namespace axom { namespace klee { namespace {

TEST(ShapeSetTest, resolvePath_noPathSet) {
    ShapeSet shapeSet;
    EXPECT_THROW(shapeSet.resolvePath("anyPath"), std::logic_error);
}

TEST(ShapeSetTest, resolvePath_startWithSimpleFileName) {
    ShapeSet shapeSet;
    shapeSet.setPath("file.yaml");
    EXPECT_EQ("newPath.txt", shapeSet.resolvePath("newPath.txt"));
    EXPECT_EQ("d1/d2/newPath.txt", shapeSet.resolvePath("d1/d2/newPath.txt"));
    EXPECT_EQ("/abs/path/newPath.txt",
            shapeSet.resolvePath("/abs/path/newPath.txt"));
}

TEST(ShapeSetTest, resolvePath_startWithRelativeFileName) {
    ShapeSet shapeSet;
    shapeSet.setPath("path/to/file.yaml");
    EXPECT_EQ("path/to/newPath.txt", shapeSet.resolvePath("newPath.txt"));
    EXPECT_EQ("path/to/d1/d2/newPath.txt",
            shapeSet.resolvePath("d1/d2/newPath.txt"));
    EXPECT_EQ("/abs/path/newPath.txt",
            shapeSet.resolvePath("/abs/path/newPath.txt"));
}

TEST(ShapeSetTest, resolvePath_startWithAbsoluteFileName) {
    ShapeSet shapeSet;
    shapeSet.setPath("/path/to/file.yaml");
    EXPECT_EQ("/path/to/newPath.txt", shapeSet.resolvePath("newPath.txt"));
    EXPECT_EQ("/path/to/d1/d2/newPath.txt",
            shapeSet.resolvePath("d1/d2/newPath.txt"));
    EXPECT_EQ("/other/abs/path/newPath.txt",
            shapeSet.resolvePath("/other/abs/path/newPath.txt"));
}

}}}
