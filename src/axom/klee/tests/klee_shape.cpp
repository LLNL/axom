// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Shape.hpp"

#include <stdexcept>

#include "gtest/gtest.h"

namespace axom { namespace klee { namespace {

TEST(ShapeTest, replaces_no_lists_given) {
  Shape shape;
  EXPECT_TRUE(shape.replaces("some_material"));
}

TEST(ShapeTest, replaces_replacement_list_given) {
  Shape shape;
  shape.setMaterialsReplaced({"mat1", "mat2"});
  EXPECT_TRUE(shape.replaces("mat1"));
  EXPECT_TRUE(shape.replaces("mat2"));
  EXPECT_FALSE(shape.replaces("mat3"));
}

TEST(ShapeTest, replaces_non_replacement_list_given) {
  Shape shape;
  shape.setMaterialsNotReplaced({"mat1", "mat2"});
  EXPECT_FALSE(shape.replaces("mat1"));
  EXPECT_FALSE(shape.replaces("mat2"));
  EXPECT_TRUE(shape.replaces("mat3"));
}

TEST(ShapeTest, set_not_replaced_list_after_replaced_list) {
    Shape shape;
    shape.setMaterialsReplaced({"mat1", "mat2"});
    EXPECT_THROW(shape.setMaterialsNotReplaced({"mat1", "mat2"}),
            std::logic_error);
}

TEST(ShapeTest, set_replaced_list_after_not_replaced_list) {
    Shape shape;
    shape.setMaterialsNotReplaced({"mat1", "mat2"});
    EXPECT_THROW(shape.setMaterialsReplaced({"mat1", "mat2"}),
            std::logic_error);
}

}}}
