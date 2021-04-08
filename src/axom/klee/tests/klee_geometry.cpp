// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"

#include "axom/klee/tests/KleeTestUtils.hpp"

#include "gtest/gtest.h"

#include <memory>

namespace axom
{
namespace klee
{
using test::MockOperator;
using ::testing::Return;

TEST(GeometryTest, dimensions_noOperators)
{
  TransformableGeometryProperties startProperties {Dimensions::Three,
                                                   LengthUnit::mils};
  Geometry geometry {startProperties, "test format", "test path", nullptr};
  EXPECT_EQ(startProperties, geometry.getStartProperties());
  EXPECT_EQ(startProperties, geometry.getEndProperties());
}

TEST(GeometryTest, dimensions_dimensionPreservingOperator)
{
  TransformableGeometryProperties startProperties {Dimensions::Two,
                                                   LengthUnit::mils};
  TransformableGeometryProperties endProperties {Dimensions::Three,
                                                 LengthUnit::cm};
  auto mockOperator = std::make_shared<MockOperator>(startProperties);
  Geometry geometry {startProperties, "test format", "test path", mockOperator};

  ON_CALL(*mockOperator, getEndProperties()).WillByDefault(Return(endProperties));
  EXPECT_CALL(*mockOperator, getEndProperties());

  EXPECT_EQ(startProperties, geometry.getStartProperties());
  EXPECT_EQ(endProperties, geometry.getEndProperties());
}

}  // namespace klee
}  // namespace axom
