// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"

#include "axom/klee/tests/KleeTestUtils.hpp"

#include "gtest/gtest.h"

#include <memory>

namespace axom { namespace klee {

using test::MockOperator;
using ::testing::Return;

TEST(GeometryTest, dimensions_noOperators) {
    Geometry geometry;
    geometry.setInitialDimensions(Dimensions::Three);
    EXPECT_EQ(Dimensions::Three, geometry.getInitialDimensions());
    EXPECT_EQ(Dimensions::Three, geometry.getDimensions());
}

TEST(GeometryTest, dimensions_dimensionPreservingOperator) {
    Geometry geometry;
    geometry.setInitialDimensions(Dimensions::Three);

    auto mockOperator = std::make_shared<MockOperator>();
    ON_CALL(*mockOperator, startDims()).WillByDefault(Return(Dimensions::Three));
    ON_CALL(*mockOperator, endDims()).WillByDefault(Return(Dimensions::Two));
    geometry.setGeometryOperator(mockOperator);
    EXPECT_CALL(*mockOperator, endDims());

    EXPECT_EQ(Dimensions::Three, geometry.getInitialDimensions());
    EXPECT_EQ(Dimensions::Two, geometry.getDimensions());
}

}}
