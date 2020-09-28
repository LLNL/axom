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
    geometry.setInitialDimensions(3);
    EXPECT_EQ(3, geometry.getInitialDimensions());
    EXPECT_EQ(3, geometry.getDimensions());
}

TEST(GeometryTest, dimensions_dimensionPreservingOperator) {
    Geometry geometry;
    geometry.setInitialDimensions(3);

    auto mockOperator = std::make_shared<MockOperator>();
    ON_CALL(*mockOperator, startDims()).WillByDefault(Return(3));
    ON_CALL(*mockOperator, endDims()).WillByDefault(Return(2));
    geometry.setGeometryOperator(mockOperator);
    EXPECT_CALL(*mockOperator, endDims());

    EXPECT_EQ(3, geometry.getInitialDimensions());
    EXPECT_EQ(2, geometry.getDimensions());
}

}}
