// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IOUtil.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "conduit.hpp"

#include "KleeMatchers.hpp"

namespace axom { namespace klee { namespace internal {

using test::AlmostEqArray;
using ::testing::ElementsAre;

conduit::Node parseValue(const std::string &value) {
    conduit::Node node;
    std::string input = "value: ";
    input += value;
    node.parse(input, "yaml");
    return node["value"];
}

TEST(io_util, toDouble){
    EXPECT_DOUBLE_EQ(123.456, toDouble(parseValue("123.456")));
    EXPECT_DOUBLE_EQ(123, toDouble(parseValue("123")));
    EXPECT_THROW(toDouble(parseValue("abc")), std::invalid_argument);
    EXPECT_THROW(toDouble(parseValue("[1, 2]")), std::invalid_argument);
}

TEST(io_util, toDoubleVector){
    EXPECT_THAT(toDoubleVector(parseValue("[1.2, 3.4]"), 2),
            ElementsAre(1.2, 3.4));
    EXPECT_THAT(toDoubleVector(parseValue("[1, 2]"), 2),
            ElementsAre(1.0, 2.0));
    EXPECT_THROW(toDoubleVector(parseValue("[1, 2]"), 3),
            std::invalid_argument) << "Wrong length";
    EXPECT_THROW(toDoubleVector(parseValue("[a, b]"), 2),
            std::invalid_argument) << "Wrong type";
}

}}}
