// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include <string>

#include "gtest/gtest.h"

#include "axom/sina/include/CppBridge.hpp"

namespace axom
{
namespace sina
{
namespace internal
{
namespace testing
{
namespace
{


TEST(CppBridge, make_unique_noParams) {
    std::unique_ptr<std::string> ptr = make_unique<std::string>();
    EXPECT_TRUE(ptr->empty());
}

TEST(CppBridge, make_unique_withParam) {
    std::unique_ptr<std::string> ptr = make_unique<std::string>("foo");
    EXPECT_EQ("foo", *ptr);
}

}  // end nameless namespace
}  // end testing namespace
}  // end internal namespace
}  // end sina namespace
}  // end axom namespace
