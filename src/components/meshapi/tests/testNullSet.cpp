/*
 * \file testNullSet.cpp
 *
 * Unit tests for the NullSet class
 */

#include "gtest/gtest.h"

#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

TEST(gtest_meshapi_set,construct_nullset)
{
    asctoolkit::meshapi::Set* s = new asctoolkit::meshapi::NullSet();
    delete s;

    EXPECT_TRUE( true );
}

TEST(gtest_meshapi_set,subscript_fails_nullset)
{
    std::cout<<"\n****** Testing subscipt access on NullSet -- code is expected to assert and die." << std::endl;

    typedef asctoolkit::meshapi::Set::SizeType SizeType;
    asctoolkit::meshapi::NullSet n;

    EXPECT_EQ(n.size(), SizeType()) <<"size of null set is defined to be zero";

    // add this line to avoid a warning in the output about thread safety
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(n[0],"") << "subscript operator on null set asserts";
}
