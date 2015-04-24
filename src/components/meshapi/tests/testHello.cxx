#include "gtest/gtest.h"

#include "meshapi/Hello.hpp"
//------------------------------------------------------------------------------

TEST(gtest_meshapi_hello,basic_assert_example)
{
    EXPECT_TRUE( true );
}

TEST(gtest_meshapi_hello,test_calling_hello)
{
    asctoolkit::meshapi::Hello h;
    h.print();
    EXPECT_TRUE( true );
}
