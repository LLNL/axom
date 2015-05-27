/*
 * testSet.cxx
 *
 *  Created on: Apr 23, 2015
 *      Author: weiss27
 */

#include "gtest/gtest.h"

#include "meshapi/Set.hpp"
#include "meshapi/OrderedSet.hpp"

static unsigned int const NUM_ELEMS = 5;

TEST(gtest_meshapi_set,construct_set)
{
    asctoolkit::meshapi::Set* s = new asctoolkit::meshapi::OrderedSet(NUM_ELEMS);
    delete s;

    EXPECT_TRUE( true );
}
