#include "gtest/gtest.h"

#include "common/Utilities.hpp"

//------------------------------------------------------------------------------

// Putting a test here to get us started, should add other functions too.
// -- Aaron
TEST(utilities,atk_macros_pass)
{

    int value = 5;
    ATK_ASSERT( value == 5 );
    ATK_ASSERT_MSG( value == 5, "ATK_ASSERT passed, shouldn't see this");

    ATK_ASSERT_MSG( value != 5, "ATK_ASSERT should have failed check of " << value << " != 5" );
}
