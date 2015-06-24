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

#ifdef ATK_DEBUG
    // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode

    // add this line to avoid a warning in the output about thread safety
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH(
        ATK_ASSERT_MSG( value != 5, "ATK_ASSERT should have failed check of " << value << " != 5" )
        , ""
        ) << "Testing AKT_ASSERT feature.";
#endif
}
