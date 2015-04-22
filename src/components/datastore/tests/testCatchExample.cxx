#define CATCH_CONFIG_MAIN
#include "catch.hpp"


//------------------------------------------------------------------------------

TEST_CASE( "This basic catch test should pass") {

  int one = 1;
  REQUIRE_FALSE( one == 2 );

}
