#define CATCH_CONFIG_MAIN
#include "catch.hpp"


TEST_CASE( "This should fail") {

  int one = 1;
  REQUIRE( one == 2 );

}

//------------------------------------------------------------------------------

TEST_CASE( "This should pass") {

  int one = 1;
  REQUIRE_FALSE( one == 2 );

}
