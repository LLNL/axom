/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#define CATCH_CONFIG_MAIN
#include "catch.hpp"


//------------------------------------------------------------------------------

TEST_CASE("basic_test", "This basic catch test should pass")
{

  int one = 1;
  REQUIRE_FALSE( one == 2 );

}
