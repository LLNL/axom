/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "sidre/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_smoke,create_datastore)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datastore_delete(ds);
  EXPECT_TRUE( true );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
   int result = 0;

   ::testing::InitGoogleTest(&argc, argv);

   UnitTestLogger logger;  // create & initialize test logger,
                       // finalized when exiting main scope

   result = RUN_ALL_TESTS();

   return result;
}
