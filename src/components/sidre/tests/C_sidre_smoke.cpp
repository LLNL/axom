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

//------------------------------------------------------------------------------

TEST(sidre_smoke,valid_invalid)
{
  ATK_datastore * ds = ATK_datastore_new();

  ATK_IndexType idx = 3;
  EXPECT_TRUE(idx != ATK_InvalidIndex);

  const char * name = "foo";
  EXPECT_TRUE(ATK_is_name_valid(name)); 

  ATK_datagroup * root = ATK_datastore_get_root(ds);

  const char * gp_name = ATK_datagroup_get_group_name(root, idx);
  EXPECT_TRUE(gp_name == NULL);
  EXPECT_TRUE(gp_name == ATK_InvalidName);
  EXPECT_FALSE(ATK_is_name_valid(gp_name));
  EXPECT_TRUE(ATK_datagroup_get_group_index(root, name) == ATK_InvalidIndex);

  ATK_datastore_delete(ds);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
