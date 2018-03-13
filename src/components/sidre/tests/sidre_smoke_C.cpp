/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "sidre/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_smoke,create_datastore)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_datastore_delete(ds);
  EXPECT_TRUE( true );
}

//------------------------------------------------------------------------------

TEST(sidre_smoke,valid_invalid)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();

  SIDRE_IndexType idx = 3;
  EXPECT_TRUE(idx != SIDRE_InvalidIndex);

  const char* name = "foo";
  EXPECT_TRUE(SIDRE_name_is_valid(name));

  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  const char* gp_name = SIDRE_group_get_group_name(root, idx);
  EXPECT_TRUE(gp_name == NULL);
  EXPECT_TRUE(gp_name == SIDRE_InvalidName);
  EXPECT_FALSE(SIDRE_name_is_valid(gp_name));
  EXPECT_TRUE(SIDRE_group_get_group_index(root,
                                          name) == SIDRE_InvalidIndex);

  SIDRE_datastore_delete(ds);
}
