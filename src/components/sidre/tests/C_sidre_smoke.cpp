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
