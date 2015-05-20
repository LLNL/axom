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

#include "sidre/sidre.hpp"


#include "conduit/conduit.h"

using asctoolkit::sidre::DataStore;

//------------------------------------------------------------------------------

TEST(sidre_smoke,create_datastore)
{
    DataStore *ds = new DataStore();
    delete ds;
    EXPECT_TRUE( true );
}

//------------------------------------------------------------------------------

TEST(sidre_smoke,conduit_in_sidre_smoke)
{
    // make sure we are linking with conduit ok. 
    conduit::Node n;
    n["field"] = 100;
    EXPECT_EQ(n["field"].to_index_t(),100u);
}

