#include "gtest/gtest.h"

#include "sidre/c-datastore.h"

//------------------------------------------------------------------------------

TEST(C_sidre_smoke,create_datastore)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datastore_delete(ds);
    EXPECT_TRUE( true );
}

//------------------------------------------------------------------------------
