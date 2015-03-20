#include "gtest/gtest.h"

#include "datastore/datastore.h"

#include "conduit/conduit.h"

using DataStoreNS::DataStore;

//------------------------------------------------------------------------------

TEST(datastore_smoke,create_datastore)
{
    DataStore *ds = new DataStore();
    delete ds;
    EXPECT_TRUE( true );
}

//------------------------------------------------------------------------------

TEST(datastore_smoke,conduit_in_datastore_smoke)
{
    // make sure we are linking with conduit ok. 
    conduit::Node n;
    n["field"] = 100;
    EXPECT_EQ(n["field"].to_index_t(),100);
}

