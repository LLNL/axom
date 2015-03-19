#include "gtest/gtest.h"

#include "datastore/DataStore.hpp"

using DataStoreNS::DataStore;

//------------------------------------------------------------------------------

TEST(datastore_smoke,create_datastore)
{
    DataStore *ds = new DataStore();
    delete ds;
    EXPECT_TRUE( true );
}
