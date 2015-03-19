#include "gtest/gtest.h"

// do we need one include that brings all of these in?
#include "datastore/DataStore.hpp"
#include "datastore/DataBuffer.hpp"

using DataStoreNS::DataStore;
using DataStoreNS::DataBuffer;

//------------------------------------------------------------------------------

TEST(datastore_smoke,create_datastore)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->CreateDataBuffer();
    EXPECT_EQ(dbuff->GetUID(),1);
    

//    delete ds;
}
