#include "gtest/gtest.h"

// do we need one include that brings all of these in?
#include "datastore/DataStore.hpp"
#include "datastore/DataBufferhpp"

using DataStoreNS::DataStore;
using DataStoreNS::DataBuffer;

//------------------------------------------------------------------------------

TEST(datastore_smoke,create_datastore)
{
    DataStore *ds = new DataStore();

    EXPECT_TRUE( true );
    
    
    DataBuffer *dbuff = ds->CreateDataBuffer();
    
    delete ds;
}
