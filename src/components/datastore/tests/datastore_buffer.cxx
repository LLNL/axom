#include "gtest/gtest.h"

// do we need one include that brings all of these in?
#include "datastore/DataStore.hpp"
#include "datastore/DataBuffer.hpp"

using DataStoreNS::DataStore;
using DataStoreNS::DataBuffer;

//------------------------------------------------------------------------------

TEST(datastore_buffer,create_buffers)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff_0 = ds->CreateDataBuffer();
    DataBuffer *dbuff_1 = ds->CreateDataBuffer();
    
    EXPECT_EQ(dbuff_0->GetUID(),0);
    EXPECT_EQ(dbuff_1->GetUID(),1);
    ds->DeleteDataBuffer(0);
    
    DataBuffer *dbuff_3 = ds->CreateDataBuffer();
    EXPECT_EQ(dbuff_3->GetUID(),0);
    delete ds;
}
