#include "gtest/gtest.h"

// do we need one include that brings all of these in?
#include "datastore/datastore.h"

using DataStoreNS::DataStore;
using DataStoreNS::DataBuffer;

using namespace conduit;

//------------------------------------------------------------------------------

TEST(datastore_buffer,create_buffers)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff_0 = ds->CreateBuffer();
    DataBuffer *dbuff_1 = ds->CreateBuffer();
    
    EXPECT_EQ(dbuff_0->GetUID(),0);
    EXPECT_EQ(dbuff_1->GetUID(),1);
    ds->DestroyBuffer(0);
    
    DataBuffer *dbuff_3 = ds->CreateBuffer();
    EXPECT_EQ(dbuff_3->GetUID(),0);
    delete ds;
}


TEST(datastore_buffer,alloc_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->CreateBuffer();

    dbuff->SetDescriptor(DataType::Arrays::uint32(10));
    dbuff->Allocate();
    dbuff->ApplyDescriptor();
    
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->GetDescriptor().total_bytes());

    delete ds;
    
}


TEST(datastore_buffer,init_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->CreateBuffer();

    dbuff->Init(DataType::Arrays::uint32(10));
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->GetDescriptor().total_bytes());
    delete ds;
    
}

