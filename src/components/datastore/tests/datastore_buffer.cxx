#include "gtest/gtest.h"

// do we need one include that brings all of these in?
#include "datastore/datastore.h"

using sidre::DataStore;
using sidre::DataBuffer;

using namespace conduit;

//------------------------------------------------------------------------------

TEST(datastore_buffer,create_buffers)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff_0 = ds->createBuffer();
    DataBuffer *dbuff_1 = ds->createBuffer();
    
    EXPECT_EQ(dbuff_0->GetUID(),0);
    EXPECT_EQ(dbuff_1->GetUID(),1);
    ds->destroyBuffer(0);
    
    DataBuffer *dbuff_3 = ds->createBuffer();
    EXPECT_EQ(dbuff_3->GetUID(),0);
    ds->print();
    delete ds;
}


TEST(datastore_buffer,alloc_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->Declare(DataType::uint32(10));
    dbuff->Allocate();
    
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->getDescriptor().total_bytes());
  
    ds->print();
    delete ds;
    
}


TEST(datastore_buffer,init_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->Allocate(DataType::uint32(10));
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->getDescriptor().total_bytes());

    ds->print();
    delete ds;
    
}

