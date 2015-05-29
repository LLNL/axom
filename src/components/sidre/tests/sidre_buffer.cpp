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

using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataBuffer;

using namespace conduit;

//------------------------------------------------------------------------------

TEST(sidre_buffer,create_buffers)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff_0 = ds->createBuffer();
    DataBuffer *dbuff_1 = ds->createBuffer();
    
    EXPECT_EQ(dbuff_0->getUID(), 0);
    EXPECT_EQ(dbuff_1->getUID(), 1);
    ds->destroyBuffer(0);
    
    DataBuffer *dbuff_3 = ds->createBuffer();
    EXPECT_EQ(dbuff_3->getUID(), 0);
    ds->print();
    delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,alloc_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->createBuffer();

    //dbuff->declare(DataType::uint32(10));
    dbuff->declare(asctoolkit::sidre::UINT32_T, 10);
    dbuff->allocate();
    
    uint32 *data_ptr = dbuff->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
    {
        data_ptr[i] = i*i;
    }

    dbuff->getNode().print_detailed();

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              dbuff->getSchema().total_bytes());
  
    ds->print();
    delete ds;
    
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,init_buffer_for_uint32_array)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->allocate(DataType::uint32(10));
    uint32 *data_ptr = dbuff->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
    {
        data_ptr[i] = i*i;
    }

    dbuff->getNode().print_detailed();

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              dbuff->getSchema().total_bytes());

    ds->print();
    delete ds;
    
}


//------------------------------------------------------------------------------

TEST(sidre_buffer,realloc_buffer)
{
    DataStore *ds = new DataStore();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->allocate(DataType::int64(5));
    

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              sizeof(int64)*5);
    
    int64 *data_ptr = dbuff->getNode().as_int64_ptr();
    
    for(int i=0;i<5;i++)
    {
        data_ptr[i] = 5;
    }

    dbuff->getNode().print_detailed();

    dbuff->reallocate(DataType::int64(10));

    // data buffer changes 
    data_ptr = dbuff->getNode().as_int64_ptr();

    for(int i=0;i<5;i++)
    {
        EXPECT_EQ(data_ptr[i],5);
    }

    for(int i=5;i<10;i++)
    {
        data_ptr[i] = 10;
    }


    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              sizeof(int64)*10);
    
    dbuff->getNode().print_detailed();
    
    ds->print();
    delete ds;
    
}



