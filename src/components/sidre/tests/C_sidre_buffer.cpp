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

#include "sidre/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,create_buffers)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_databuffer *dbuff_0 = ATK_datastore_create_buffer(ds);
    ATK_databuffer *dbuff_1 = ATK_datastore_create_buffer(ds);
    
    EXPECT_EQ(ATK_databuffer_get_uid(dbuff_0), 0);
    EXPECT_EQ(ATK_databuffer_get_uid(dbuff_1), 1);
    ATK_datastore_destroy_buffer(ds, 0);
    
    ATK_databuffer *dbuff_3 = ATK_datastore_create_buffer(ds);
    EXPECT_EQ(ATK_databuffer_get_uid(dbuff_3), 0);
    //    ds->print();
    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,alloc_buffer_for_uint32_array)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_databuffer *dbuff = ATK_datastore_create_buffer(ds);

    ATK_databuffer_declare(dbuff, ATK_UINT32_T, 10);
    ATK_databuffer_allocate(dbuff);
#if 0
    
    uint32 *data_ptr = dbuff->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
    {
        data_ptr[i] = i*i;
    }

    dbuff->getNode().print_detailed();

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              dbuff->getSchema().total_bytes());
  
    ds->print();
#endif
    ATK_datastore_delete(ds);
    
}

#if 0
//------------------------------------------------------------------------------

TEST(C_sidre_buffer,init_buffer_for_uint32_array)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_databuffer *dbuff = ATK_datastore_create_buffer(ds);

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
    ATK_datastore_delete(ds);
    
}


//------------------------------------------------------------------------------

TEST(C_sidre_buffer,realloc_buffer)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_databuffer *dbuff = ATK_datastore_create_buffer(ds);

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
    ATK_datastore_delete(ds);
    
}



#endif
