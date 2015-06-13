/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include <stdint.h>

#include "gtest/gtest.h"

#include "sidre/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,create_buffers)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_databuffer * dbuff_0 = ATK_datastore_create_buffer(ds);
  ATK_databuffer * dbuff_1 = ATK_datastore_create_buffer(ds);

  EXPECT_EQ(ATK_databuffer_get_index(dbuff_0), 0);
  EXPECT_EQ(ATK_databuffer_get_index(dbuff_1), 1);
  ATK_datastore_destroy_buffer(ds, 0);

  ATK_databuffer * dbuff_3 = ATK_datastore_create_buffer(ds);
  EXPECT_EQ(ATK_databuffer_get_index(dbuff_3), 0);
  //    ds->print();
  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,alloc_buffer_for_int_array)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_databuffer * dbuff = ATK_datastore_create_buffer(ds);

  ATK_databuffer_declare(dbuff, ATK_C_INT_T, 10);
  ATK_databuffer_allocate_existing(dbuff);

  //    uint32_t *data_ptr = ATK_databuffer_get_data(dbuff);
  int * data_ptr = (int *) ATK_databuffer_get_data(dbuff);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

#if 0
  dbuff->getNode().print_detailed();

  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
            dbuff->getSchema().total_bytes());

  ds->print();
#endif
  ATK_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,init_buffer_for_int_array)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_databuffer * dbuff = ATK_datastore_create_buffer(ds);

  ATK_databuffer_allocate_from_type(dbuff, ATK_C_INT_T, 10);
  int * data_ptr = (int *) ATK_databuffer_get_data(dbuff);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

#ifdef XXX
  dbuff->getNode().print_detailed();

  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
            dbuff->getSchema().total_bytes());
#endif
  ATK_datastore_print(ds);
  ATK_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,realloc_buffer)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_databuffer * dbuff = ATK_datastore_create_buffer(ds);

  ATK_databuffer_declare(dbuff, ATK_C_LONG_T, 5);
  ATK_databuffer_allocate_existing(dbuff);

#if 0
  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
            sizeof(long)*5);
#endif

  long * data_ptr = (long *) ATK_databuffer_get_data(dbuff);

  for(int i=0 ; i<5 ; i++)
  {
    data_ptr[i] = 5;
  }

#ifdef XXX
  dbuff->getNode().print_detailed();
#endif

  ATK_databuffer_reallocate(dbuff, ATK_C_LONG_T, 10);

  // data buffer changes
  data_ptr = (long *) ATK_databuffer_get_data(dbuff);

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(data_ptr[i],5);
  }

  for(int i=5 ; i<10 ; i++)
  {
    data_ptr[i] = 10;
  }

#ifdef XXX
  EXPECT_EQ(ATK_databuffer_get_total_bytes(dbuff), sizeof(long)*10);
  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
            sizeof(long)*10);

  dbuff->getNode().print_detailed();

  ATK_datastore_print(ds);
#endif
  ATK_datastore_delete(ds);

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
   int result = 0;

   ::testing::InitGoogleTest(&argc, argv);

   UnitTestLogger logger;  // create & initialize test logger,
                       // finalized when exiting main scope

   result = RUN_ALL_TESTS();

   return result;
}
