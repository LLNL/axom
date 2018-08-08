/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <stdint.h>

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,create_buffers)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_buffer* dbuff_0 = SIDRE_datastore_create_buffer_empty(ds);
  SIDRE_buffer* dbuff_1 = SIDRE_datastore_create_buffer_empty(ds);

  EXPECT_EQ(SIDRE_buffer_get_index(dbuff_0), 0);
  EXPECT_EQ(SIDRE_buffer_get_index(dbuff_1), 1);
  SIDRE_datastore_destroy_buffer(ds, 0);

  SIDRE_buffer* dbuff_3 = SIDRE_datastore_create_buffer_empty(ds);
  EXPECT_EQ(SIDRE_buffer_get_index(dbuff_3), 0);

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_buffer,alloc_buffer_for_int_array)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_buffer* dbuff = SIDRE_datastore_create_buffer_from_type(ds,
                                                                SIDRE_INT_ID,
                                                                10);

//  SIDRE_buffer_declare(dbuff, SIDRE_INT_ID, 10);
  SIDRE_buffer_allocate_existing(dbuff);

  EXPECT_EQ(SIDRE_buffer_get_type_id(dbuff), SIDRE_INT_ID);
  EXPECT_EQ(SIDRE_buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_buffer_get_bytes_per_element(dbuff), sizeof(int));
  EXPECT_EQ(SIDRE_buffer_get_total_bytes(dbuff), sizeof(int) * 10);

  int* data_ptr = (int*) SIDRE_buffer_get_void_ptr(dbuff);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  SIDRE_buffer_print(dbuff);

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_buffer,init_buffer_for_int_array)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_buffer* dbuff = SIDRE_datastore_create_buffer_empty(ds);

  SIDRE_buffer_allocate_from_type(dbuff, SIDRE_INT_ID, 10);

  EXPECT_EQ(SIDRE_buffer_get_type_id(dbuff), SIDRE_INT_ID);
  EXPECT_EQ(SIDRE_buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_buffer_get_bytes_per_element(dbuff), sizeof(int));
  EXPECT_EQ(SIDRE_buffer_get_total_bytes(dbuff), sizeof(int) * 10);

  int* data_ptr = (int*) SIDRE_buffer_get_void_ptr(dbuff);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  SIDRE_buffer_print(dbuff);

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_buffer,realloc_buffer)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_buffer* dbuff = SIDRE_datastore_create_buffer_from_type(ds,
                                                                SIDRE_LONG_ID,
                                                                5);

  SIDRE_buffer_allocate_existing(dbuff);

  EXPECT_EQ(SIDRE_buffer_get_type_id(dbuff), SIDRE_LONG_ID);
  EXPECT_EQ(SIDRE_buffer_get_num_elements(dbuff), 5u);
  EXPECT_EQ(SIDRE_buffer_get_bytes_per_element(dbuff), sizeof(long));
  EXPECT_EQ(SIDRE_buffer_get_total_bytes(dbuff), sizeof(long) * 5);

  long* data_ptr = (long*) SIDRE_buffer_get_void_ptr(dbuff);

  for(int i=0 ; i<5 ; i++)
  {
    data_ptr[i] = 5;
  }

  SIDRE_buffer_print(dbuff);

  SIDRE_buffer_reallocate(dbuff, 10);

  EXPECT_EQ(SIDRE_buffer_get_type_id(dbuff), SIDRE_LONG_ID);
  EXPECT_EQ(SIDRE_buffer_get_num_elements(dbuff), 10u);
  EXPECT_EQ(SIDRE_buffer_get_bytes_per_element(dbuff), sizeof(long));
  EXPECT_EQ(SIDRE_buffer_get_total_bytes(dbuff), sizeof(long) * 10);

  // data buffer changes
  data_ptr = (long*) SIDRE_buffer_get_void_ptr(dbuff);

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(data_ptr[i],5);
  }

  for(int i=5 ; i<10 ; i++)
  {
    data_ptr[i] = 10;
  }

  SIDRE_buffer_print(dbuff);

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);
}
