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

#include "conduit/conduit.hpp"

using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataBuffer;

//------------------------------------------------------------------------------

TEST(sidre_buffer,create_buffers)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff_0 = ds->createBuffer();
  DataBuffer * dbuff_1 = ds->createBuffer();

  EXPECT_EQ(dbuff_0->getIndex(), 0);
  EXPECT_EQ(dbuff_1->getIndex(), 1);
  ds->destroyBuffer(0);

  DataBuffer * dbuff_3 = ds->createBuffer();
  EXPECT_EQ(dbuff_3->getIndex(), 0);

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,alloc_buffer_for_int_array)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  dbuff->allocate(CONDUIT_NATIVE_INT_DATATYPE_ID, 10);
  dbuff->allocate();

  EXPECT_EQ(dbuff->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dbuff->getNumberOfElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), sizeof(int) * 10);

  int * data_ptr = static_cast<int *>(dbuff->getData());

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dbuff->print();

  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,init_buffer_for_int_array)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  dbuff->allocate(CONDUIT_NATIVE_INT_DATATYPE_ID, 10);

  EXPECT_EQ(dbuff->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dbuff->getNumberOfElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), sizeof(int) * 10);

  int * data_ptr = static_cast<int *>(dbuff->getData());

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dbuff->print();

  ds->print();
  delete ds;
}


//------------------------------------------------------------------------------

TEST(sidre_buffer,realloc_buffer)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  dbuff->allocate(CONDUIT_NATIVE_LONG_DATATYPE_ID, 5);

  EXPECT_EQ(dbuff->getTypeID(), CONDUIT_NATIVE_LONG_DATATYPE_ID);
  EXPECT_EQ(dbuff->getNumberOfElements(), 5u);
  EXPECT_EQ(dbuff->getTotalBytes(), sizeof(long) * 5);

  long * data_ptr = static_cast<long *>(dbuff->getData());

  for(int i=0 ; i<5 ; i++)
  {
    data_ptr[i] = 5;
  }

  dbuff->print();

  dbuff->reallocate(CONDUIT_NATIVE_LONG_DATATYPE_ID, 10);

  EXPECT_EQ(dbuff->getTypeID(), CONDUIT_NATIVE_LONG_DATATYPE_ID);
  EXPECT_EQ(dbuff->getNumberOfElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(), sizeof(long) * 10);

  // data buffer changes
  data_ptr = static_cast<long *>(dbuff->getData());

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(data_ptr[i],5);
  }

  for(int i=5 ; i<10 ; i++)
  {
    data_ptr[i] = 10;
  }

  dbuff->print();

  ds->print();
  delete ds;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
