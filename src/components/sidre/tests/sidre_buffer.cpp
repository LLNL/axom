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
using asctoolkit::sidre::DataType;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataView;
using asctoolkit::sidre::IndexType;

//------------------------------------------------------------------------------


// id's are created in sequence.
// After detaching or destroying a buffer, the next createBuffer will recycle the id.

TEST(sidre_buffer,create_buffers)
{
  DataStore * ds = new DataStore();
  EXPECT_EQ(0u, ds->getNumBuffers());

  // Create two buffers
  DataBuffer * dbuff_0 = ds->createBuffer();
  EXPECT_EQ(1u, ds->getNumBuffers());
  EXPECT_EQ(0, dbuff_0->getIndex());

  DataBuffer * dbuff_1 = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff_1->getIndex());

  //----------
  // destroy by index
  ds->destroyBuffer(0);
  EXPECT_EQ(1u, ds->getNumBuffers());

  DataBuffer * dbuff_0b = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(0, dbuff_0b->getIndex());

  //----------
  // destroy by pointer
  ds->destroyBuffer(dbuff_1);
  EXPECT_EQ(1u, ds->getNumBuffers());

  DataBuffer * dbuff_1b = ds->createBuffer();
  EXPECT_EQ(2u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff_1b->getIndex());

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,create_buffer_with_description)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer(asctoolkit::sidre::INT_ID, 10);
  dbuff->allocate();

  EXPECT_EQ(dbuff->getTypeID(), asctoolkit::sidre::INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 10));

  int * data_ptr = dbuff->getData();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dbuff->print();

  ds->print();

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,alloc_buffer_for_int_array)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  dbuff->allocate(asctoolkit::sidre::INT_ID, 10);
  // Should be a warning and no-op, buffer is already allocated, we don't want
  // to re-allocate and leak memory.
  dbuff->allocate();

  EXPECT_EQ(dbuff->getTypeID(), asctoolkit::sidre::INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 10));

  //  int * data_ptr = static_cast<int *>(dbuff->getData());
  int * data_ptr = dbuff->getData();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dbuff->print();

  ds->print();

  dbuff->deallocate();

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,init_buffer_for_int_array)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  dbuff->allocate(asctoolkit::sidre::INT_ID, 10);

  EXPECT_EQ(dbuff->getTypeID(), asctoolkit::sidre::INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 10));

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

  dbuff->allocate(asctoolkit::sidre::INT_ID, 5);

  EXPECT_EQ(dbuff->getTypeID(), asctoolkit::sidre::INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 5u);
  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 5));

  int * data_ptr = static_cast<int *>(dbuff->getVoidPtr());

  for(int i=0 ; i<5 ; i++)
  {
    data_ptr[i] = 5;
  }

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(data_ptr[i],5);
  }

//  dbuff->print();

  dbuff->reallocate(10);

  std::cerr << dbuff->getTypeID() << std::endl;
  EXPECT_EQ(dbuff->getTypeID(), asctoolkit::sidre::INT_ID);
  EXPECT_EQ(dbuff->getNumElements(), 10u);
  EXPECT_EQ(dbuff->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(sizeof(int) * 10));

  // data buffer changes
  data_ptr = static_cast<int *>(dbuff->getVoidPtr());

  for(int i=5 ; i<10 ; i++)
  {
    data_ptr[i] = 10;
  }

  for(int i=0 ; i<10 ; i++)
  {
    int value = 5;
    if ( i > 4)
    {
      value = 10;
    }
    EXPECT_EQ(data_ptr[i],value);
  }

//  dbuff->print();

//  ds->print();
  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_buffer,with_multiple_views)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff;
  DataView * dv1, * dv2;

  dbuff = ds->createBuffer();
  IndexType idx = dbuff->getIndex();

  dv1 = root->createView("e1", dbuff);
  dv2 = root->createView("e2", dbuff);
  EXPECT_TRUE(dv1->hasBuffer());
  EXPECT_TRUE(dv2->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 1u);
  EXPECT_EQ(dbuff->getNumViews(), 2);

  // Detach buffer from first view will not detach from datastore.
  dv1->attachBuffer(NULL);
  EXPECT_FALSE(dv1->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 1u);
  EXPECT_EQ(dbuff->getNumViews(), 1);

  // Detach buffer from second view will detach from datastore.
  dv2->attachBuffer(NULL);
  EXPECT_FALSE(dv2->hasBuffer());
  EXPECT_EQ(ds->getNumBuffers(), 0u);

  // Buffer has been destroyed since there are no more attached views
  EXPECT_TRUE(ds->getBuffer(idx) == NULL);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_buffer,move_buffer)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataBuffer * dbuff, * dbuff2;
  DataView * dv1, * dv2;

  dbuff = ds->createBuffer();

  dv1 = root->createView("e1", dbuff);

  dbuff2 = dv1->detachBuffer();
  EXPECT_FALSE(dv1->hasBuffer());
  EXPECT_EQ(dbuff2, dbuff);

  dv2 = root->createView("e2");
  dv2->attachBuffer(dbuff2);

  EXPECT_TRUE(dv2->hasBuffer());
  EXPECT_EQ(1u, ds->getNumBuffers());
  EXPECT_EQ(1, dbuff->getNumViews());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_buffer,destroy_all_buffers)
{
  DataStore * ds = new DataStore();

  ds->createBuffer();
  ds->createBuffer();
  ds->createBuffer();
  ds->createBuffer();
  EXPECT_EQ(4u, ds->getNumBuffers());

  ds->destroyAllBuffers();
  EXPECT_EQ(0u, ds->getNumBuffers());

  delete ds;
}
