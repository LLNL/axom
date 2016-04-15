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
using asctoolkit::sidre::IndexType;

//------------------------------------------------------------------------------

TEST(sidre_datastore,destroy_buffer)
{
  DataStore * ds = new DataStore();
  DataBuffer * dbuff = ds->createBuffer();

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  EXPECT_EQ(dbuff->getIndex(), bufferIndex);

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_TRUE( ds->getFirstValidBufferIndex() == -1 );

  // After destroy, that buffer index should be available again, and have been re-used..
  DataBuffer * dbuff2 = ds->createBuffer();
  (void)dbuff2;
  std::cerr << ds->getFirstValidBufferIndex() << std::endl;
  std::cerr << bufferIndex << std::endl;
  EXPECT_TRUE( ds->getFirstValidBufferIndex() == bufferIndex );

  // check error condition
  IndexType badBufferIndex = 9999;
  ds->destroyBuffer(badBufferIndex);

  delete ds;
}
