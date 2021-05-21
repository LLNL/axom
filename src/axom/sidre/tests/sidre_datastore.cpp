// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/core/sidre.hpp"

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::IndexType;

//------------------------------------------------------------------------------

TEST(sidre_datastore, destroy_buffer)
{
  DataStore* ds = new DataStore();
  Buffer* dbuff = ds->createBuffer();

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  EXPECT_EQ(dbuff->getIndex(), bufferIndex);

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_TRUE(ds->getFirstValidBufferIndex() == -1);

  // After destroy, that buffer index should be available again, and have been
  // re-used..
  Buffer* dbuff2 = ds->createBuffer();
  (void)dbuff2;
  std::cerr << ds->getFirstValidBufferIndex() << std::endl;
  std::cerr << bufferIndex << std::endl;
  EXPECT_TRUE(ds->getFirstValidBufferIndex() == bufferIndex);

  // check error condition
  IndexType badBufferIndex = 9999;
  ds->destroyBuffer(badBufferIndex);

  delete ds;
}
