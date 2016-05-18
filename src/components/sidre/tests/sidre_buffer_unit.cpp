/*
 * Copyright (c) 2016, Lawrence Livermore National Security, LLC.
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
using asctoolkit::sidre::DataView;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::IndexType;
using asctoolkit::sidre::InvalidIndex;

#include <map>

//------------------------------------------------------------------------------

int countMismatch(unsigned int elts, int * standard, int * undertest, bool printTest = false);

int countMismatch(unsigned int elts, int * standard, int * undertest, bool printTest)
{
  int retval = 0;

  for (unsigned int i = 0; i < elts; ++i)
  {
    if (standard[i] != undertest[i]) ++retval;
    if (printTest) std::cout << "  " << undertest[i];
  }
  if (printTest) std::cout << std::endl;

  return retval;
}
 
// Ctor tests
// Create buffer, verify its index is what is expected, verify it has zero elements 
// and zero total bytes, verify it has zero views (hasView(idx) returns false, 
// getView(idx) returns null ptr).
TEST(sidre_databuffer, buffer_create)
{
  DataStore * ds = new DataStore();

  DataBuffer * buf1 = ds->createBuffer();

  EXPECT_EQ(ATK_NULLPTR, buf1->getVoidPtr());
  // EXPECT_EQ(ATK_NULLPTR, buf1->getData());
  EXPECT_EQ(0, buf1->getIndex());
  EXPECT_EQ(0, buf1->getNumViews());
  EXPECT_EQ(0, buf1->getTotalBytes());
  EXPECT_FALSE(buf1->isAllocated());
  EXPECT_FALSE(buf1->isDescribed());

  delete ds;
}

// Test interaction of buffer deletion with views
TEST(sidre_databuffer,buffer_delete_view_detach)
{
  DataStore * ds = new DataStore();

  int vAtest[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  int vB1test[] = { 1, 3, 5, 7 };
  int vB2test[] = { 16, 15, 14, -2 };

  // Create a buffer bA, fill with dummy data, attach two views vA and vB
  DataBuffer * bA = ds->createBuffer(asctoolkit::sidre::INT32_ID, 8)->allocate();
  bA->copyBytesIntoBuffer(vAtest, 8 * sizeof(int));
  DataView *vA = ds->getRoot()->createView("vA", bA)->apply(asctoolkit::sidre::INT32_ID, 8);
  DataView *vB = ds->getRoot()->createView("vB", bA)->apply(asctoolkit::sidre::INT32_ID, 4, 1, 2);

  // Verify vA and vB are attached to bA, and we can see the data
  EXPECT_EQ(2, bA->getNumViews());
  EXPECT_EQ(bA, vA->getBuffer());
  EXPECT_EQ(0, countMismatch(8, vAtest, vA->getArray()));
  EXPECT_EQ(bA, vB->getBuffer());
  EXPECT_EQ(0, countMismatch(4, vB1test, vB->getArray(), true));

  // Detach buffer bA from view vB
  EXPECT_EQ(bA, vB->detachBuffer());

  // Verify that vB now has no buffer but vA still does, and that bA still exists
  EXPECT_EQ(1, bA->getNumViews());
  EXPECT_EQ(bA, vA->getBuffer());
  EXPECT_EQ(0, countMismatch(8, vAtest, vA->getArray()));
  EXPECT_EQ(ATK_NULLPTR, vB->getBuffer());

  // Make a new buffer bB and attach to vB; verify we can see the data in bB
  // Detach bA from vA; verify that vA nas no buffer and bA is destroyed
  // Destroy bB; verify that vB now has no buffer and cannot access data


  // TODO: test effect of DataBuffer::deallocate();

  delete ds;
}


