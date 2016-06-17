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

void verifyDescribedBuffer(DataBuffer * buf, bool isDescribed, 
			   asctoolkit::sidre::DataTypeId tid, int eltsize, int eltcount)
{
  EXPECT_FALSE(buf->isAllocated());
  EXPECT_EQ(isDescribed, buf->isDescribed());
  EXPECT_EQ(tid, buf->getTypeID());
  EXPECT_EQ(eltsize * eltcount, buf->getTotalBytes());
  EXPECT_EQ(eltcount, buf->getNumElements());
  EXPECT_EQ(ATK_NULLPTR, buf->getVoidPtr());
}

void verifyAllocatedBuffer(DataBuffer * buf, asctoolkit::sidre::DataTypeId tid, int eltsize, int eltcount)
{
  EXPECT_TRUE(buf->isAllocated());
  EXPECT_TRUE(buf->isDescribed());
  EXPECT_EQ(tid, buf->getTypeID());
  EXPECT_EQ(eltsize * eltcount, buf->getTotalBytes());
  EXPECT_EQ(eltcount, buf->getNumElements());
  EXPECT_NE(ATK_NULLPTR, buf->getVoidPtr());
}

// Test describe methods
TEST(sidre_databuffer,buffer_describe)
{
  DataStore * ds = new DataStore();

  DataBuffer * buf1 = ds->createBuffer();
  {
    SCOPED_TRACE("created undescribed");
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }

  // invalid number of elements (not yet described)
  buf1->describe(asctoolkit::sidre::INT32_ID, -3);
  {
    SCOPED_TRACE("wrong desc, not yet described, no change");
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }

  // Describe to some arbitrary specification
  buf1->describe(asctoolkit::sidre::INT32_ID, 200);
  {
    SCOPED_TRACE("described");
    verifyDescribedBuffer(buf1, true, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  // invalid number of elements (already described)
  buf1->describe(asctoolkit::sidre::INT32_ID, -3);
  {
    SCOPED_TRACE("wrong desc, already described, no change");
    verifyDescribedBuffer(buf1, true, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  // change something legal; it should change (because we haven't allocated the buffer yet)
  buf1->describe(asctoolkit::sidre::FLOAT64_ID, 27);
  {
    SCOPED_TRACE("legal change");
    verifyDescribedBuffer(buf1, true, asctoolkit::sidre::FLOAT64_ID, 8, 27);
  }

  //  Repeat all tests on buffer created with some arbitrary description
  DataBuffer * buf2 = ds->createBuffer(asctoolkit::sidre::UINT16_ID, 45);
  {
    SCOPED_TRACE("created described");
    verifyDescribedBuffer(buf2, true, asctoolkit::sidre::UINT16_ID, 2, 45);
  }

  // invalid number of elements
  buf2->describe(asctoolkit::sidre::INT32_ID, -1);
  {
    SCOPED_TRACE("wrong desc 2, no change");
    verifyDescribedBuffer(buf2, true, asctoolkit::sidre::UINT16_ID, 2, 45);
  }

  // change something legal; it should change (because we haven't allocated the buffer yet)
  buf2->describe(asctoolkit::sidre::FLOAT32_ID, 56);
  {
    SCOPED_TRACE("legal change 2");
    verifyDescribedBuffer(buf2, true, asctoolkit::sidre::FLOAT32_ID, 4, 56);
  }

  delete ds;
}

// Test allocate methods
TEST(sidre_databuffer,buffer_allocate)
{
  DataStore * ds = new DataStore();

  DataBuffer * buf1 = ds->createBuffer();
  verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);

  // Call allocate and reallocate on a non-described buffer: should be no-op
  {
    SCOPED_TRACE("not described: no-op allocate");
    buf1->allocate();
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }

  // Describe to some arbitrary specification
  {
    SCOPED_TRACE("describe, allocate");
    buf1->describe(asctoolkit::sidre::INT32_ID, 200);
    verifyDescribedBuffer(buf1, true, asctoolkit::sidre::INT32_ID, 4, 200);
    buf1->allocate();
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  {
    SCOPED_TRACE("no-op second call to describe");

    // attempt to re-describe an allocated Buffer: should be no-op.
    buf1->describe(asctoolkit::sidre::UINT32_ID, 76);
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  // Allocate a buffer that was created with a description
  {
    SCOPED_TRACE("allocate buffer created with description");
    DataBuffer * buf2 = ds->createBuffer(asctoolkit::sidre::INT16_ID, 65);
    verifyDescribedBuffer(buf2, true, asctoolkit::sidre::INT16_ID, 2, 65);
    buf2->allocate();
    verifyAllocatedBuffer(buf2, asctoolkit::sidre::INT16_ID, 2, 65);
  }

  {
    SCOPED_TRACE("alloc-describe buffer created with description");
    DataBuffer * buf3 = ds->createBuffer(asctoolkit::sidre::INT16_ID, 41);
    verifyDescribedBuffer(buf3, true, asctoolkit::sidre::INT16_ID, 2, 41);
    buf3->allocate(asctoolkit::sidre::FLOAT64_ID, 7);
    verifyAllocatedBuffer(buf3, asctoolkit::sidre::FLOAT64_ID, 4, 7);
  }

  delete ds;
}

// Test reallocate methods
TEST(sidre_databuffer,buffer_reallocate)
{
  DataStore * ds = new DataStore();

  DataBuffer * buf1 = ds->createBuffer();
  verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);

  {
    SCOPED_TRACE("not described: no-op reallocate(45)");
    buf1->reallocate(45);
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }
  {
    SCOPED_TRACE("not described: no-op reallocate(-2)");
    buf1->reallocate(-2);
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }
  {
    SCOPED_TRACE("not described: no-op reallocate(0)");
    buf1->reallocate(0);
    verifyDescribedBuffer(buf1, false, asctoolkit::sidre::NO_TYPE_ID, 0, 0);
  }

  // Describe to some arbitrary specification
  {
    SCOPED_TRACE("describe, allocate");
    buf1->describe(asctoolkit::sidre::INT32_ID, 200);
    verifyDescribedBuffer(buf1, true, asctoolkit::sidre::INT32_ID, 4, 200);
    buf1->allocate();
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  {
    SCOPED_TRACE("reallocate with negative count");
    // reallocation with negative number of elements should not work
    buf1->reallocate(-1);
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 200);
  }

  // but with zero or more, reallocation should work
  {
    SCOPED_TRACE("reallocate(45)");
    buf1->reallocate(45);
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 45);
  }
  {
    SCOPED_TRACE("reallocate(0)");  // ???
    buf1->reallocate(0);
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 0);
  }
  {
    SCOPED_TRACE("reallocate(3)");
    buf1->reallocate(3);
    verifyAllocatedBuffer(buf1, asctoolkit::sidre::INT32_ID, 4, 3);
  }
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
  // EXPECT_EQ(0, countMismatch(4, vB1test, vB->getArray(), true));

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


