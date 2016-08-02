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
using asctoolkit::sidre::DataTypeId;

using asctoolkit::sidre::getTypeID;

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
			   DataTypeId tid, int eltsize, int eltcount)
{
  EXPECT_FALSE(buf->isAllocated());
  EXPECT_EQ(isDescribed, buf->isDescribed());
  EXPECT_EQ(tid, buf->getTypeID());
  EXPECT_EQ(eltsize * eltcount, buf->getTotalBytes());
  EXPECT_EQ(eltcount, buf->getNumElements());
  EXPECT_EQ(ATK_NULLPTR, buf->getVoidPtr());
}

void verifyAllocatedBuffer(DataBuffer * buf, DataTypeId tid, int eltsize, int eltcount)
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

  bool isDescribed = false;
  DataTypeId tid = getTypeID(SIDRE_NO_TYPE_ID);
  int eltsize = 0;
  int eltcount = 0;

  DataBuffer * buf1 = ds->createBuffer();
  {
    SCOPED_TRACE("created undescribed");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // invalid number of elements (not yet described)
  buf1->describe(getTypeID(SIDRE_INT_ID), -3);
  {
    SCOPED_TRACE("wrong desc, not yet described, no change");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // Describe to some arbitrary specification
  isDescribed = true;
  tid = getTypeID(SIDRE_INT_ID);
  eltsize = sizeof(int);
  eltcount = 200;
  buf1->describe(tid, eltcount);
  {
    SCOPED_TRACE("described");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // invalid number of elements (already described)
  buf1->describe(getTypeID(SIDRE_INT_ID), -3);
  {
    SCOPED_TRACE("wrong desc, already described, no change");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // change something legal; it should change (because we haven't allocated the buffer yet)
  tid = getTypeID(SIDRE_DOUBLE_ID);
  eltsize = sizeof(double);
  eltcount = 27;
  buf1->describe(tid, eltcount);
  {
    SCOPED_TRACE("legal change");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  //  Repeat all tests on buffer created with some arbitrary description
  tid = getTypeID(SIDRE_FLOAT_ID);
  eltsize = sizeof(float);
  eltcount = 45;
  DataBuffer * buf2 = ds->createBuffer(tid, eltcount);
  {
    SCOPED_TRACE("created described");
    verifyDescribedBuffer(buf2, isDescribed, tid, eltsize, eltcount);
  }

  // invalid number of elements
  buf2->describe(getTypeID(SIDRE_INT_ID), -1);
  {
    SCOPED_TRACE("wrong desc 2, no change");
    verifyDescribedBuffer(buf2, isDescribed, tid, eltsize, eltcount);
  }

  // change something legal; it should change (because we haven't allocated the buffer yet)
  tid = getTypeID(SIDRE_INT_ID);
  eltsize = sizeof(int);
  eltcount = 56;
  buf2->describe(tid, eltcount);
  {
    SCOPED_TRACE("legal change 2");
    verifyDescribedBuffer(buf2, isDescribed, tid, eltsize, eltcount);
  }

  delete ds;
}

// Test allocate methods
TEST(sidre_databuffer,buffer_allocate)
{
  DataStore * ds = new DataStore();

  bool isDescribed = false;
  DataTypeId tid = getTypeID(SIDRE_NO_TYPE_ID);
  int eltsize = 0;
  int eltcount = 0;

  DataBuffer * buf1 = ds->createBuffer();
  verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);

  // Call allocate and reallocate on a non-described buffer: should be no-op
  {
    SCOPED_TRACE("not described: no-op allocate");
    buf1->allocate();
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // Describe to some arbitrary specification
  isDescribed = true;
  tid = getTypeID(SIDRE_INT_ID);
  eltsize = sizeof(int);
  eltcount = 200;
  buf1->describe(tid, eltcount);
  {
    SCOPED_TRACE("describe, allocate");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
    buf1->allocate();
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }

  {
    SCOPED_TRACE("no-op second call to describe");

    // attempt to re-describe an allocated Buffer: should be no-op.
    buf1->describe(getTypeID(SIDRE_INT_ID), 76);
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }

  // Allocate a buffer that was created with a description
  tid = getTypeID(SIDRE_UINT_ID);
  eltsize = sizeof(unsigned int);
  eltcount = 65;
  DataBuffer * buf2 = ds->createBuffer(tid, eltcount);
  {
    SCOPED_TRACE("allocate buffer created with description");
    verifyDescribedBuffer(buf2, isDescribed, tid, eltsize, eltcount);
    buf2->allocate();
    verifyAllocatedBuffer(buf2, tid, eltsize, eltcount);
  }

  tid = getTypeID(SIDRE_FLOAT_ID);
  eltsize = sizeof(float);
  eltcount = 41;
  DataBuffer * buf3 = ds->createBuffer(tid, eltcount);
  {
    SCOPED_TRACE("alloc-describe buffer created with description");
    verifyDescribedBuffer(buf3, isDescribed, tid, eltsize, eltcount);
    
    tid = getTypeID(SIDRE_DOUBLE_ID);
    eltsize = sizeof(double);
    eltcount = 7;
    buf3->allocate(tid, eltcount);
    verifyAllocatedBuffer(buf3, tid, eltsize, eltcount);
  }

  tid = getTypeID(SIDRE_FLOAT_ID);
  eltsize = sizeof(float);
  eltcount = 0;
  DataBuffer * buf4 = ds->createBuffer(tid, eltcount);
  {
    SCOPED_TRACE("allocate a described zero-element buffer");
    verifyDescribedBuffer(buf4, isDescribed, tid, eltsize, eltcount);
    buf4->allocate();
    verifyAllocatedBuffer(buf4, tid, eltsize, eltcount);
  }

  isDescribed = false;
  tid = getTypeID(SIDRE_NO_TYPE_ID);
  DataBuffer * buf4a = ds->createBuffer();
  {
    SCOPED_TRACE("describe and allocate (in one step) a zero-element buffer");
    verifyDescribedBuffer(buf4a, isDescribed, tid, eltsize, eltcount);
    tid = getTypeID(SIDRE_FLOAT_ID);
    buf4a->allocate(tid, eltcount);
    verifyAllocatedBuffer(buf4a, tid, eltsize, eltcount);

    void *ptr = buf4a->getVoidPtr();
    std::cout << "Pointer from zero-allocated buffer is: " << ptr << std::endl;
  }

  eltcount = 12;
  {
    SCOPED_TRACE("REallocate a described zero-element buffer");
    buf4a->reallocate(eltcount);
    eltcount = 0;
    buf4a->reallocate(eltcount);
    verifyAllocatedBuffer(buf4a, tid, eltsize, eltcount);
  }

  delete ds;
}

// Test reallocate methods
TEST(sidre_databuffer,buffer_reallocate)
{
  DataStore * ds = new DataStore();

  bool isDescribed = false;
  DataTypeId tid = getTypeID(SIDRE_NO_TYPE_ID);
  int eltsize = 0;
  int eltcount = 0;

  DataBuffer * buf1 = ds->createBuffer();
  verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);

  {
    SCOPED_TRACE("not described: no-op reallocate(45)");
    buf1->reallocate(45);
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }
  {
    SCOPED_TRACE("not described: no-op reallocate(-2)");
    buf1->reallocate(-2);
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }
  {
    SCOPED_TRACE("not described: no-op reallocate(0)");
    buf1->reallocate(0);
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
  }

  // Describe to some arbitrary specification
  isDescribed = true;
  tid = getTypeID(SIDRE_INT_ID);
  eltsize = sizeof(int);
  eltcount = 200;
  buf1->describe(tid, eltcount);
  {
    SCOPED_TRACE("describe, allocate");
    verifyDescribedBuffer(buf1, isDescribed, tid, eltsize, eltcount);
    buf1->allocate();
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }

  {
    SCOPED_TRACE("reallocate with negative count");
    // reallocation with negative number of elements should not work
    buf1->reallocate(-1);
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }

  // but with zero or more, reallocation should work
  eltcount = 45;
  buf1->reallocate(eltcount);
  {
    SCOPED_TRACE("reallocate(45)");
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }
  eltcount = 0;
  buf1->reallocate(eltcount);
  {
    SCOPED_TRACE("reallocate(0)");
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }
  eltcount = 3;
  buf1->reallocate(eltcount);
  {
    SCOPED_TRACE("reallocate(3)");
    verifyAllocatedBuffer(buf1, tid, eltsize, eltcount);
  }
}

// Test interaction of buffer deletion with views
TEST(sidre_databuffer,buffer_delete_view_detach)
{
  DataStore * ds = new DataStore();

  int vAtest[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  int vB1test[] = { 1, 2, 3, 4, 5, 6, 7 };
  int vB2test[] = { 16, 15, 14, -2 };

  DataTypeId tid = getTypeID(SIDRE_INT_ID);
  int eltsize = sizeof(int);
  int Acount = 8;
  int Bcount = 7;

  // Create a buffer bA, fill with dummy data, attach two views vA and vB
  DataBuffer * bA = ds->createBuffer(tid, Acount)->allocate();
  bA->copyBytesIntoBuffer(vAtest, Acount * eltsize);
  DataView *vA = ds->getRoot()->createView("vA", bA)->apply(tid, Acount);
  DataView *vB = ds->getRoot()->createView("vB", bA)->apply(tid, Bcount, 1, 1);

  // Verify vA and vB are attached to bA, and we can see the data
  EXPECT_EQ(2, bA->getNumViews());
  EXPECT_EQ(bA, vA->getBuffer());
  EXPECT_EQ(0, countMismatch(Acount, vAtest, vA->getArray()));
  EXPECT_EQ(bA, vB->getBuffer());
  EXPECT_EQ(0, countMismatch(Bcount, vB1test, vB->getArray(), true));

  // Detach buffer bA from view vB
  EXPECT_EQ(bA, vB->detachBuffer());

  // Verify that vB now has no buffer but vA still does, and that bA still exists
  EXPECT_EQ(1, bA->getNumViews());
  EXPECT_EQ(bA, vA->getBuffer());
  EXPECT_EQ(0, countMismatch(8, vAtest, vA->getArray()));
  EXPECT_EQ(ATK_NULLPTR, vB->getBuffer());

  // Make a new buffer bB and attach to vB; verify we can see the data in bB
  DataBuffer *bB = ds->createBuffer(tid, Bcount)->allocate();
  bB->copyBytesIntoBuffer(vB1test, Bcount * eltsize);
  vB->attachBuffer(tid, Bcount, bB);
  
  // Detach bA from vA using attach(NULL); verify that vA nas no buffer and bA is destroyed
  IndexType aidx = bA->getIndex();
  vA->attachBuffer(ATK_NULLPTR);
  //  EXPECT_FALSE(vA->isAttached());
  EXPECT_FALSE(ds->hasBuffer(aidx));
  // Destroy bB; verify that vB now has no buffer and cannot access data


  // TODO: test effect of DataBuffer::deallocate();

  delete ds;
}


