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
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::IndexType;
using asctoolkit::sidre::InvalidIndex;

#include <map>
#include <random>

//------------------------------------------------------------------------------

// This is a group of unit tests to run through the Sidre DataStore API, and
// test for pre and post conditions.  The tests will only use public APIs unless
// I'm told differently.

// Public APIs to exercise:
//
//   DataStore();
//   ~DataStore();
//   DataBuffer * createBuffer();
//   DataBuffer * createBuffer( TypeID type, SidreLength num_elems );
//   void destroyBuffer( DataBuffer * buff );
//   void destroyBuffer( IndexType idx );
//   void destroyAllBuffers();

// Public APIs for getting state:
//
//   DataGroup * getRoot()
//   size_t getNumBuffers() const
//   bool hasBuffer( IndexType idx ) const
//   DataBuffer * getBuffer( IndexType idx ) const;
//   IndexType getFirstValidBufferIndex() const;
//   IndexType getNextValidBufferIndex(IndexType idx) const;
//   void info(Node& n) const;

// I will assume print() works, and I don't plan to use it.

void verifyEmptyGroupNamed( DataGroup * dg, std::string name )
{
  EXPECT_EQ( name, dg->getName() );
 
  EXPECT_EQ( 0UL, dg->getNumGroups() );
  EXPECT_FALSE( dg->hasGroup( -1 ) );
  EXPECT_FALSE( dg->hasGroup(  0 ) );
  EXPECT_FALSE( dg->hasGroup(  1 ) );
  EXPECT_FALSE( dg->hasGroup("some_name") );
  EXPECT_EQ( InvalidIndex, dg->getGroupIndex("some_other_name") );
  EXPECT_EQ( InvalidIndex, dg->getFirstValidGroupIndex() );
  EXPECT_EQ( InvalidIndex, dg->getNextValidGroupIndex(0) );
  EXPECT_EQ( InvalidIndex, dg->getNextValidGroupIndex(4) );

  EXPECT_EQ( 0UL, dg->getNumViews() );
  EXPECT_FALSE( dg->hasView( -1 ) );
  EXPECT_FALSE( dg->hasView(  0 ) );
  EXPECT_FALSE( dg->hasView(  1 ) );
  EXPECT_FALSE( dg->hasView("some_name") );
  EXPECT_EQ( InvalidIndex, dg->getViewIndex("some_other_name") );
  EXPECT_EQ( InvalidIndex, dg->getFirstValidViewIndex() );
  EXPECT_EQ( InvalidIndex, dg->getNextValidViewIndex(0) );
  EXPECT_EQ( InvalidIndex, dg->getNextValidViewIndex(4) );
}

void verifyBufferIdentity(DataStore * ds, std::map<IndexType, DataBuffer *> & bs)
{
  int bufcount = bs.size();

  // Does ds contain the number of buffers we expect?
  EXPECT_EQ(bufcount, ds->getNumBuffers());

  // Does ds contain the buffer IDs and pointers we expect?
  int iteratedCount = 0;
  IndexType idx = ds->getFirstValidBufferIndex();
  while(InvalidIndex != idx && iteratedCount < bufcount) {
    EXPECT_TRUE(bs.count(idx) == 1);
    if (bs.count(idx) == 1) {
      EXPECT_EQ(bs[idx], ds->getBuffer(idx));
    }
    idx = ds->getNextValidBufferIndex(idx);
    iteratedCount += 1;
  }

  // Have we iterated over exactly the number of buffers we expect, finishing on InvalidIndex?
  EXPECT_EQ(bufcount, iteratedCount);
  EXPECT_EQ(InvalidIndex, idx);
}

TEST(sidre_datastore,default_ctor)
{
  DataStore * ds = new DataStore();

  // After construction, the DataStore should contain no buffers.

  EXPECT_EQ( 0, ds->getNumBuffers() );
  EXPECT_FALSE( ds->hasBuffer(-15) );
  EXPECT_FALSE( ds->hasBuffer(-1) );
  EXPECT_FALSE( ds->hasBuffer(0) );
  EXPECT_FALSE( ds->hasBuffer(1) );
  EXPECT_FALSE( ds->hasBuffer(8) );

  EXPECT_EQ( InvalidIndex, ds->getFirstValidBufferIndex() );
  EXPECT_EQ( InvalidIndex, ds->getNextValidBufferIndex(0) );
  EXPECT_EQ( InvalidIndex, ds->getNextValidBufferIndex(4) );

  // The new DataStore should contain exactly one group, the root group.
  // The root group should be named "/" and should contain no views and no groups.

  DataGroup * dg = ds->getRoot();
  
  ASSERT_NE( ATK_NULLPTR, dg );
  EXPECT_EQ( ATK_NULLPTR, dg->getParent() );
  EXPECT_EQ( ds, dg->getDataStore() );

  verifyEmptyGroupNamed(dg, "/");

  delete ds;
}

// The dtor destroys all buffers and deletes the root group.
// I suppose I can get pointers to the root group and some buffers, but how do I
// test to see if the pointer is invalid?  

TEST(sidre_datastore,create_destroy_buffers)
{
  DataStore * ds = new DataStore();
  EXPECT_EQ( 0, ds->getNumBuffers() );

  DataBuffer * dbuff = ds->createBuffer();
  EXPECT_EQ( 1, ds->getNumBuffers() );

  IndexType bufferIndex = ds->getFirstValidBufferIndex();
  EXPECT_EQ( dbuff->getIndex(), bufferIndex);
  EXPECT_EQ( InvalidIndex, ds->getNextValidBufferIndex(bufferIndex) );

  // Do we get the buffer we expect?
  EXPECT_EQ(dbuff, ds->getBuffer(bufferIndex));
  IndexType badBufferIndex = 9999;
  EXPECT_EQ(ATK_NULLPTR, ds->getBuffer(badBufferIndex));

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_EQ( 0, ds->getNumBuffers() );
  EXPECT_EQ( InvalidIndex, ds->getFirstValidBufferIndex() );
  EXPECT_FALSE(ds->hasBuffer(bufferIndex));
  EXPECT_EQ(ATK_NULLPTR, ds->getBuffer(bufferIndex));
  EXPECT_EQ(ATK_NULLPTR, ds->getBuffer(badBufferIndex));

  // After destroy, that buffer index should be available again, and have been re-used..
  DataBuffer * dbuff2 = ds->createBuffer( asctoolkit::sidre::FLOAT32_ID, 16 );
  IndexType d2Index = dbuff2->getIndex();
  EXPECT_EQ(bufferIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(bufferIndex, d2Index);
  EXPECT_TRUE(ds->hasBuffer(d2Index));
  EXPECT_EQ(dbuff2, ds->getBuffer(bufferIndex));
  EXPECT_EQ(ATK_NULLPTR, ds->getBuffer(badBufferIndex));

  DataBuffer * dbuff3 = ds->createBuffer();
  IndexType d3Index = dbuff3->getIndex();
  EXPECT_EQ( 2, ds->getNumBuffers() );
  EXPECT_TRUE(ds->hasBuffer(d3Index));

  // Try destroying the first one; see if we have the correct count and indices
  ds->destroyBuffer(ds->getFirstValidBufferIndex());
  EXPECT_EQ( 1, ds->getNumBuffers() );
  EXPECT_FALSE(ds->hasBuffer(d2Index));
  EXPECT_TRUE(ds->hasBuffer(d3Index));

  // Add some more buffers, then try destroying the second one; see if we have
  // the correct count and indices
  DataBuffer * dbuff4 = ds->createBuffer();
  DataBuffer * dbuff5 = ds->createBuffer();
  IndexType d4Index = dbuff4->getIndex();
  IndexType d5Index = dbuff5->getIndex();
  EXPECT_EQ( 3, ds->getNumBuffers() );
  EXPECT_TRUE(ds->hasBuffer(d3Index));
  EXPECT_TRUE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));

  ds->destroyBuffer(d4Index);
  EXPECT_EQ( 2, ds->getNumBuffers() );
  EXPECT_TRUE(ds->hasBuffer(d3Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));

  // Can we destroy all buffers?
  ds->destroyAllBuffers();
  EXPECT_EQ( 0, ds->getNumBuffers() );
  EXPECT_FALSE(ds->hasBuffer(d2Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_FALSE(ds->hasBuffer(d5Index));

  // check error condition
  ds->destroyBuffer(badBufferIndex);

  delete ds;
}

// Test iteration through buffers, as well as proper index and buffer behavior
// while buffers are created and deleted
TEST(sidre_datastore,iterate_buffers)
{
  DataStore * ds = new DataStore();
  EXPECT_EQ( 0, ds->getNumBuffers() );

  IndexType badBufferIndex = 9999;
  // Do we get sidre::InvalidIndex for several queries with no buffers?
  EXPECT_EQ(InvalidIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(0));
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(badBufferIndex));
  EXPECT_EQ(InvalidIndex, ds->getNextValidBufferIndex(InvalidIndex));

  std::map<IndexType, DataBuffer *> bs;
  int bufcount = 0;

  {
    SCOPED_TRACE("simple 1");
    ds->destroyAllBuffers();
    bs.clear();
    bufcount = 20;

    for (int i = 0; i < bufcount; ++i)
    {
      DataBuffer *b = ds->createBuffer(asctoolkit::sidre::FLOAT64_ID, 400*i);
      IndexType idx = b->getIndex();
      bs[idx] = b;
    }

    verifyBufferIdentity(ds, bs);
  }

  {
    SCOPED_TRACE("simple 2");
    ds->destroyAllBuffers();
    bs.clear();
    bufcount = 50;

    for (int i = 0; i < bufcount; ++i)
    {
      DataBuffer *b = ds->createBuffer(asctoolkit::sidre::FLOAT64_ID, 400*i % 10000);
      IndexType idx = b->getIndex();
      bs[idx] = b;
    }

    int i = 0;
    std::map<IndexType, DataBuffer *> nbs;
    std::map<IndexType, DataBuffer *>::iterator bsit = bs.begin(), bsend = bs.end();
    for (; bsit != bsend; ++bsit)
    {
      if (i % 5 && i % 7)
      {
        nbs[bsit->first] = bsit->second;
      }
      else
      {
	ds->destroyBuffer(bsit->first);
      }
      i += 1;
    }

    verifyBufferIdentity(ds, nbs);
  }

  {
    SCOPED_TRACE("random ");
    ds->destroyAllBuffers();
    bs.clear();
    std::vector<int> idxlist;
    int initbufcount = 50;

    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::poisson_distribution<int> poi(3);

    for (int i = 0; i < initbufcount; ++i)
    {
      DataBuffer *b = ds->createBuffer(asctoolkit::sidre::FLOAT64_ID, 400*i % 10000);
      IndexType idx = b->getIndex();
      bs[idx] = b;
      idxlist.push_back(idx);
    }

    int totalrounds = 100;
    for (int round = 0; round < totalrounds; ++round)
    {
      SCOPED_TRACE(round);

      {
        SCOPED_TRACE("remove");
        // Remove a few 
        int rmvcount = poi(rng);
        if (rmvcount > bs.size()) rmvcount = bs.size();
        for (int i = 0; i < rmvcount; ++i)
        {
          std::uniform_int_distribution<int> uni(0, idxlist.size()-1);
          int rmvidx = uni(rng);
          int rmvid = idxlist[rmvidx];
          EXPECT_TRUE(ds->hasBuffer(rmvid));
          EXPECT_TRUE(bs.count(rmvid) == 1);
          ds->destroyBuffer(rmvid);
          bs.erase(rmvid);
          idxlist.erase(idxlist.begin() + rmvidx);
          EXPECT_FALSE(ds->hasBuffer(rmvid));
          EXPECT_FALSE(bs.count(rmvid) == 1);
        }
      }

      {
        SCOPED_TRACE("add");
        // Add a few 
        int addcount = poi(rng);
        for (int i = 0; i < addcount; ++i)
        {
          DataBuffer *buf = ds->createBuffer(asctoolkit::sidre::FLOAT64_ID, 400);
          int addid = buf->getIndex();
          EXPECT_TRUE(ds->hasBuffer(addid));
          EXPECT_TRUE(bs.count(addid) < 1);
          bs[addid] = buf;
          idxlist.push_back(addid);
        }
      }

      verifyBufferIdentity(ds, bs);
    }
  }

  delete ds;
}

