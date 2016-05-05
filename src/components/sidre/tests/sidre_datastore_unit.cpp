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
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getGroupIndex("some_other_name") );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getFirstValidGroupIndex() );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getNextValidGroupIndex(0) );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getNextValidGroupIndex(4) );

  EXPECT_EQ( 0UL, dg->getNumViews() );
  EXPECT_FALSE( dg->hasView( -1 ) );
  EXPECT_FALSE( dg->hasView(  0 ) );
  EXPECT_FALSE( dg->hasView(  1 ) );
  EXPECT_FALSE( dg->hasView("some_name") );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getViewIndex("some_other_name") );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getFirstValidViewIndex() );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getNextValidViewIndex(0) );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, dg->getNextValidViewIndex(4) );
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

  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, ds->getFirstValidBufferIndex() );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, ds->getNextValidBufferIndex(0) );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, ds->getNextValidBufferIndex(4) );

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
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, ds->getNextValidBufferIndex(bufferIndex) );

  ds->destroyBuffer(bufferIndex);
  // should be no buffers
  EXPECT_EQ( 0, ds->getNumBuffers() );
  EXPECT_EQ( asctoolkit::sidre::InvalidIndex, ds->getFirstValidBufferIndex() );

  // After destroy, that buffer index should be available again, and have been re-used..
  DataBuffer * dbuff2 = ds->createBuffer( asctoolkit::sidre::FLOAT32_ID, 16 );
  IndexType d2Index = dbuff2->getIndex();
  EXPECT_EQ(bufferIndex, ds->getFirstValidBufferIndex());
  EXPECT_EQ(bufferIndex, d2Index);
  EXPECT_TRUE(ds->hasBuffer(d2Index));

  DataBuffer * dbuff3 = ds->createBuffer();
  IndexType d3Index = dbuff3->getIndex();
  EXPECT_EQ( 2, ds->getNumBuffers() );
  EXPECT_TRUE(ds->hasBuffer(d3Index));

  ds->destroyBuffer(ds->getFirstValidBufferIndex());
  EXPECT_EQ( 1, ds->getNumBuffers() );
  EXPECT_FALSE(ds->hasBuffer(d3Index));

  DataBuffer * dbuff4 = ds->createBuffer();
  DataBuffer * dbuff5 = ds->createBuffer();
  IndexType d4Index = dbuff4->getIndex();
  IndexType d5Index = dbuff5->getIndex();
  EXPECT_EQ( 3, ds->getNumBuffers() );
  EXPECT_TRUE(ds->hasBuffer(d4Index));
  EXPECT_TRUE(ds->hasBuffer(d5Index));

  ds->destroyAllBuffers();
  EXPECT_EQ( 0, ds->getNumBuffers() );
  EXPECT_FALSE(ds->hasBuffer(d2Index));
  EXPECT_FALSE(ds->hasBuffer(d4Index));
  EXPECT_FALSE(ds->hasBuffer(d5Index));

  // check error condition
  IndexType badBufferIndex = 9999;
  ds->destroyBuffer(badBufferIndex);

  delete ds;
}
