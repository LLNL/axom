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
#include "sidre/SidreMalloc.hpp"

using asctoolkit::sidre::SidreLength;
//using asctoolkit::sidre::TypeID;
//using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
//using asctoolkit::sidre::IndexType;
//using asctoolkit::sidre::InvalidIndex;
//using asctoolkit::sidre::isNameValid;
//using asctoolkit::sidre::indexIsValid;
using asctoolkit::sidre::DataType;
//using asctoolkit::slic::setAbortOnError;
//using asctoolkit::slic::setAbortOnAssert;

//----------------------------------------------------------------------
//------------------------------------------------------------------------------

TEST(sidre_static,int_buffer_from_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  int buffer[10];
  for(int i=0 ; i<10 ; i++)
  {
    buffer[i] = i*i;
  }

  DataView * dv = registerStaticNode(root, "snode",
				     static_cast<void *>(buffer),
				     CONDUIT_NATIVE_INT_DATATYPE_ID, 10);

  //  EXPECT_EQ(dv->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dv->getNumberOfElements(), 10u);
  //  int * data_ptr = dv->getValue();
  int * data_ptr = (int *) dv->getDataPointer();

  EXPECT_EQ(buffer, data_ptr);

  //  dv->print();

  //  EXPECT_EQ(dv->getTotalBytes(), sizeof(int) * 10);
  delete ds;

}

//----------------------------------------------------------------------
//------------------------------------------------------------------------------

// Add static array to datastore using overloaded function to
// imply type.

TEST(sidre_static,int_buffer_from_view_overload)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  int buffer[10];
  for(int i=0 ; i<10 ; i++)
  {
    buffer[i] = i*i;
  }

  DataView * dv = registerStaticNode(root, "snode", buffer, 10);

  //  EXPECT_EQ(dv->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dv->getNumberOfElements(), 10u);
  //  int * data_ptr = dv->getValue();
  int * data_ptr = (int *) dv->getDataPointer();

  EXPECT_EQ(buffer, data_ptr);

  //  dv->print();

  //  EXPECT_EQ(dv->getTotalBytes(), sizeof(int) * 10);
  delete ds;

}

//----------------------------------------------------------------------
//------------------------------------------------------------------------------

// Add static scalar to datastore by implying type and length
TEST(sidre_static,int_buffer_from_view_overload_scalar)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  int buffer;
  buffer = 2;

  DataView * dv = registerStaticNode(root, "snode", &buffer);

  //  EXPECT_EQ(dv->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dv->getNumberOfElements(), 1u);
  //  int * data_ptr = dv->getValue();
  int * data_ptr = (int *) dv->getDataPointer();

  EXPECT_EQ(&buffer, data_ptr);

  //  dv->print();

  //  EXPECT_EQ(dv->getTotalBytes(), sizeof(int) * 10);
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
