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

TEST(sidre_malloc,int_buffer_from_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv = registerMallocNode(root, "snode");

  dv->allocate(CONDUIT_NATIVE_INT_DATATYPE_ID, 10);
  //  EXPECT_EQ(dv->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  //  int * data_ptr = dv->getValue();
  int * data_ptr = (int *) dv->getDataPointer();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  dv->print();

  //  EXPECT_EQ(dv->getTotalBytes(), sizeof(int) * 10);
  delete ds;

}

//----------------------------------------------------------------------

// Create a Malloc MetaBuffer node
// Register with Datastore
// Extract address from Datastore
TEST(sidre_simple,int_array)
{
  SidreLength nitems = 10;

  DataStore * ds = new DataStore();

  DataGroup * root = ds->getRoot();

  DataView * view = registerMallocNode(root, "snode", CONDUIT_NATIVE_INT_DATATYPE_ID, nitems);
  EXPECT_TRUE(view != ATK_NULLPTR);
  
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
