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
#include "sidre/SidreVector.hpp"
#include "sidre/SidreConduit.hpp"

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

TEST(sidre_vector,int_buffer_from_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  std::vector<int> var;

  for (int i=0; i < 10; i++)
  {
      var.push_back(i);
  }
  DataView * dv = registerVectorNode(root, "snode", &var);

  //  EXPECT_EQ(dv->getTypeID(), CONDUIT_NATIVE_INT_DATATYPE_ID);
  EXPECT_EQ(dv->getNumberOfElements(), 10u);
  //  int * data_ptr = dv->getValue();
  int * data_ptr = (int *) dv->getDataPointer();

  for(int i=0 ; i<10 ; i++)
  {
      EXPECT_EQ(var[i], data_ptr[i]);
  }

  //  printView(dv);
  dv->print();

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
