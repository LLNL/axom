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
#include "sidre/SidreConduit.hpp"

#include "conduit/conduit.hpp"

#include <iostream>

//using asctoolkit::sidre::SidreLength;
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
//----------------------------------------------------------------------

// Create a Condit node
// Register with Datastore
// Extract address from Datastore
TEST(sidre_conduit,int_array)
{
  // Create Conduit Node
  int    int_av[6]   = {-2,-4,-8,-16,-32,-64};
  conduit::int_array   int_av_a(int_av,conduit::DataType::c_int(6));
  conduit::Node n;
  n.set(int_av_a);
  n.schema().print();
  int *int_ptr = n.as_int_ptr();
  for(int i=0;i<6;i++)
  {
      EXPECT_EQ(int_ptr[i],int_av[i]);
      std::cout << int_ptr[i] << " ";
  }
  std::cout << std::endl;

  // Add Node to DataStore
  DataStore * ds = new DataStore();

  DataGroup * root = ds->getRoot();

  DataView * view = registerConduitNode(root, "cnode", &n);
  EXPECT_TRUE(view != ATK_NULLPTR);

  int num_elements = view->getNumberOfElements();
  EXPECT_EQ(num_elements, 6);

  int * iptr = (int *) view->getDataPointer();
  for (int i=0; i < 6; i++)
  {
      EXPECT_EQ(int_av[i], iptr[i]);
  }

  delete ds;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Create array in Datastore
// Extract as a Condit node
TEST(sidre_conduit,extract_conduit)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * dv = root->createViewAndBuffer("u0",DataType::c_int(10));
  int * data_ptr = dv->getValue();

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  //  dv->print();


  conduit::Node *node = createConduitNode(dv);
  EXPECT_TRUE(node == ATK_NULLPTR);


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
