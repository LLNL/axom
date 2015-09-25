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
//using asctoolkit::sidre::DataType;
//using asctoolkit::slic::setAbortOnError;
//using asctoolkit::slic::setAbortOnAssert;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

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
  EXPECT_FALSE(view);
#if 0
  DataGroup * group = root->createGroup("test");

  EXPECT_TRUE(group->getName() == std::string("test") );

  DataGroup * group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == ATK_NULLPTR);
#endif
  EXPECT_TRUE(true);

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
