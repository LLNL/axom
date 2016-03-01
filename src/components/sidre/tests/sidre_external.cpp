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

#include <vector>

#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
using asctoolkit::sidre::SidreLength;
using asctoolkit::sidre::DOUBLE_ID;
using asctoolkit::sidre::INT_ID;

//------------------------------------------------------------------------------
// Test DataGroup::createView() -- external
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  const SidreLength len = 11;

  int * idata = new int[len];
  double * ddata = new double[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  DataView * iview = root->createView("idata", idata)->apply(INT_ID, len);
  DataView * dview = root->createView("ddata", ddata)->apply(DOUBLE_ID, len);
  EXPECT_EQ(root->getNumViews(), 2u);

  EXPECT_EQ(iview->getNumElements(), len);
  EXPECT_TRUE(iview->isApplied());
  EXPECT_TRUE(iview->isAllocated());
  EXPECT_TRUE(iview->isExternal());
  EXPECT_FALSE(iview->isOpaque());

  EXPECT_EQ(dview->getNumElements(), len);
  EXPECT_TRUE(dview->isApplied());
  EXPECT_TRUE(dview->isAllocated());
  EXPECT_TRUE(dview->isExternal());
  EXPECT_FALSE(dview->isOpaque());

  iview->print();
  dview->print();

  int * idata_chk = iview->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = dview->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  delete ds;
  delete [] idata;
  delete [] ddata;
}

//------------------------------------------------------------------------------
// Test DataGroup::createView() -- external, null pointer
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view_null)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  int * idata = ATK_NULLPTR;

  DataView * view = root->createView("null", idata)->apply(INT_ID, 0);

  EXPECT_TRUE(view->isApplied());
  EXPECT_FALSE(view->isAllocated());
  EXPECT_TRUE(view->isExternal());
  EXPECT_FALSE(view->isOpaque());
  EXPECT_EQ(view->getNumElements(), 0u);
  EXPECT_EQ(view->getTotalBytes(), 0u);

  void * ptr = view->getVoidPtr();
  EXPECT_EQ(ptr, ATK_NULLPTR);

  // getData will not work since the address is NULL
  //  int * idata_chk = view->getData();
  //EXPECT_EQ(idata_chk, ATK_NULLPTR);

  view->print();

  delete ds;
}
#if 0
//------------------------------------------------------------------------------
// Test DataGroup::save(), DataGroup::load() with described external views
// TODO - The save/load functionality needs to be fixed.
//------------------------------------------------------------------------------
TEST(sidre_external, save_load_external_view)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  const SidreLength len = 11;

  int * idata = new int[len];
  double * ddata = new double[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  DataView * iview = root->createView("idata", idata)->apply(INT_ID, len);
  DataView * dview = root->createView("ddata", ddata)->apply(DOUBLE_ID, len);
  EXPECT_EQ(root->getNumViews(), 2u);

//  iview->print();
//  dview->print();

  ds->getRoot()->save("out_sidre_external_save_restore_external_view",
                      "conduit");

//  ds->print();


  DataStore * ds2 = new DataStore();

  ds2->getRoot()->load("out_sidre_external_save_restore_external_view",
                       "conduit");

//  ds2->print();

  DataGroup * root2 = ds2->getRoot();

  EXPECT_EQ(root2->getNumViews(), 2u);

  int * idata_chk = iview->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = dview->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  delete ds;
  delete ds2;
  delete [] idata;
  delete [] ddata;
}
#endif
