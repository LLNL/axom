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

#include "common/CommonTypes.hpp"
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
  const int ndims = 1;
  SidreLength shape[] = { len };

  int * idata = new int[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
  }

  for (unsigned int i=0 ; i < 8 ; i++)
  {
    DataView * view = ATK_NULLPTR;

    switch (i)
    {
    case 0:
      view = root->createView("data0", INT_ID, len, idata);
      break;
    case 1:
      view = root->createView("data1", INT_ID, len)
             ->setExternalDataPtr(idata);
      break;
    case 2:
      view = root->createView("data2")
             ->setExternalDataPtr(INT_ID, len, idata);
      break;
    case 3:
      view = root->createView("data3", idata)
             ->apply(INT_ID, len);
      break;

    case 4:
      view = root->createView("data4", INT_ID, ndims, shape, idata);
      break;
    case 5:
      view = root->createView("data5", INT_ID, ndims, shape)
             ->setExternalDataPtr(idata);
      break;
    case 6:
      view = root->createView("data6")
             ->setExternalDataPtr(INT_ID, ndims, shape, idata);
      break;
    case 7:
      view = root->createView("data7", idata)
             ->apply(INT_ID, ndims, shape);
      break;
    }

    EXPECT_EQ(root->getNumViews(), i + 1);

    EXPECT_TRUE(view->isDescribed());
    EXPECT_TRUE(view->isAllocated());
    EXPECT_TRUE(view->isApplied());

    EXPECT_TRUE(view->isExternal());
    EXPECT_FALSE(view->isOpaque());

    EXPECT_EQ(view->getTypeID(), INT_ID);
    EXPECT_EQ(view->getNumElements(), len);

    view->print();

    int * idata_chk = view->getData();
    for (int ii = 0 ; ii < len ; ++ii)
    {
      EXPECT_EQ(idata_chk[ii], idata[ii]);
    }

  }
  delete ds;
  delete [] idata;
}

//------------------------------------------------------------------------------
// Test DataGroup::createView() -- external, null pointer
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view_null)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();
  int * idata = ATK_NULLPTR;

  for (int i=0 ; i < 2 ; i++)
  {
    DataView * view = ATK_NULLPTR;

    switch (i)
    {
    case 0:
      view = root->createView("null0", idata)->apply(INT_ID, 0);
      break;
    case 1:
      view = root->createView("null1", INT_ID, 0)->setExternalDataPtr(idata);
      break;
    case 2:
      //	  view = root->createView("null2", INT_ID, 0, idata);
      break;
    }

    EXPECT_TRUE(view->isDescribed());
    EXPECT_FALSE(view->isAllocated());
    EXPECT_FALSE(view->isApplied());

    EXPECT_TRUE(view->isEmpty());
    EXPECT_FALSE(view->isExternal());
    EXPECT_FALSE(view->isOpaque());

    EXPECT_EQ(view->getTypeID(), INT_ID);
    EXPECT_EQ(view->getNumElements(), 0u);
    EXPECT_EQ(view->getTotalBytes(), 0u);

    void * ptr = view->getVoidPtr();
    EXPECT_EQ(ptr, static_cast<void *>(ATK_NULLPTR));

    // getData will not work since the address is NULL
    //  int * idata_chk = view->getData();
    //EXPECT_EQ(idata_chk, ATK_NULLPTR);

    view->print();

  }
  delete ds;
}

//------------------------------------------------------------------------------
// Test EXTERNAL-> EMPTY
//------------------------------------------------------------------------------
TEST(sidre_external, transition_external_view_to_empty)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();
  const SidreLength len = 11;
  int idata[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
  }

  DataView * view = root->createView("data0", INT_ID, len)
                    ->setExternalDataPtr(idata);
  EXPECT_TRUE(view->isExternal());

  // Transition from EXTERNAL to EMPTY
  view->setExternalDataPtr(NULL);

  EXPECT_TRUE(view->isDescribed());
  EXPECT_FALSE(view->isAllocated());
  EXPECT_FALSE(view->isApplied());

  EXPECT_TRUE(view->isEmpty());
  EXPECT_FALSE(view->isExternal());
  EXPECT_FALSE(view->isOpaque());

  EXPECT_EQ(view->getTypeID(), INT_ID);
  EXPECT_EQ(view->getNumElements(), len);

  void * ptr = view->getVoidPtr();
  EXPECT_EQ(ptr, static_cast<void *>(ATK_NULLPTR));

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
