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
using asctoolkit::sidre::DataType;

//------------------------------------------------------------------------------
// Test DataGroup::createView() -- external
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view)
{
  DataStore * ds   = new DataStore();
  DataGroup * root = ds->getRoot();

  const int len = 11;

  int * idata = new int[len];
  double * ddata = new double[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  (void) root->createView("idata", idata)->apply(DataType::c_int(len));
  (void) root->createView("ddata", ddata)->apply(DataType::c_double(len));
  EXPECT_EQ(root->getNumViews(), 2u);

  root->getView("idata")->print();
  root->getView("ddata")->print();

  int * idata_chk = root->getView("idata")->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = root->getView("ddata")->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  delete ds;
  delete [] idata;
  delete [] ddata;
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

  const int len = 11;

  int * idata = new int[len];
  double * ddata = new double[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  (void) root->createView("idata", idata)->apply(DataType::c_int(len));
  (void) root->createView("ddata", ddata)->apply(DataType::c_double(len));
  EXPECT_EQ(root->getNumViews(), 2u);

//  root->getView("idata")->print();
//  root->getView("ddata")->print();

  ds->getRoot()->save("out_sidre_external_save_restore_external_view",
                      "conduit");

//  ds->print();


  DataStore * ds2 = new DataStore();

  ds2->getRoot()->load("out_sidre_external_save_restore_external_view",
                       "conduit");

//  ds2->print();

  DataGroup * root2 = ds2->getRoot();

  EXPECT_EQ(root2->getNumViews(), 2u);

  int * idata_chk = root2->getView("idata")->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = root2->getView("ddata")->getData();
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
