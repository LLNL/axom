// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <vector>

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/sidre/core/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::DOUBLE_ID;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::View;

//------------------------------------------------------------------------------
// Test Group::createView() -- external
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const IndexType len = 11;
  const int ndims = 1;
  IndexType shape[] = {len};

  int* idata = new int[len];

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
  }

  for(unsigned int i = 0; i < 8; i++)
  {
    View* view = nullptr;

    switch(i)
    {
    case 0:
      view = root->createView("data0", INT_ID, len, idata);
      break;
    case 1:
      view = root->createView("data1", INT_ID, len)->setExternalDataPtr(idata);
      break;
    case 2:
      view = root->createView("data2")->setExternalDataPtr(INT_ID, len, idata);
      break;
    case 3:
      view = root->createView("data3", idata)->apply(INT_ID, len);
      break;

    case 4:
      view = root->createView("data4", INT_ID, ndims, shape, idata);
      break;
    case 5:
      view =
        root->createView("data5", INT_ID, ndims, shape)->setExternalDataPtr(idata);
      break;
    case 6:
      view =
        root->createView("data6")->setExternalDataPtr(INT_ID, ndims, shape, idata);
      break;
    case 7:
      view = root->createView("data7", idata)->apply(INT_ID, ndims, shape);
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

    EXPECT_EQ(idata, view->getVoidPtr());

    view->print();

    int* idata_chk = view->getData();
    for(int ii = 0; ii < len; ++ii)
    {
      EXPECT_EQ(idata_chk[ii], idata[ii]);
    }
  }
  delete ds;
  delete[] idata;
}

//------------------------------------------------------------------------------
// Test Group::createView() -- external, null pointer
//------------------------------------------------------------------------------
TEST(sidre_external, create_external_view_null)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  int* idata = nullptr;

  for(int i = 0; i < 2; i++)
  {
    View* view = nullptr;

    switch(i)
    {
    case 0:
      view = root->createView("null0", idata)->apply(INT_ID, 0);
      break;
    case 1:
      view = root->createView("null1", INT_ID, 0)->setExternalDataPtr(idata);
      break;
    case 2:
      //    view = root->createView("null2", INT_ID, 0, idata);
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

    void* ptr = view->getVoidPtr();
    EXPECT_EQ(static_cast<void*>(nullptr), ptr);

    // getData will not work since the address is NULL
    //  int * idata_chk = view->getData();
    //EXPECT_EQ(idata_chk, nullptr);

    view->print();
  }
  delete ds;
}

//------------------------------------------------------------------------------
// Test EXTERNAL-> EMPTY
//------------------------------------------------------------------------------
TEST(sidre_external, transition_external_view_to_empty)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  const IndexType len = 11;
  int idata[len];

  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii;
  }

  View* view = root->createView("data0", INT_ID, len)->setExternalDataPtr(idata);
  EXPECT_TRUE(view->isExternal());
  EXPECT_EQ(idata, view->getVoidPtr());

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

  void* ptr = view->getVoidPtr();
  EXPECT_EQ(static_cast<void*>(nullptr), ptr);

  view->print();

  delete ds;
}

TEST(sidre_external, verify_external_layout)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  const int SZ = 11;
  int extData[SZ];

  SLIC_INFO(
    "Tests that the hierarchy generated by Group::"
    << "createExternalLayout() consists only of groups and views"
    << " associated with external pointers (described or undescribed).");

  // Create some internal views
  root->createViewAndAllocate("int/desc/bufferview", axom::sidre::INT_ID, SZ);
  root->createViewScalar("int/scalar/scalarview", SZ);
  root->createViewString("int/string/stringview", "A string");

  // Sanity check on internal views
  EXPECT_FALSE(root->getView("int/desc/bufferview")->isExternal());
  EXPECT_TRUE(root->getView("int/desc/bufferview")->isDescribed());
  EXPECT_FALSE(root->getView("int/scalar/scalarview")->isExternal());
  EXPECT_TRUE(root->getView("int/scalar/scalarview")->isDescribed());
  EXPECT_FALSE(root->getView("int/string/stringview")->isExternal());
  EXPECT_TRUE(root->getView("int/string/stringview")->isDescribed());

  // Check that the generated external layout is empty
  {
    axom::sidre::Node emptyNode;
    root->createExternalLayout(emptyNode);

    SLIC_INFO("External node's json before adding external data."
              << " (should be empty): \n\t" << emptyNode.to_json());

    EXPECT_EQ(0, emptyNode.number_of_children());
  }

  // Create some external views
  root->createView("ext/desc/external_desc", axom::sidre::INT_ID, SZ, extData);
  root->createView("ext/undesc/external_opaque")->setExternalDataPtr(extData);

  // Sanity check on the external views
  EXPECT_TRUE(root->getView("ext/desc/external_desc")->isExternal());
  EXPECT_TRUE(root->getView("ext/desc/external_desc")->isDescribed());

  EXPECT_TRUE(root->getView("ext/undesc/external_opaque")->isExternal());
  EXPECT_FALSE(root->getView("ext/undesc/external_opaque")->isDescribed());

  // Initialize the data so we can test it later
  int* bufData = root->getView("int/desc/bufferview")->getData();
  for(int i = 0; i < SZ; ++i)
  {
    extData[i] = i;
    bufData[i] = 100 + i;
  }

  // Check that the generated external layout matches expectation
  {
    /// Copy the external layout into a conduit node and test the layout
    axom::sidre::Node node;
    root->createExternalLayout(node);

    SLIC_INFO("External node's json: \n\t" << node.to_json());

    EXPECT_TRUE(node.number_of_children() == 1);
    EXPECT_TRUE(node.has_path("ext"));
    EXPECT_TRUE(node["ext"].number_of_children() == 2);

    // Described external views are present and we can access the data
    EXPECT_TRUE(node.has_path("ext/desc"));
    EXPECT_TRUE(node.has_path("ext/desc/external_desc"));
    int* extLayoutData = node["ext/desc/external_desc"].value();
    for(int i = 0; i < SZ; ++i)
    {
      EXPECT_EQ(extData[i], extLayoutData[i]);
    }

    // Opaque views are not present, but their containing groups are present
    EXPECT_TRUE(node.has_path("ext/undesc"));
    EXPECT_FALSE(node.has_path("ext/undesc/external_opaque"));

    // Buffer and scalar views are not present in the external layout
    EXPECT_FALSE(node.has_path("int"));
  }
}

#if 0
//------------------------------------------------------------------------------
// Test Group::save(), Group::load() with described external views
// TODO - The save/load functionality needs to be fixed.
//------------------------------------------------------------------------------
TEST(sidre_external, save_load_external_view)
{
  DataStore* ds   = new DataStore();
  Group* root = ds->getRoot();

  const IndexType len = 11;

  int* idata = new int[len];
  double* ddata = new double[len];

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  View* iview = root->createView("idata", idata)->apply(INT_ID, len);
  View* dview = root->createView("ddata", ddata)->apply(DOUBLE_ID, len);
  EXPECT_EQ(root->getNumViews(), 2u);

//  iview->print();
//  dview->print();

  ds->getRoot()->save("out_sidre_external_save_restore_external_view",
                      "conduit");

//  ds->print();


  DataStore* ds2 = new DataStore();

  ds2->getRoot()->load("out_sidre_external_save_restore_external_view",
                       "conduit");

//  ds2->print();

  Group* root2 = ds2->getRoot();

  EXPECT_EQ(root2->getNumViews(), 2u);

  int* idata_chk = iview->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = dview->getData();
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
