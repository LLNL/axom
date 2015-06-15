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

#include "sidre/sidre.h"

#include <stdlib.h>

//------------------------------------------------------------------------------
// Test DataBuffer::declareExternal()
//------------------------------------------------------------------------------

TEST(C_sidre_external, declare_external_buffer)
{
  ATK_datastore * ds = ATK_datastore_new();

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  ATK_databuffer * dbuff_0 = ATK_datastore_create_buffer(ds);
  ATK_databuffer * dbuff_1 = ATK_datastore_create_buffer(ds);
  ATK_databuffer * dbuff_2 = ATK_datastore_create_buffer(ds);

  ATK_databuffer_allocate_from_type(dbuff_0, ATK_C_DOUBLE_T, len);
  ATK_databuffer_declare_external(dbuff_1, idata, ATK_C_INT_T, len);
  ATK_databuffer_declare_external(dbuff_2, ddata, ATK_C_DOUBLE_T, len);

  EXPECT_EQ(ATK_databuffer_is_external(dbuff_0), false);
  EXPECT_EQ(ATK_databuffer_is_external(dbuff_1), true);
  EXPECT_EQ(ATK_databuffer_is_external(dbuff_2), true);

  EXPECT_EQ(ATK_databuffer_get_total_bytes(dbuff_0), sizeof(double)*len);
  EXPECT_EQ(ATK_databuffer_get_total_bytes(dbuff_1), sizeof(int)*len);
  EXPECT_EQ(ATK_databuffer_get_total_bytes(dbuff_2), sizeof(double)*len);

  ATK_datastore_print(ds);

  ATK_datastore_delete(ds);
  free(idata);
  free(ddata);
}

//------------------------------------------------------------------------------
// Test DataGroup::create_external_view()
//------------------------------------------------------------------------------
TEST(C_sidre_external, create_external_view)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  ATK_dataview *iview = ATK_datagroup_create_external_view(root, "idata", idata, ATK_C_INT_T, len);
  ATK_dataview *dview = ATK_datagroup_create_external_view(root, "ddata", ddata, ATK_C_DOUBLE_T, len);
  EXPECT_EQ(ATK_datagroup_get_num_views(root), 2u);

#ifdef XXX
  root->getView("idata")->getNode().print_detailed();
  root->getView("ddata")->getNode().print_detailed();
#endif

  int * idata_chk = (int *) ATK_dataview_get_data_buffer(iview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * idata_chk = (double *) ATK_dataview_get_data_buffer(dview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  ATK_datastore_delete(ds);
  free(idata);
  free(ddata);
}

//------------------------------------------------------------------------------
// Test DataGroup::save(), DataGroup::load() with external buffers
//------------------------------------------------------------------------------
TEST(C_sidre_external, save_load_external_view)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  ATK_dataview *iview = ATK_datagroup_create_external_view(root, "idata", idata, ATK_C_INT_T, len);
  ATK_dataview *dview = ATK_datagroup_create_external_view(root, "ddata", ddata, ATK_C_DOUBLE_T, len);
  EXPECT_EQ(ATK_datagroup_get_num_views(root), 2u);
  ATK_databuffer *tmpbuf;
  tmpbuf = ATK_dataview_get_buffer(iview)
  EXPECT_EQ(ATK_databuffer_is_external(), true);
  tmpbuf = ATK_dataview_get_buffer(dview)
  EXPECT_EQ(ATK_databuffer_is_external(), true);

#ifdef XXX
  ATK_dataview *iview = ATK_datagroup_get_view(root, "idata");
  ATK_dataview *dview = ATK_datagroup_get_view(root, "idata");
  root->getView("idata")->getNode().print_detailed();
  root->getView("ddata")->getNode().print_detailed();
#endif

  ATK_datagroup_save(root, "out_sidre_external_save_restore_external_view", "conduit");

  ATK_databuffer_print(ds);


  ATK_datastore * ds2 = ATK_datastore_new();
  ATK_datagroup * root2 = ATK_datastore_get_root(ds);

  ATK_datagroup_load(root, "out_sidre_external_save_restore_external_view","conduit");

  ATK_datagroup_print(ds2);

  ATK_dataview iview2 = ATK_datagroup_getView(root2, "idata");
  ATK_dataview dview2 = ATK_datagroup_getView(root2, "ddata");

  EXPECT_EQ(ATK_datagroup_get_num_views(root2), 2u);
  ATK_databuffer tmpbuf;
  tmpbuf = ATK_dataview_get_buffer(iview2);
  EXPECT_EQ(ATK_databuffer_is_external(tmpbuf), false);
  tmpbuf = ATK_dataview_get_buffer(dview2);
  EXPECT_EQ(ATK_databuffer_is_external(tmpbuf), false);

  int * idata_chk = (int *) ATK_dataview_get_data_buffer(iview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * idata_chk = (double *) ATK_dataview_get_data_buffer(dview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  ATK_datastore_delete(ds);
  ATK_datastore_delete(ds2);
  free(idata);
  free(ddata);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
   int result = 0;

   ::testing::InitGoogleTest(&argc, argv);

   UnitTestLogger logger;  // create & initialize test logger,
                       // finalized when exiting main scope

   result = RUN_ALL_TESTS();

   return result;
}
