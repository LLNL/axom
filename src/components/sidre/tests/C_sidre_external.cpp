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
  SIDRE_datastore * ds = SIDRE_datastore_new();

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_databuffer * dbuff_0 = SIDRE_datastore_create_buffer(ds);
  SIDRE_databuffer * dbuff_1 = SIDRE_datastore_create_buffer(ds);
  SIDRE_databuffer * dbuff_2 = SIDRE_datastore_create_buffer(ds);

  SIDRE_databuffer_allocate_from_type(dbuff_0, SIDRE_DOUBLE_ID, len);
  SIDRE_databuffer_declare(dbuff_1, SIDRE_INT_ID, len);
  SIDRE_databuffer_set_external_data(dbuff_1, idata);
  SIDRE_databuffer_declare(dbuff_2, SIDRE_DOUBLE_ID, len);
  SIDRE_databuffer_set_external_data(dbuff_2, ddata);

  EXPECT_EQ(SIDRE_databuffer_is_external(dbuff_0), false);
  EXPECT_EQ(SIDRE_databuffer_is_external(dbuff_1), true);
  EXPECT_EQ(SIDRE_databuffer_is_external(dbuff_2), true);

  EXPECT_EQ(SIDRE_databuffer_get_total_bytes(dbuff_0), sizeof(double)*len);
  EXPECT_EQ(SIDRE_databuffer_get_total_bytes(dbuff_1), sizeof(int)*len);
  EXPECT_EQ(SIDRE_databuffer_get_total_bytes(dbuff_2), sizeof(double)*len);

  SIDRE_datastore_print(ds);

  SIDRE_datastore_delete(ds);
  free(idata);
  free(ddata);
}

//------------------------------------------------------------------------------
// Test DataGroup::create_external_view()
//------------------------------------------------------------------------------
TEST(C_sidre_external, create_external_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_dataview * iview = SIDRE_datagroup_create_external_view(root, "idata", idata, SIDRE_INT_ID, len);
  SIDRE_dataview * dview = SIDRE_datagroup_create_external_view(root, "ddata", ddata, SIDRE_DOUBLE_ID, len);
  EXPECT_EQ(SIDRE_datagroup_get_num_views(root), 2u);

  SIDRE_dataview_print(iview);
  SIDRE_dataview_print(dview);

  int * idata_chk = (int *) SIDRE_dataview_get_data_pointer(iview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = (double *) SIDRE_dataview_get_data_pointer(dview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_datastore_delete(ds);
  free(idata);
  free(ddata);
}

//------------------------------------------------------------------------------
// Test DataGroup::save(), DataGroup::load() with external buffers
//------------------------------------------------------------------------------
TEST(C_sidre_external, save_load_external_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  const int len = 11;

  int * idata = (int *) malloc(sizeof(int) * len);
  double * ddata = (double *) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_dataview * iview = SIDRE_datagroup_create_external_view(root, "idata", idata, SIDRE_INT_ID, len);
  SIDRE_dataview * dview = SIDRE_datagroup_create_external_view(root, "ddata", ddata, SIDRE_DOUBLE_ID, len);
  EXPECT_EQ(SIDRE_datagroup_get_num_views(root), 2u);
  SIDRE_databuffer * tmpbuf;
  tmpbuf = SIDRE_dataview_get_buffer(iview);
  EXPECT_EQ(SIDRE_databuffer_is_external(tmpbuf), true);
  tmpbuf = SIDRE_dataview_get_buffer(dview);
  EXPECT_EQ(SIDRE_databuffer_is_external(tmpbuf), true);

  SIDRE_dataview_print(iview);
  SIDRE_dataview_print(dview);

  SIDRE_datagroup_save(root, "out_sidre_external_save_restore_external_view", "conduit");

  SIDRE_datastore_print(ds);


  SIDRE_datastore * ds2 = SIDRE_datastore_new();
  SIDRE_datagroup * root2 = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup_load(root, "out_sidre_external_save_restore_external_view","conduit");

  SIDRE_datastore_print(ds2);

  SIDRE_dataview * iview2 = SIDRE_datagroup_get_view_from_name(root2, "idata");
  SIDRE_dataview * dview2 = SIDRE_datagroup_get_view_from_name(root2, "ddata");

  EXPECT_EQ(SIDRE_datagroup_get_num_views(root2), 2u);
  tmpbuf = SIDRE_dataview_get_buffer(iview2);
  EXPECT_EQ(SIDRE_databuffer_is_external(tmpbuf), false);
  tmpbuf = SIDRE_dataview_get_buffer(dview2);
  EXPECT_EQ(SIDRE_databuffer_is_external(tmpbuf), false);

  int * idata_chk = (int *) SIDRE_dataview_get_data_pointer(iview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double * ddata_chk = (double *) SIDRE_dataview_get_data_pointer(dview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);
  free(idata);
  free(ddata);
}
