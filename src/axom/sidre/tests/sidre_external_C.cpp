/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

#include <stdlib.h>

//------------------------------------------------------------------------------
// Test Group::create_external_view()
//------------------------------------------------------------------------------
TEST(C_sidre_external, create_external_view)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf;
  SIDRE_view iview_buf, dview_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  const int len = 11;

  int* idata = (int*) malloc(sizeof(int) * len);
  double* ddata = (double*) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_view* iview =
    SIDRE_group_create_view_external(root, "idata", idata, &iview_buf);
  SIDRE_view_apply_type_nelems(iview, SIDRE_INT_ID, len);
  SIDRE_view* dview =
    SIDRE_group_create_view_external(root, "ddata", ddata, &dview_buf);
  SIDRE_view_apply_type_nelems(dview, SIDRE_DOUBLE_ID, len);
  EXPECT_EQ(SIDRE_group_get_num_views(root), 2u);

  SIDRE_view_print(iview);
  SIDRE_view_print(dview);

  int* idata_chk = (int*) SIDRE_view_get_void_ptr(iview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = (double*) SIDRE_view_get_void_ptr(dview);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_datastore_delete(ds);
  free(idata);
  free(ddata);
}

#if 0
//------------------------------------------------------------------------------
// Test Group::save(), Group::load() with external buffers
//------------------------------------------------------------------------------
TEST(C_sidre_external, save_load_external_view)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  const int len = 11;

  int* idata = (int*) malloc(sizeof(int) * len);
  double* ddata = (double*) malloc(sizeof(double) * len);

  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii;
    ddata[ii] = idata[ii] * 2.0;
  }

  SIDRE_view* iview =
    SIDRE_group_create_view_external(root, "idata", idata);
  SIDRE_view_apply_type_nelems(iview, SIDRE_INT_ID, len);
  SIDRE_view* dview =
    SIDRE_group_create_view_external(root, "ddata", ddata);
  SIDRE_view_apply_type_nelems(dview, SIDRE_DOUBLE_ID, len);

  EXPECT_EQ(SIDRE_group_get_num_views(root), 2u);

  SIDRE_view_print(iview);
  SIDRE_view_print(dview);

  SIDRE_group_save(root, "out_sidre_external_save_restore_external_view",
                   "conduit");

  SIDRE_datastore_print(ds);


  SIDRE_datastore* ds2 = SIDRE_datastore_new();
  SIDRE_group* root2 = SIDRE_datastore_get_root(ds);

  SIDRE_group_load(root, "out_sidre_external_save_restore_external_view",
                   "conduit");

  SIDRE_datastore_print(ds2);

  SIDRE_view* iview2 = SIDRE_group_get_view_from_name(root2, "idata");
  SIDRE_view* dview2 = SIDRE_group_get_view_from_name(root2, "ddata");

  EXPECT_EQ(SIDRE_group_get_num_views(root2), 2u);

  int* idata_chk = (int*) SIDRE_view_get_void_ptr(iview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(idata_chk[ii], idata[ii]);
  }

  double* ddata_chk = (double*) SIDRE_view_get_void_ptr(dview2);
  for (int ii = 0 ; ii < len ; ++ii)
  {
    EXPECT_EQ(ddata_chk[ii], ddata[ii]);
  }

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);
  free(idata);
  free(ddata);
}

#endif
