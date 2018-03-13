/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#include "sidre/sidre.h"
#include "slic/wrapSLIC.h"

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_name)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* group = SIDRE_group_create_group(root, "test");

  //    EXPECT_TRUE(group->getName() == std::string("test") );
  EXPECT_TRUE(strcmp(SIDRE_group_get_name(group), "test") == 0);

  SIDRE_group* group2 = SIDRE_group_get_group_from_name(root, "foo");
  EXPECT_TRUE(group2 == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_parent)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_group* child = SIDRE_group_create_group(parent, "child");

  EXPECT_TRUE( SIDRE_group_get_parent(child) == parent );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_datastore)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* group = SIDRE_group_create_group(root, "parent");

  EXPECT_TRUE( SIDRE_group_get_data_store(group) == ds );

  SIDRE_datastore const* const_ds = SIDRE_group_get_data_store(group);
  EXPECT_TRUE( const_ds == ds );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_group* child = SIDRE_group_create_group(parent, "child");
  EXPECT_TRUE( SIDRE_group_get_parent(child) == parent );

  EXPECT_TRUE( SIDRE_group_get_group_from_name(parent, "child") == child);
  EXPECT_TRUE( SIDRE_group_get_group_from_index(parent, 0) == child);

  // check error condition
  EXPECT_TRUE( SIDRE_group_get_group_from_name(parent,
                                               "non-existant group") == NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");

  SIDRE_view* view = SIDRE_group_create_view_empty(parent, "view");

  EXPECT_TRUE( SIDRE_group_get_view_from_name(parent, "view") == view );
  EXPECT_TRUE( SIDRE_group_get_view_from_index(parent, 0) == view );

  // check error condition
  EXPECT_TRUE( SIDRE_group_get_view_from_name(parent,
                                              "non-existant view") ==
               NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view_names_and_indicies)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_view* view1 = SIDRE_group_create_view_empty(parent, "view1");
  SIDRE_view* view2 = SIDRE_group_create_view_empty(parent, "view2");

  EXPECT_EQ(SIDRE_group_get_num_views(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_view_index(parent, "view1");
  SIDRE_IndexType idx2 = SIDRE_group_get_view_index(parent, "view2");

  const char* name1 = SIDRE_group_get_view_name(parent, idx1);
  const char* name2 = SIDRE_group_get_view_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "view1") == 0);
  EXPECT_TRUE(strcmp(SIDRE_view_get_name(view1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "view2") == 0);
  EXPECT_TRUE(strcmp(SIDRE_view_get_name(view2), name2) == 0);

  // check error conditions
  SIDRE_IndexType idx3 = SIDRE_group_get_view_index(parent, "view3");
  EXPECT_TRUE(idx3 == SIDRE_InvalidIndex);

  const char* name3 = SIDRE_group_get_view_name(parent, idx3);
  EXPECT_TRUE(name3 == NULL);
  EXPECT_FALSE(SIDRE_name_is_valid(name3));

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getFirstValidGroupIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_group_index)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_group* group1 = SIDRE_group_create_group(parent, "group1");
  SIDRE_group* group2 = SIDRE_group_create_group(parent, "group2");
  EXPECT_EQ(SIDRE_group_get_num_groups(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_first_valid_group_index(parent);
  SIDRE_IndexType idx2 = SIDRE_group_get_next_valid_group_index(parent, idx1);
  SIDRE_IndexType idx3 = SIDRE_group_get_next_valid_group_index(parent, idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(SIDRE_InvalidIndex, idx3);

  SIDRE_group* group1out = SIDRE_group_get_group_from_index(parent, idx1);
  SIDRE_group* group2out = SIDRE_group_get_group_from_index(parent, idx2);
  EXPECT_EQ(group1, group1out);
  EXPECT_EQ(group2, group2out);

  // check error conditions
  SIDRE_group* emptyGroup =
    SIDRE_group_create_group(root, "emptyGroup");

  SIDRE_IndexType badidx1 = SIDRE_group_get_first_valid_group_index(
    emptyGroup);
  SIDRE_IndexType badidx2 = SIDRE_group_get_next_valid_group_index(
    emptyGroup, badidx1);

  EXPECT_EQ(SIDRE_InvalidIndex, badidx1);
  EXPECT_EQ(SIDRE_InvalidIndex, badidx2);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_view_index)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_view* view1 = SIDRE_group_create_view_empty(parent, "view1");
  SIDRE_view* view2 = SIDRE_group_create_view_empty(parent, "view2");
  EXPECT_EQ(SIDRE_group_get_num_views(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_first_valid_view_index(parent);
  SIDRE_IndexType idx2 = SIDRE_group_get_next_valid_view_index(parent, idx1);
  SIDRE_IndexType idx3 = SIDRE_group_get_next_valid_view_index(parent, idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(SIDRE_InvalidIndex, idx3);

  SIDRE_view* view1out = SIDRE_group_get_view_from_index(parent, idx1);
  SIDRE_view* view2out = SIDRE_group_get_view_from_index(parent, idx2);
  EXPECT_EQ(view1, view1out);
  EXPECT_EQ(view2, view2out);

  // check error conditions
  SIDRE_group* emptyGroup =
    SIDRE_group_create_group(root, "emptyGroup");

  SIDRE_IndexType badidx1 = SIDRE_group_get_first_valid_view_index(
    emptyGroup);
  SIDRE_IndexType badidx2 = SIDRE_group_get_next_valid_view_index(
    emptyGroup, badidx1);

  EXPECT_EQ(SIDRE_InvalidIndex, badidx1);
  EXPECT_EQ(SIDRE_InvalidIndex, badidx2);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group_name_index)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent");
  SIDRE_group* group1 = SIDRE_group_create_group(parent, "group1");
  SIDRE_group* group2 = SIDRE_group_create_group(parent, "group2");

  EXPECT_EQ(SIDRE_group_get_num_groups(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_group_index(parent, "group1");
  SIDRE_IndexType idx2 = SIDRE_group_get_group_index(parent, "group2");

  const char* name1 = SIDRE_group_get_group_name(parent, idx1);
  const char* name2 = SIDRE_group_get_group_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "group1") == 0);
  EXPECT_TRUE(strcmp(SIDRE_group_get_name(group1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "group2") == 0);
  EXPECT_TRUE(strcmp(SIDRE_group_get_name(group2), name2) == 0);

  // check error conditions
  SIDRE_IndexType idx3 = SIDRE_group_get_group_index(parent, "group3");
  EXPECT_TRUE(idx3 == SIDRE_InvalidIndex);

  const char* name3 = SIDRE_group_get_group_name(parent, idx3);
  EXPECT_TRUE(name3 == NULL);
  EXPECT_TRUE(SIDRE_name_is_valid(name3) == 0);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// createView()
// createViewAndAllocate()
// destroyView()
// destroyViewAndData()
// hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_view)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* group = SIDRE_group_create_group(root, "parent");

  SIDRE_view* view = SIDRE_group_create_view_empty(group, "view");
  EXPECT_TRUE( SIDRE_group_get_parent(group) == root );
  EXPECT_FALSE( SIDRE_view_has_buffer(view) );

  EXPECT_TRUE( SIDRE_group_has_view(group, "view") );
  // try creating view again, should be a no-op.
  EXPECT_TRUE( SIDRE_group_create_view_empty(group, "view") == NULL );

  SIDRE_group_destroy_view(group, "view");
  // destroy already destroyed group.  Should be a no-op, not a failure
  SIDRE_group_destroy_view(group, "view");

  EXPECT_FALSE( SIDRE_group_has_view(group, "view") );

  // try api call that specifies specific type and length
  SIDRE_group_create_view_and_allocate_nelems( group,
                                               "viewWithLength1",
                                               SIDRE_FLOAT_ID,
                                               50);

  // error condition check - try again with duplicate name, should be a no-op
  //XXX EXPECT_TRUE( SIDRE_group_create_view_and_allocate( group,
  // "viewWithLength1", SIDRE_FLOAT64_ID, 50) );
  SIDRE_group_destroy_view_and_data_name( group, "viewWithLength1");
  EXPECT_FALSE( SIDRE_group_has_view( group, "viewWithLength1") );

  EXPECT_TRUE( SIDRE_group_create_view_and_allocate_nelems( group,
                                                            "viewWithLengthBadLen",
                                                            SIDRE_FLOAT64_ID,
                                                            -1) == NULL );

  // try api call that specifies data type in another way
  SIDRE_group_create_view_and_allocate_nelems( group, "viewWithLength2",
                                               SIDRE_FLOAT64_ID, 50 );
  EXPECT_TRUE( SIDRE_group_create_view_and_allocate_nelems( group,
                                                            "viewWithLength2",
                                                            SIDRE_FLOAT64_ID,
                                                            50 ) == NULL );
  // destroy this view using index
  SIDRE_group_destroy_view_and_data_index( group, SIDRE_group_get_first_valid_view_index(
                                             group) );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_group)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* group = SIDRE_group_create_group(root, "group");
  EXPECT_TRUE( SIDRE_group_get_parent(group) == root );

  EXPECT_TRUE( SIDRE_group_has_group(root, "group") );

  SIDRE_group_destroy_group_name(root, "group");
  EXPECT_FALSE( SIDRE_group_has_group(root, "group") );

  SIDRE_group* group2 = SIDRE_group_create_group(root, "group2");
  // shut up compiler about unused variable
  (void)group2;
  SIDRE_group_destroy_group_index( root, SIDRE_group_get_first_valid_group_index(
                                     root) );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");
  SIDRE_group_create_view_empty(flds, "a");

  EXPECT_TRUE(SIDRE_group_has_view(flds, "a"));

  // attempt to create duplicate group name

  SIDRE_group* badGroup = SIDRE_group_create_group(root, "fields");
  EXPECT_TRUE( badGroup == NULL );

  // check error condition
  // attempt to create duplicate view name.
  EXPECT_TRUE(SIDRE_group_create_view_empty(flds, "a") == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
// Restore this after copy_move is working. ATK-667
#if 0
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");

  SIDRE_group_create_view_scalar_int(flds, "i0", 1);
  SIDRE_group_create_view_scalar_float(flds, "f0", 100.0);
  SIDRE_group_create_view_scalar_double(flds, "d0", 3000.0);

  EXPECT_TRUE(SIDRE_group_has_view(flds, "i0"));
  EXPECT_TRUE(SIDRE_group_has_view(flds, "f0"));
  EXPECT_TRUE(SIDRE_group_has_view(flds, "d0"));

  // test moving a view form feds7 to sub
  // flds->createGroup("sub")->moveView(flds->getView("d0"));
  SIDRE_group* sub = SIDRE_group_create_group(flds, "sub");
  SIDRE_group_move_view(sub,
                        SIDRE_group_get_view_from_name(flds, "d0"));
  SIDRE_group_print(flds);
  EXPECT_FALSE(SIDRE_group_has_view(flds, "d0"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "sub"));
  EXPECT_TRUE(SIDRE_group_has_view(sub, "d0"));

  // check the data value
  // double *d0_data =  flds->getGroup("sub")
  //                        ->getView("d0")
  //                        ->getNode().as_double_ptr();
  double* d0_data;
  {
    SIDRE_view* tmpview = SIDRE_group_get_view_from_name(sub, "d0");
    SIDRE_buffer* tmpbuf = SIDRE_view_get_buffer(tmpview);
    d0_data = (double*) SIDRE_buffer_get_void_ptr(tmpbuf);
  }
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds top sub
  SIDRE_group_copy_view(sub,
                        SIDRE_group_get_view_from_name(flds, "i0"));

  SIDRE_group_print(flds);

  EXPECT_TRUE(SIDRE_group_has_view(flds, "i0"));
  EXPECT_TRUE(SIDRE_group_has_view(sub, "i0"));

#endif

#ifdef XXX
  // we expect the actual data  pointers to be the same
  EXPECT_EQ(SIDRE_group_get_view(flds, "i0")->getNode().data_pointer(),
            SIDRE_group_get_group_from_name("sub")->get_view(
              "i0")->getNode().data_pointer());
#endif

//  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a");
  SIDRE_group* gb = SIDRE_group_create_group(flds, "b");
  SIDRE_group* gc = SIDRE_group_create_group(flds, "c");

  SIDRE_view* i0_view =
    SIDRE_group_create_view_and_allocate_nelems(ga, "i0",
                                                SIDRE_INT_ID, 1);
  SIDRE_view* f0_view =
    SIDRE_group_create_view_and_allocate_nelems(gb, "f0",
                                                SIDRE_FLOAT_ID, 1);
  SIDRE_view* d0_view =
    SIDRE_group_create_view_and_allocate_nelems(gc, "d0",
                                                SIDRE_DOUBLE_ID, 1);

  SIDRE_view_set_scalar_int(i0_view, 1);
  SIDRE_view_set_scalar_float(f0_view, 100.0);
  SIDRE_view_set_scalar_double(d0_view, 3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

  //move "b" to a child of "sub"
  SIDRE_group* sub = SIDRE_group_create_group(flds, "sub");
  SIDRE_group_move_group(sub, gb);

  SIDRE_group_print(flds);

  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "sub"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

  SIDRE_group* tmpgrp = SIDRE_group_get_group_from_name(flds, "sub");
  EXPECT_EQ(SIDRE_group_get_group_from_name(tmpgrp, "b"), gb);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  SIDRE_datastore* const ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* const grp = SIDRE_group_create_group(root, "grp");

  const char* viewName1 = "viewBuffer1";
  const char* viewName2 = "viewBuffer2";

  // XXX const
  SIDRE_view* view1 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName1,
                                                SIDRE_INT_ID, 1);
  SIDRE_view* view2 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName2,
                                                SIDRE_FLOAT_ID, 1);

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName1));
  EXPECT_EQ( SIDRE_group_get_view_from_name(grp, viewName1), view1 );

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName2));
  EXPECT_EQ( SIDRE_group_get_view_from_name(grp, viewName2), view2 );

  SIDRE_buffer* tmpbuf = SIDRE_view_get_buffer(view1);
  SIDRE_IndexType bufferId1 = SIDRE_buffer_get_index(tmpbuf);

  SIDRE_group_destroy_view_and_data_name(grp, viewName1);


  EXPECT_FALSE(SIDRE_group_has_view(grp, viewName1));
  EXPECT_EQ(SIDRE_datastore_get_num_buffers(ds), 1u);

  SIDRE_buffer* buffer1 = SIDRE_datastore_get_buffer(ds, bufferId1);
  bool buffValid = true;
  if( buffer1 == NULL )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_alloc_view_and_buffer)
{
  SIDRE_datastore* const ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* const grp = SIDRE_group_create_group(root, "grp");

  const char* viewName1 = "viewBuffer1";

  // use create + alloc convenience methods
  // this one is the DataType & method
  SIDRE_view* const view1 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName1,
                                                SIDRE_INT_ID, 10);

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName1));
  EXPECT_EQ( SIDRE_group_get_view_from_name(grp, viewName1), view1 );


#ifdef XXX
  int* v1_vals = (int*) SIDRE_view_get_void_ptr(view1);
  double* v2_vals = (double*)  SIDRE_view_get_void_ptr(view1);

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
    v2_vals[i] = i * 3.1415;
  }
#endif

  EXPECT_EQ(SIDRE_view_get_num_elements(view1), 10u);
  EXPECT_EQ(SIDRE_view_get_total_bytes(view1), 10 * sizeof(int));

  SIDRE_group_destroy_view_and_data_name(grp, viewName1);

  SIDRE_datastore_delete(ds);
}

#ifdef XXX
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_view_of_buffer_with_datatype)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  // use create + alloc convenience methods
  // this one is the DataType & method
  SIDRE_view* base =
    SIDRE_group_create_view_and_allocate_nelems(root, "base",
                                                SIDRE_INT_ID, 10);
#ifdef XXX
  int* base_vals = (int*) SIDRE_view_get_void_ptr(base);
  for(int i=0 ; i<10 ; i++)
  {
    if(i < 5)
    {
      base_vals[i] = 10;
    }
    else
    {
      base_vals[i] = 20;
    }
  }
#endif

  SIDRE_buffer* base_buff = SIDRE_view_get_buffer(base);
  // create two views into this buffer
  // view for the first 5 values
  SIDRE_group_create_view(root, "sub_a", base_buff, SIDRE_C_INT_T, 5);
  // view for the second 5 values

  int* sub_a_vals = (int*) SIDRE_view_get_void_ptr(SIDRE_group_get_view_from_name(
                                                     root,
                                                     "sub_a"));

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
  }

  SIDRE_datastore_delete(ds);
}
#endif

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_simple)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a");

  SIDRE_group_create_view_scalar_int(ga, "i0", 1);

  EXPECT_TRUE(SIDRE_group_has_group(root, "fields"));
  EXPECT_TRUE(SIDRE_group_has_group(SIDRE_group_get_group_from_name(root,
                                                                    "fields"),
                                    "a"));
  EXPECT_TRUE(SIDRE_group_has_view(SIDRE_group_get_group_from_name(
                                     SIDRE_group_get_group_from_name(root,
                                                                     "fields"),
                                     "a"), "i0"));

// TODO - fix wrapping, change to datastore save call.
//  SIDRE_group_save(root, "C_out_sidre_group_save_restore_simple","conduit");

  SIDRE_datastore_print(ds);

  // Doesn't work yet.
#if 0
  SIDRE_datastore* ds2 = SIDRE_datastore_new();

  SIDRE_group_load(SIDRE_datastore_get_root(
                     ds2), "C_out_sidre_group_save_restore_simple",
                   "conduit");

  SIDRE_datastore_print(ds2);

  root = SIDRE_datastore_get_root(ds2);
  flds = SIDRE_group_get_group_from_name(root, "fields");

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  ga = SIDRE_group_get_group_from_name(flds, "a");
  i0_view = SIDRE_group_get_view_from_name(ga, "i0");
  EXPECT_EQ(SIDRE_view_get_data_int(i0_view), 1);

  SIDRE_datastore_print(ds2);

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);
#endif
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,rename_group)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* child1 = SIDRE_group_create_group(root, "g_a");
  SIDRE_group* child2 = SIDRE_group_create_group(root, "g_b");
  SIDRE_group* child3 = SIDRE_group_create_group(root, "g_c");

  bool success = SIDRE_group_rename(child1, "g_r");
  EXPECT_TRUE( success );
  EXPECT_TRUE( strcmp(SIDRE_group_get_name(child1), "g_r") == 0 );
  EXPECT_TRUE( SIDRE_group_has_group(root, "g_r") );
  EXPECT_FALSE( SIDRE_group_has_group(root, "g_a") );

  success = SIDRE_group_rename(child2, "fields/g_s");
  EXPECT_FALSE( success );
  EXPECT_TRUE( strcmp(SIDRE_group_get_name(child2), "g_b") == 0 );

  success = SIDRE_group_rename(child3, "g_b");
  EXPECT_FALSE( success );
  EXPECT_TRUE( strcmp(SIDRE_group_get_name(child3), "g_c") == 0 );

}

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a");
  SIDRE_group* gb = SIDRE_group_create_group(flds, "b");
  SIDRE_group* gc = SIDRE_group_create_group(flds, "c");

  SIDRE_group_create_view_scalar_int(ga, "i0", 1);
  SIDRE_group_create_view_scalar_float(gb, "f0", 100.0);
  SIDRE_group_create_view_scalar_double(gc, "d0", 3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

// TODO - Fix wrapping, change save call to use datastore class
//  SIDRE_group_save(root, "C_out_sidre_group_save_restore_complex",
//                       "conduit");

  SIDRE_datastore_print(ds);

  SIDRE_datastore* ds2 = SIDRE_datastore_new();
  root = SIDRE_datastore_get_root(ds2);

#if 0
  SIDRE_group_load(root, "C_out_sidre_group_save_restore_complex",
                   "conduit");

  flds = SIDRE_group_get_group_from_name(root, "fields");
  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

  ga = SIDRE_group_get_group_from_name(flds, "a");
  gb = SIDRE_group_get_group_from_name(flds, "b");
  gc = SIDRE_group_get_group_from_name(flds, "c");

  i0_view = SIDRE_group_get_view_from_name(ga, "i0");
  f0_view = SIDRE_group_get_view_from_name(gb, "f0");
  d0_view = SIDRE_group_get_view_from_name(gc, "d0");

  EXPECT_EQ(SIDRE_view_get_data_int(i0_view), 1);
  EXPECT_NEAR(SIDRE_view_get_data_float(f0_view), 100.0, 1e-12);
  EXPECT_NEAR(SIDRE_view_get_data_double(d0_view), 3000.0, 1e-12);

  SIDRE_datastore_print(ds2);

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);
#endif
}
