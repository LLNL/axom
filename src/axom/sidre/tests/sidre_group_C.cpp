// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"
#include "axom/slic/interface/c_fortran/wrapSLIC.h"

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_name)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, group_buf, group2_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* group = SIDRE_group_create_group(root, "test", &group_buf);

  //    EXPECT_TRUE(group->getName() == std::string("test") );
  EXPECT_TRUE(strcmp(SIDRE_group_get_name(group), "test") == 0);

  SIDRE_group* group2 =
    SIDRE_group_get_group_from_name(root, "foo", &group2_buf);
  EXPECT_TRUE(group2 == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_parent)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf, child_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_group* child = SIDRE_group_create_group(parent, "child", &child_buf);

  SIDRE_group tmp_buf;
  SIDRE_group* tmp = SIDRE_group_get_parent(child, &tmp_buf);
  EXPECT_TRUE( tmp->addr == parent->addr );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_datastore)
{
  SIDRE_datastore ds_buf, const_ds_buf;
  SIDRE_group root_buf, group_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* group = SIDRE_group_create_group(root, "parent", &group_buf);

  EXPECT_TRUE( SIDRE_group_get_data_store(group, &ds_buf) == ds );

  SIDRE_datastore const* const_ds = SIDRE_group_get_data_store(group,
                                                               &const_ds_buf);
  EXPECT_TRUE( const_ds->addr == ds->addr );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf, child_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_group* child = SIDRE_group_create_group(parent, "child", &child_buf);
  EXPECT_TRUE( SIDRE_group_get_parent(child, &parent_buf) == parent );

  SIDRE_group tmp_buf;
  SIDRE_group* tmp = SIDRE_group_get_group_from_name(parent, "child", &tmp_buf);
  EXPECT_TRUE(tmp->addr = child->addr);
  tmp = SIDRE_group_get_group_from_index(parent, 0, &tmp_buf);
  EXPECT_TRUE(tmp->addr = child->addr);

  // check error condition
  EXPECT_TRUE( SIDRE_group_get_group_from_name(parent,
                                               "non-existant group",
                                               &parent_buf) == NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf;
  SIDRE_view view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);

  SIDRE_view* view = SIDRE_group_create_view_empty(parent, "view", &view_buf);

  SIDRE_view tmp_buf;
  SIDRE_view* tmp = SIDRE_group_get_view_from_name(parent, "view", &tmp_buf);
  EXPECT_TRUE(tmp->addr == view->addr);
  tmp = SIDRE_group_get_view_from_index(parent, 0, &tmp_buf);
  EXPECT_TRUE(tmp->addr == view->addr);

  // check error condition
  EXPECT_TRUE( SIDRE_group_get_view_from_name(parent,
                                              "non-existant view", &view_buf) ==
               NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view_names_and_indicies)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf;
  SIDRE_view view1_buf, view2_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_view* view1 =
    SIDRE_group_create_view_empty(parent, "view1", &view1_buf);
  SIDRE_view* view2 =
    SIDRE_group_create_view_empty(parent, "view2", &view2_buf);

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf, group1_buf, group2_buf;
  SIDRE_group group1out_buf, group2out_buf, emptyGroup_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_group* group1 = SIDRE_group_create_group(parent, "group1", &group1_buf);
  SIDRE_group* group2 = SIDRE_group_create_group(parent, "group2", &group2_buf);
  EXPECT_EQ(SIDRE_group_get_num_groups(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_first_valid_group_index(parent);
  SIDRE_IndexType idx2 = SIDRE_group_get_next_valid_group_index(parent, idx1);
  SIDRE_IndexType idx3 = SIDRE_group_get_next_valid_group_index(parent, idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(SIDRE_InvalidIndex, idx3);

  SIDRE_group* group1out = SIDRE_group_get_group_from_index(parent, idx1,
                                                            &group1out_buf);
  SIDRE_group* group2out = SIDRE_group_get_group_from_index(parent, idx2,
                                                            &group2out_buf);
  EXPECT_EQ(group1->addr, group1out->addr);
  EXPECT_EQ(group2->addr, group2out->addr);

  // check error conditions
  SIDRE_group* emptyGroup =
    SIDRE_group_create_group(root, "emptyGroup", &emptyGroup_buf);

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf, emptyGroup_buf;
  SIDRE_view view1_buf, view2_buf;
  SIDRE_view view1out_buf, view2out_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_view* view1 =
    SIDRE_group_create_view_empty(parent, "view1", &view1_buf);
  SIDRE_view* view2 =
    SIDRE_group_create_view_empty(parent, "view2", &view2_buf);
  EXPECT_EQ(SIDRE_group_get_num_views(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_group_get_first_valid_view_index(parent);
  SIDRE_IndexType idx2 = SIDRE_group_get_next_valid_view_index(parent, idx1);
  SIDRE_IndexType idx3 = SIDRE_group_get_next_valid_view_index(parent, idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(SIDRE_InvalidIndex, idx3);

  SIDRE_view* view1out = SIDRE_group_get_view_from_index(parent, idx1,
                                                         &view1out_buf);
  SIDRE_view* view2out = SIDRE_group_get_view_from_index(parent, idx2,
                                                         &view2out_buf);
  EXPECT_EQ(view1->addr, view1out->addr);
  EXPECT_EQ(view2->addr, view2out->addr);

  // check error conditions
  SIDRE_group* emptyGroup =
    SIDRE_group_create_group(root, "emptyGroup", &emptyGroup_buf);

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, parent_buf, group1_buf, group2_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);

  SIDRE_group* parent = SIDRE_group_create_group(root, "parent", &parent_buf);
  SIDRE_group* group1 = SIDRE_group_create_group(parent, "group1", &group1_buf);
  SIDRE_group* group2 = SIDRE_group_create_group(parent, "group2", &group2_buf);

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, group_buf;
  SIDRE_view view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* group = SIDRE_group_create_group(root, "parent", &group_buf);

  SIDRE_view* view = SIDRE_group_create_view_empty(group, "view", &view_buf);
  EXPECT_TRUE( SIDRE_group_get_parent(group, &root_buf) == root );
  EXPECT_FALSE( SIDRE_view_has_buffer(view) );

  EXPECT_TRUE( SIDRE_group_has_view(group, "view") );
  // try creating view again, should be a no-op.
  EXPECT_TRUE( SIDRE_group_create_view_empty(group, "view",
                                             &view_buf) == NULL );

  SIDRE_group_destroy_view(group, "view");
  // destroy already destroyed group.  Should be a no-op, not a failure
  SIDRE_group_destroy_view(group, "view");

  EXPECT_FALSE( SIDRE_group_has_view(group, "view") );

  // try api call that specifies specific type and length
  SIDRE_group_create_view_and_allocate_nelems( group,
                                               "viewWithLength1",
                                               SIDRE_FLOAT_ID,
                                               50, &view_buf);

  // error condition check - try again with duplicate name, should be a no-op
  //XXX EXPECT_TRUE( SIDRE_group_create_view_and_allocate( group,
  // "viewWithLength1", SIDRE_FLOAT64_ID, 50) );
  SIDRE_group_destroy_view_and_data_name( group, "viewWithLength1");
  EXPECT_FALSE( SIDRE_group_has_view( group, "viewWithLength1") );

  EXPECT_TRUE( SIDRE_group_create_view_and_allocate_nelems( group,
                                                            "viewWithLengthBadLen",
                                                            SIDRE_FLOAT64_ID,
                                                            -1,
                                                            &view_buf) ==
               NULL );

  // try api call that specifies data type in another way
  SIDRE_group_create_view_and_allocate_nelems( group, "viewWithLength2",
                                               SIDRE_FLOAT64_ID, 50,
                                               &view_buf );
  EXPECT_TRUE( SIDRE_group_create_view_and_allocate_nelems( group,
                                                            "viewWithLength2",
                                                            SIDRE_FLOAT64_ID,
                                                            50,
                                                            &view_buf ) ==
               NULL );
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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, group_buf, group2_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* group = SIDRE_group_create_group(root, "group", &group_buf);
  EXPECT_TRUE( SIDRE_group_get_parent(group, &root_buf) == root );

  EXPECT_TRUE( SIDRE_group_has_group(root, "group") );

  SIDRE_group_destroy_group_name(root, "group");
  EXPECT_FALSE( SIDRE_group_has_group(root, "group") );

  SIDRE_group* group2 = SIDRE_group_create_group(root, "group2", &group2_buf);
  // shut up compiler about unused variable
  (void)group2;
  SIDRE_group_destroy_group_index( root, SIDRE_group_get_first_valid_group_index(
                                     root) );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, flds_buf, badGroup_buf;
  SIDRE_view view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields", &flds_buf);
  SIDRE_group_create_view_empty(flds, "a", &view_buf);

  EXPECT_TRUE(SIDRE_group_has_view(flds, "a"));

  // attempt to create duplicate group name

  SIDRE_group* badGroup =
    SIDRE_group_create_group(root, "fields", &badGroup_buf);
  EXPECT_TRUE( badGroup == NULL );

  // check error condition
  // attempt to create duplicate view name.
  EXPECT_TRUE(SIDRE_group_create_view_empty(flds, "a", &view_buf) == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
// Restore this after copy_move is working. ATK-667
#if 0
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, sub_buf;
  SIDRE_view view_buf

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields");

  SIDRE_group_create_view_scalar_int(flds, "i0", 1, &view_buf);
  SIDRE_group_create_view_scalar_float(flds, "f0", 100.0, &view_buf);
  SIDRE_group_create_view_scalar_double(flds, "d0", 3000.0, &view_buf);

  EXPECT_TRUE(SIDRE_group_has_view(flds, "i0"));
  EXPECT_TRUE(SIDRE_group_has_view(flds, "f0"));
  EXPECT_TRUE(SIDRE_group_has_view(flds, "d0"));

  // test moving a view form feds7 to sub
  // flds->createGroup("sub")->moveView(flds->getView("d0"));
  SIDRE_group* sub = SIDRE_group_create_group(flds, "sub", &sub_buf);
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
    SIDRE_view tmpview_buf;
    SIDRE_buffer tmpbuf_buf;
    SIDRE_view* tmpview =
      SIDRE_group_get_view_from_name(sub, "d0", &tmpview_buf);
    SIDRE_buffer* tmpbuf = SIDRE_view_get_buffer(tmpview, &tmpbuf_buf);
    d0_data = (double*) SIDRE_buffer_get_void_ptr(tmpbuf);
  }
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds top sub
  SIDRE_group_copy_view(sub,
                        SIDRE_group_get_view_from_name(flds, "i0", &view_buf));

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, flds_buf, sub_buf;
  SIDRE_group ga_buf, gb_buf, gc_buf;
  SIDRE_view i0_view_buf, f0_view_buf, d0_view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields", &flds_buf);

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a", &ga_buf);
  SIDRE_group* gb = SIDRE_group_create_group(flds, "b", &gb_buf);
  SIDRE_group* gc = SIDRE_group_create_group(flds, "c", &gc_buf);

  SIDRE_view* i0_view =
    SIDRE_group_create_view_and_allocate_nelems(ga, "i0",
                                                SIDRE_INT_ID, 1, &i0_view_buf);
  SIDRE_view* f0_view =
    SIDRE_group_create_view_and_allocate_nelems(gb, "f0",
                                                SIDRE_FLOAT_ID, 1,
                                                &f0_view_buf);
  SIDRE_view* d0_view =
    SIDRE_group_create_view_and_allocate_nelems(gc, "d0",
                                                SIDRE_DOUBLE_ID, 1,
                                                &d0_view_buf);

  SIDRE_view_set_scalar_int(i0_view, 1);
  SIDRE_view_set_scalar_float(f0_view, 100.0);
  SIDRE_view_set_scalar_double(d0_view, 3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

  // move "b" to a child of "sub"
  SIDRE_group* sub = SIDRE_group_create_group(flds, "sub", &sub_buf);
  SIDRE_group_move_group(sub, gb, &gb_buf);

  SIDRE_group_print(flds);

  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "sub"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

  SIDRE_group tmpgrp_buf, tmpgrp2_buf;
  SIDRE_group* tmpgrp =
    SIDRE_group_get_group_from_name(flds, "sub", &tmpgrp_buf);
  SIDRE_group* tmpgrp2 =
    SIDRE_group_get_group_from_name(tmpgrp, "b", &tmpgrp2_buf);
  EXPECT_EQ(tmpgrp2->addr, gb->addr);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, grp_buf;
  SIDRE_view view_buf, view1_buf, view2_buf;
  SIDRE_buffer buffer1_buf, tmpbuf_buf;

  SIDRE_datastore* const ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* const grp = SIDRE_group_create_group(root, "grp", &grp_buf);

  const char* viewName1 = "viewBuffer1";
  const char* viewName2 = "viewBuffer2";

  // XXX const
  SIDRE_view* view1 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName1,
                                                SIDRE_INT_ID, 1, &view1_buf);
  SIDRE_view* view2 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName2,
                                                SIDRE_FLOAT_ID, 1, &view2_buf);

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName1));
  SIDRE_view* view = SIDRE_group_get_view_from_name(grp, viewName1, &view_buf);
  EXPECT_EQ(view->addr, view1->addr);

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName2));
  view = SIDRE_group_get_view_from_name(grp, viewName2, &view_buf);
  EXPECT_EQ(view->addr, view2->addr);

  SIDRE_buffer* tmpbuf = SIDRE_view_get_buffer(view1, &tmpbuf_buf);
  SIDRE_IndexType bufferId1 = SIDRE_buffer_get_index(tmpbuf);

  SIDRE_group_destroy_view_and_data_name(grp, viewName1);


  EXPECT_FALSE(SIDRE_group_has_view(grp, viewName1));
  EXPECT_EQ(SIDRE_datastore_get_num_buffers(ds), 1u);

  SIDRE_buffer* buffer1 =
    SIDRE_datastore_get_buffer(ds, bufferId1, &buffer1_buf);
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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, grp_buf;
  SIDRE_view view1_buf;

  SIDRE_datastore* const ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* const grp = SIDRE_group_create_group(root, "grp", &grp_buf);

  const char* viewName1 = "viewBuffer1";

  // use create + alloc convenience methods
  // this one is the DataType & method
  SIDRE_view* const view1 =
    SIDRE_group_create_view_and_allocate_nelems(grp, viewName1,
                                                SIDRE_INT_ID, 10, &view1_buf);

  EXPECT_TRUE(SIDRE_group_has_view(grp, viewName1));
  SIDRE_view viewtmp_buf;
  SIDRE_view* viewtmp = SIDRE_group_get_view_from_name(grp, viewName1,
                                                       &viewtmp_buf);
  EXPECT_EQ(viewtmp->addr, view1->addr);


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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf, grp_buf;
  SIDRE_group flds_buf, ga_buf;
  SIDRE_view view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields", &flds_buf);

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a", &ga_buf);

  SIDRE_group_create_view_scalar_int(ga, "i0", 1, &view_buf);

  EXPECT_TRUE(SIDRE_group_has_group(root, "fields"));
  EXPECT_TRUE(SIDRE_group_has_group(SIDRE_group_get_group_from_name(root,
                                                                    "fields",
                                                                    &grp_buf),
                                    "a"));
  EXPECT_TRUE(SIDRE_group_has_view(SIDRE_group_get_group_from_name(
                                     SIDRE_group_get_group_from_name(root,
                                                                     "fields",
                                                                     &grp_buf),
                                     "a", &grp_buf), "i0"));

// TODO - fix wrapping, change to datastore save call.
//  SIDRE_group_save(root, "C_out_sidre_group_save_restore_simple","conduit");

  SIDRE_datastore_print(ds);

  // Doesn't work yet.
#if 0
  SIDRE_datastore* ds2 = SIDRE_datastore_new(ds2_buf);

  SIDRE_group_load(SIDRE_datastore_get_root(
                     ds2, &root_buf), "C_out_sidre_group_save_restore_simple",
                   "conduit");

  SIDRE_datastore_print(ds2);

  root = SIDRE_datastore_get_root(ds2, &root_buf);
  flds = SIDRE_group_get_group_from_name(root, "fields", &flds_buf);

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
  SIDRE_datastore ds_buf;
  SIDRE_group root_buf;
  SIDRE_group child1_buf, child2_buf, child3_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* child1 = SIDRE_group_create_group(root, "g_a", &child1_buf);
  SIDRE_group* child2 = SIDRE_group_create_group(root, "g_b", &child2_buf);
  SIDRE_group* child3 = SIDRE_group_create_group(root, "g_c", &child3_buf);

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
  SIDRE_datastore ds_buf, ds2_buf;
  SIDRE_group root_buf, flds_buf;
  SIDRE_group ga_buf, gb_buf, gc_buf;
  SIDRE_view view_buf;

  SIDRE_datastore* ds = SIDRE_datastore_new(&ds_buf);
  SIDRE_group* root = SIDRE_datastore_get_root(ds, &root_buf);
  SIDRE_group* flds = SIDRE_group_create_group(root, "fields", &flds_buf);

  SIDRE_group* ga = SIDRE_group_create_group(flds, "a", &ga_buf);
  SIDRE_group* gb = SIDRE_group_create_group(flds, "b", &gb_buf);
  SIDRE_group* gc = SIDRE_group_create_group(flds, "c", &gc_buf);

  SIDRE_group_create_view_scalar_int(ga, "i0", 1, &view_buf);
  SIDRE_group_create_view_scalar_float(gb, "f0", 100.0, &view_buf);
  SIDRE_group_create_view_scalar_double(gc, "d0", 3000.0, &view_buf);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_group_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_group_has_group(flds, "c"));

// TODO - Fix wrapping, change save call to use datastore class
//  SIDRE_group_save(root, "C_out_sidre_group_save_restore_complex",
//                       "conduit");

  SIDRE_datastore_print(ds);

  SIDRE_datastore* ds2 = SIDRE_datastore_new(&ds2_buf);
  root = SIDRE_datastore_get_root(ds2, &root_buf);

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
