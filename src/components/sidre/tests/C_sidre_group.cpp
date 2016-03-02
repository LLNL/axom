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
#include "slic/wrapSLIC.h"

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_name)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * group = SIDRE_datagroup_create_group(root, "test");

  //    EXPECT_TRUE(group->getName() == std::string("test") );
  EXPECT_TRUE(strcmp(SIDRE_datagroup_get_name(group), "test") == 0);

  SIDRE_datagroup * group2 = SIDRE_datagroup_get_group(root, "foo");
  EXPECT_TRUE(group2 == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_parent)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");
  SIDRE_datagroup * child = SIDRE_datagroup_create_group(parent, "child");

  EXPECT_TRUE( SIDRE_datagroup_get_parent(child) == parent );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_datastore)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * group = SIDRE_datagroup_create_group(root, "parent");

  EXPECT_TRUE( SIDRE_datagroup_get_data_store(group) == ds );

  SIDRE_datastore const * const_ds = SIDRE_datagroup_get_data_store(group);
  EXPECT_TRUE( const_ds == ds );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");
  SIDRE_datagroup * child = SIDRE_datagroup_create_group(parent, "child");
  EXPECT_TRUE( SIDRE_datagroup_get_parent(child) == parent );

  EXPECT_TRUE( SIDRE_datagroup_has_group(parent, "child") );
  // check error condition
  EXPECT_TRUE( SIDRE_datagroup_get_group(parent,
                                         "non-existant group") == NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");

  SIDRE_dataview * view = SIDRE_datagroup_create_view_empty(parent, "view");

  EXPECT_TRUE( SIDRE_datagroup_get_view_from_name(parent, "view") == view );

  // check error condition
  EXPECT_TRUE( SIDRE_datagroup_get_view_from_name(parent,
                                                  "non-existant view") ==
               NULL );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view_names_and_indicies)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");
  SIDRE_dataview * view1 = SIDRE_datagroup_create_view_empty(parent, "view1");
  SIDRE_dataview * view2 = SIDRE_datagroup_create_view_empty(parent, "view2");

  EXPECT_EQ(SIDRE_datagroup_get_num_views(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_datagroup_get_view_index(parent, "view1");
  SIDRE_IndexType idx2 = SIDRE_datagroup_get_view_index(parent, "view2");

  const char * name1 = SIDRE_datagroup_get_view_name(parent, idx1);
  const char * name2 = SIDRE_datagroup_get_view_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "view1") == 0);
  EXPECT_TRUE(strcmp(SIDRE_dataview_get_name(view1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "view2") == 0);
  EXPECT_TRUE(strcmp(SIDRE_dataview_get_name(view2), name2) == 0);

  // check error conditions
  SIDRE_IndexType idx3 = SIDRE_datagroup_get_view_index(parent, "view3");
  EXPECT_TRUE(idx3 == SIDRE_InvalidIndex);

  const char * name3 = SIDRE_datagroup_get_view_name(parent, idx3);
  EXPECT_TRUE(name3 == NULL);
  EXPECT_FALSE(SIDRE_name_is_valid(name3));

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_view_index)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");
  SIDRE_dataview * view1 = SIDRE_datagroup_create_view_empty(parent, "view1");
  SIDRE_dataview * view2 = SIDRE_datagroup_create_view_empty(parent, "view2");

  SIDRE_datagroup * emptyGroup =
    SIDRE_datagroup_create_group(root, "emptyGroup");

  EXPECT_EQ(SIDRE_datagroup_get_num_views(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_datagroup_get_first_valid_view_index(parent);
  SIDRE_IndexType idx2 =
    SIDRE_datagroup_get_next_valid_view_index(parent, idx1);

  const char * name1 = SIDRE_datagroup_get_view_name(parent, idx1);
  const char * name2 = SIDRE_datagroup_get_view_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "view1") == 0);
  EXPECT_TRUE(strcmp(SIDRE_dataview_get_name(view1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "view2") == 0);
  EXPECT_TRUE(strcmp(SIDRE_dataview_get_name(view2), name2) == 0);

  // check error conditions
  SIDRE_IndexType badidx1 = SIDRE_datagroup_get_first_valid_view_index(
    emptyGroup);
  SIDRE_IndexType badidx2 = SIDRE_datagroup_get_next_valid_view_index(
    emptyGroup, badidx1);

  EXPECT_TRUE(badidx1 == SIDRE_InvalidIndex);
  EXPECT_TRUE(badidx2 == SIDRE_InvalidIndex);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group_name_index)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_datagroup * parent = SIDRE_datagroup_create_group(root, "parent");
  SIDRE_datagroup * group1 = SIDRE_datagroup_create_group(parent, "group1");
  SIDRE_datagroup * group2 = SIDRE_datagroup_create_group(parent, "group2");

  EXPECT_EQ(SIDRE_datagroup_get_num_groups(parent), 2u);

  SIDRE_IndexType idx1 = SIDRE_datagroup_get_group_index(parent, "group1");
  SIDRE_IndexType idx2 = SIDRE_datagroup_get_group_index(parent, "group2");

  const char * name1 = SIDRE_datagroup_get_group_name(parent, idx1);
  const char * name2 = SIDRE_datagroup_get_group_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "group1") == 0);
  EXPECT_TRUE(strcmp(SIDRE_datagroup_get_name(group1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "group2") == 0);
  EXPECT_TRUE(strcmp(SIDRE_datagroup_get_name(group2), name2) == 0);

  // check error conditions
  SIDRE_IndexType idx3 = SIDRE_datagroup_get_group_index(parent, "group3");
  EXPECT_TRUE(idx3 == SIDRE_InvalidIndex);

  const char * name3 = SIDRE_datagroup_get_group_name(parent, idx3);
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
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * group = SIDRE_datagroup_create_group(root, "parent");

  SIDRE_dataview * view = SIDRE_datagroup_create_view_empty(group, "view");
  EXPECT_TRUE( SIDRE_datagroup_get_parent(group) == root );
  EXPECT_FALSE( SIDRE_dataview_has_buffer(view) );

  EXPECT_TRUE( SIDRE_datagroup_has_view(group, "view") );
  // try creating view again, should be a no-op.
  EXPECT_TRUE( SIDRE_datagroup_create_view_empty(group, "view") == NULL );

  SIDRE_datagroup_destroy_view(group, "view");
  // destroy already destroyed group.  Should be a no-op, not a failure
  SIDRE_datagroup_destroy_view(group, "view");

  EXPECT_FALSE( SIDRE_datagroup_has_view(group, "view") );

  // try api call that specifies specific type and length
  SIDRE_datagroup_create_view_and_allocate_nelems( group,
                                                   "viewWithLength1",
                                                   SIDRE_FLOAT_ID,
                                                   50);

  // error condition check - try again with duplicate name, should be a no-op
  //XXX EXPECT_TRUE( SIDRE_datagroup_create_view_and_allocate( group, "viewWithLength1", SIDRE_FLOAT64_ID, 50) );
  SIDRE_datagroup_destroy_view_and_data_name( group, "viewWithLength1");
  EXPECT_FALSE( SIDRE_datagroup_has_view( group, "viewWithLength1") );

  EXPECT_TRUE( SIDRE_datagroup_create_view_and_allocate_nelems( group,
                                                                "viewWithLengthBadLen",
                                                                SIDRE_FLOAT64_ID,
                                                                -1) == NULL );

  // try api call that specifies data type in another way
  SIDRE_datagroup_create_view_and_allocate_nelems( group, "viewWithLength2",
                                                   SIDRE_FLOAT64_ID, 50 );
  EXPECT_TRUE( SIDRE_datagroup_create_view_and_allocate_nelems( group,
                                                                "viewWithLength2",
                                                                SIDRE_FLOAT64_ID,
                                                                50 ) == NULL );
  // destroy this view using index
  SIDRE_datagroup_destroy_view_and_data_index( group, SIDRE_datagroup_get_first_valid_view_index(
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
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * group = SIDRE_datagroup_create_group(root, "group");
  EXPECT_TRUE( SIDRE_datagroup_get_parent(group) == root );

  EXPECT_TRUE( SIDRE_datagroup_has_group(root, "group") );

  SIDRE_datagroup_destroy_group_name(root, "group");
  EXPECT_FALSE( SIDRE_datagroup_has_group(root, "group") );

  SIDRE_datagroup * group2 = SIDRE_datagroup_create_group(root, "group2");
  // shut up compiler about unused variable
  (void)group2;
  SIDRE_datagroup_destroy_group_index( root, SIDRE_datagroup_get_first_valid_group_index(
                                         root) );

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * flds = SIDRE_datagroup_create_group(root, "fields");
  SIDRE_datagroup_create_view_empty(flds, "a");

  EXPECT_TRUE(SIDRE_datagroup_has_view(flds, "a"));

  // attempt to create duplicate group name

  SIDRE_datagroup * badGroup = SIDRE_datagroup_create_group(root, "fields");
  EXPECT_TRUE( badGroup == NULL );

  // check error condition
  // attempt to create duplicate view name.
  EXPECT_TRUE(SIDRE_datagroup_create_view_empty(flds, "a") == NULL);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * flds = SIDRE_datagroup_create_group(root, "fields");

  SIDRE_dataview * i0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(flds, "i0",
                                                    SIDRE_INT_ID, 1);
  SIDRE_dataview * f0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(flds, "f0",
                                                    SIDRE_FLOAT_ID, 1);
  SIDRE_dataview * d0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(flds, "d0",
                                                    SIDRE_DOUBLE_ID, 1);

  SIDRE_dataview_set_scalar_int(i0_view, 1);
  SIDRE_dataview_set_scalar_float(f0_view, 100.0);
  SIDRE_dataview_set_scalar_double(d0_view, 3000.0);

  EXPECT_TRUE(SIDRE_datagroup_has_view(flds, "i0"));
  EXPECT_TRUE(SIDRE_datagroup_has_view(flds, "f0"));
  EXPECT_TRUE(SIDRE_datagroup_has_view(flds, "d0"));

  // test moving a view form feds7 to sub
  // flds->createGroup("sub")->moveView(flds->getView("d0"));
  SIDRE_datagroup * sub = SIDRE_datagroup_create_group(flds, "sub");
  SIDRE_datagroup_move_view(sub,
                            SIDRE_datagroup_get_view_from_name(flds, "d0"));
  SIDRE_datagroup_print(flds);
  EXPECT_FALSE(SIDRE_datagroup_has_view(flds, "d0"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "sub"));
  EXPECT_TRUE(SIDRE_datagroup_has_view(sub, "d0"));

  // check the data value
  // double *d0_data =  flds->getGroup("sub")
  //                        ->getView("d0")
  //                        ->getNode().as_double_ptr();
  double * d0_data;
  {
    SIDRE_dataview * tmpview = SIDRE_datagroup_get_view_from_name(sub, "d0");
    SIDRE_databuffer * tmpbuf = SIDRE_dataview_get_buffer(tmpview);
    d0_data = (double *) SIDRE_databuffer_get_void_ptr(tmpbuf);
  }
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds top sub
  SIDRE_datagroup_copy_view(sub,
                            SIDRE_datagroup_get_view_from_name(flds, "i0"));

  SIDRE_datagroup_print(flds);

  EXPECT_TRUE(SIDRE_datagroup_has_view(flds, "i0"));
  EXPECT_TRUE(SIDRE_datagroup_has_view(sub, "i0"));

#ifdef XXX
  // we expect the actual data  pointers to be the same
  EXPECT_EQ(SIDRE_datagroup_get_view(flds, "i0")->getNode().data_pointer(),
            SIDRE_datagroup_get_group("sub")->get_view(
              "i0")->getNode().data_pointer());
#endif

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * flds = SIDRE_datagroup_create_group(root, "fields");

  SIDRE_datagroup * ga = SIDRE_datagroup_create_group(flds, "a");
  SIDRE_datagroup * gb = SIDRE_datagroup_create_group(flds, "b");
  SIDRE_datagroup * gc = SIDRE_datagroup_create_group(flds, "c");

  SIDRE_dataview * i0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(ga, "i0",
                                                    SIDRE_INT_ID, 1);
  SIDRE_dataview * f0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(gb, "f0",
                                                    SIDRE_FLOAT_ID, 1);
  SIDRE_dataview * d0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(gc, "d0",
                                                    SIDRE_DOUBLE_ID, 1);

  SIDRE_dataview_set_scalar_int(i0_view, 1);
  SIDRE_dataview_set_scalar_float(f0_view, 100.0);
  SIDRE_dataview_set_scalar_double(d0_view, 3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "c"));

  //move "b" to a child of "sub"
  SIDRE_datagroup * sub = SIDRE_datagroup_create_group(flds, "sub");
  SIDRE_datagroup_move_group(sub, gb);

  SIDRE_datagroup_print(flds);

  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "sub"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "c"));

  SIDRE_datagroup * tmpgrp = SIDRE_datagroup_get_group(flds, "sub");
  EXPECT_EQ(SIDRE_datagroup_get_group(tmpgrp, "b"), gb);

  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  SIDRE_datastore * const ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * const grp = SIDRE_datagroup_create_group(root, "grp");

  const char * viewName1 = "viewBuffer1";
  const char * viewName2 = "viewBuffer2";

  // XXX const
  SIDRE_dataview * view1 =
    SIDRE_datagroup_create_view_and_allocate_nelems(grp, viewName1,
                                                    SIDRE_INT_ID, 1);
  SIDRE_dataview * view2 =
    SIDRE_datagroup_create_view_and_allocate_nelems(grp, viewName2,
                                                    SIDRE_FLOAT_ID, 1);

  EXPECT_TRUE(SIDRE_datagroup_has_view(grp, viewName1));
  EXPECT_EQ( SIDRE_datagroup_get_view_from_name(grp, viewName1), view1 );

  EXPECT_TRUE(SIDRE_datagroup_has_view(grp, viewName2));
  EXPECT_EQ( SIDRE_datagroup_get_view_from_name(grp, viewName2), view2 );

  SIDRE_databuffer * tmpbuf = SIDRE_dataview_get_buffer(view1);
  SIDRE_IndexType bufferId1 = SIDRE_databuffer_get_index(tmpbuf);

  SIDRE_datagroup_destroy_view_and_data_name(grp, viewName1);


  EXPECT_FALSE(SIDRE_datagroup_has_view(grp, viewName1));
  EXPECT_EQ(SIDRE_datastore_get_num_buffers(ds), 1u);

  SIDRE_databuffer * buffer1 = SIDRE_datastore_get_buffer(ds, bufferId1);
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
  SIDRE_datastore * const ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * const grp = SIDRE_datagroup_create_group(root, "grp");

  const char * viewName1 = "viewBuffer1";

  // use create + alloc convenience methods
  // this one is the DataType & method
  SIDRE_dataview * const view1 =
    SIDRE_datagroup_create_view_and_allocate_nelems(grp, viewName1,
                                                    SIDRE_INT_ID, 10);

  EXPECT_TRUE(SIDRE_datagroup_has_view(grp, viewName1));
  EXPECT_EQ( SIDRE_datagroup_get_view_from_name(grp, viewName1), view1 );


#ifdef XXX
  int * v1_vals = (int *) SIDRE_dataview_get_void_ptr(view1);
  double * v2_vals = (double *)  SIDRE_dataview_get_void_ptr(view1);

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
    v2_vals[i] = i * 3.1415;
  }
#endif

  EXPECT_EQ(SIDRE_dataview_get_num_elements(view1), 10u);
  EXPECT_EQ(SIDRE_dataview_get_total_bytes(view1), 10 * sizeof(int));

  SIDRE_datagroup_destroy_view_and_data_name(grp, viewName1);

  SIDRE_datastore_delete(ds);
}

#ifdef XXX
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_view_of_buffer_with_datatype)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  // use create + alloc convenience methods
  // this one is the DataType & method
  SIDRE_dataview * base =
    SIDRE_datagroup_create_view_and_allocate_nelems(root, "base",
                                                    SIDRE_INT_ID, 10);
#ifdef XXX
  int * base_vals = (int *) SIDRE_dataview_get_void_ptr(base);
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

  SIDRE_databuffer * base_buff = SIDRE_dataview_get_buffer(base);
  // create two views into this buffer
  // view for the first 5 values
  SIDRE_datagroup_create_view(root, "sub_a", base_buff, SIDRE_C_INT_T, 5);
  // view for the second 5 values

  int * sub_a_vals = (int *) SIDRE_dataview_get_void_ptr(SIDRE_datagroup_get_view_from_name(
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
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * flds = SIDRE_datagroup_create_group(root, "fields");

  SIDRE_datagroup * ga = SIDRE_datagroup_create_group(flds, "a");

  SIDRE_dataview * i0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(ga, "i0",
                                                    SIDRE_INT_ID, 1);
  SIDRE_dataview_set_scalar_int(i0_view, 1);

  EXPECT_TRUE(SIDRE_datagroup_has_group(root, "fields"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(SIDRE_datagroup_get_group(root,
                                                                  "fields"),
                                        "a"));
  EXPECT_TRUE(SIDRE_datagroup_has_view(SIDRE_datagroup_get_group(
                                         SIDRE_datagroup_get_group(root,
                                                                   "fields"),
                                         "a"), "i0"));


  SIDRE_datagroup_save(root, "C_out_sidre_group_save_restore_simple","conduit");

  SIDRE_datastore_print(ds);

  SIDRE_datastore * ds2 = SIDRE_datastore_new();

  SIDRE_datagroup_load(SIDRE_datastore_get_root(
                         ds2), "C_out_sidre_group_save_restore_simple",
                       "conduit");

  SIDRE_datastore_print(ds2);

  root = SIDRE_datastore_get_root(ds2);
  flds = SIDRE_datagroup_get_group(root, "fields");

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "a"));
  ga = SIDRE_datagroup_get_group(flds, "a");
  i0_view = SIDRE_datagroup_get_view_from_name(ga, "i0");
  EXPECT_EQ(SIDRE_dataview_get_data_int(i0_view), 1);

  SIDRE_datastore_print(ds2);

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);

}

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_datagroup * flds = SIDRE_datagroup_create_group(root, "fields");

  SIDRE_datagroup * ga = SIDRE_datagroup_create_group(flds, "a");
  SIDRE_datagroup * gb = SIDRE_datagroup_create_group(flds, "b");
  SIDRE_datagroup * gc = SIDRE_datagroup_create_group(flds, "c");

  SIDRE_dataview * i0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(ga, "i0",
                                                    SIDRE_INT_ID, 1);
  SIDRE_dataview_set_scalar_int(i0_view, 1);

  SIDRE_dataview * f0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(gb,
                                                    "f0",
                                                    SIDRE_FLOAT_ID,
                                                    1);
  SIDRE_dataview_set_scalar_float(f0_view, 100.0);

  SIDRE_dataview * d0_view =
    SIDRE_datagroup_create_view_and_allocate_nelems(gc,
                                                    "d0",
                                                    SIDRE_DOUBLE_ID,
                                                    1);
  SIDRE_dataview_set_scalar_double(d0_view, 3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "c"));

  SIDRE_datagroup_save(root, "C_out_sidre_group_save_restore_complex",
                       "conduit");

  SIDRE_datastore_print(ds);

  SIDRE_datastore * ds2 = SIDRE_datastore_new();
  root = SIDRE_datastore_get_root(ds2);

  SIDRE_datagroup_load(root, "C_out_sidre_group_save_restore_complex",
                       "conduit");

  flds = SIDRE_datagroup_get_group(root, "fields");
  // check that all sub groups exist
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(SIDRE_datagroup_has_group(flds, "c"));

  ga = SIDRE_datagroup_get_group(flds, "a");
  gb = SIDRE_datagroup_get_group(flds, "b");
  gc = SIDRE_datagroup_get_group(flds, "c");

  i0_view = SIDRE_datagroup_get_view_from_name(ga, "i0");
  f0_view = SIDRE_datagroup_get_view_from_name(gb, "f0");
  d0_view = SIDRE_datagroup_get_view_from_name(gc, "d0");

  EXPECT_EQ(SIDRE_dataview_get_data_int(i0_view), 1);
  EXPECT_NEAR(SIDRE_dataview_get_data_float(f0_view), 100.0, 1e-12);
  EXPECT_NEAR(SIDRE_dataview_get_data_double(d0_view), 3000.0, 1e-12);

  SIDRE_datastore_print(ds2);

  SIDRE_datastore_delete(ds);
  SIDRE_datastore_delete(ds2);

}
