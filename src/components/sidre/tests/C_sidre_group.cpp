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


// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_name)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * group = ATK_datagroup_create_group(root, "test");

  //    EXPECT_TRUE(group->getName() == std::string("test") );
  EXPECT_TRUE(strcmp(ATK_datagroup_get_name(group), "test") == 0);

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_parent)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * parent = ATK_datagroup_create_group(root, "parent");
  ATK_datagroup * child = ATK_datagroup_create_group(parent, "child");

  EXPECT_TRUE( ATK_datagroup_get_parent(child) == parent );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_datastore)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * group = ATK_datagroup_create_group(root, "parent");

  EXPECT_TRUE( ATK_datagroup_get_data_store(group) == ds );

  ATK_datastore const * const_ds = ATK_datagroup_get_data_store(group);
  EXPECT_TRUE( const_ds == ds );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify hasGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,has_group)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * parent = ATK_datagroup_create_group(root, "parent");
  ATK_datagroup * child = ATK_datagroup_create_group(parent, "child");
  EXPECT_TRUE( ATK_datagroup_get_parent(child) == parent );

  EXPECT_TRUE( ATK_datagroup_has_group(parent, "child") );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,has_view)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * parent = ATK_datagroup_create_group(root, "parent");

  ATK_dataview * view = ATK_datagroup_create_view_and_buffer_simple(parent, "view");

  EXPECT_TRUE( ATK_dataview_get_owning_group(view) == parent );

  EXPECT_TRUE( ATK_datagroup_has_view(parent, "view") );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view_name_index)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * parent = ATK_datagroup_create_group(root, "parent");
  ATK_dataview * view1 = ATK_datagroup_create_view_and_buffer_simple(parent, "view1");
  ATK_dataview * view2 = ATK_datagroup_create_view_and_buffer_simple(parent, "view2");

  EXPECT_EQ(ATK_datagroup_get_num_views(parent), 2u);

  ATK_IndexType idx1 = ATK_datagroup_get_view_index(parent, "view1");
  ATK_IndexType idx2 = ATK_datagroup_get_view_index(parent, "view2");

  const char * name1 = ATK_datagroup_get_view_name(parent, idx1);
  const char * name2 = ATK_datagroup_get_view_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "view1") == 0);
  EXPECT_TRUE(strcmp(ATK_dataview_get_name(view1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "view2") == 0);
  EXPECT_TRUE(strcmp(ATK_dataview_get_name(view2), name2) == 0);

  ATK_IndexType idx3 = ATK_datagroup_get_view_index(parent, "view3");
  EXPECT_TRUE(idx3 == ATK_InvalidID);

#if 0 // C API needs some additions to make this work...
  const char * name3 = ATK_datagroup_get_view_name(parent, idx3);
  EXPECT_TRUE(name3 == NULL);
  EXPECT_TRUE(ATK_isNameValid(name3) == 0);
#endif

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group_name_index)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * parent = ATK_datagroup_create_group(root, "parent");
  ATK_datagroup * group1 = ATK_datagroup_create_group(parent, "group1");
  ATK_datagroup * group2 = ATK_datagroup_create_group(parent, "group2");

  EXPECT_EQ(ATK_datagroup_get_num_groups(parent), 2u);

  ATK_IndexType idx1 = ATK_datagroup_get_group_index(parent, "group1");
  ATK_IndexType idx2 = ATK_datagroup_get_group_index(parent, "group2");

  const char * name1 = ATK_datagroup_get_group_name(parent, idx1);
  const char * name2 = ATK_datagroup_get_group_name(parent, idx2);

  EXPECT_TRUE(strcmp(name1, "group1") == 0);
  EXPECT_TRUE(strcmp(ATK_datagroup_get_name(group1), name1) == 0);

  EXPECT_TRUE(strcmp(name2, "group2") == 0);
  EXPECT_TRUE(strcmp(ATK_datagroup_get_name(group2), name2) == 0);

  ATK_IndexType idx3 = ATK_datagroup_get_group_index(parent, "group3");
  EXPECT_TRUE(idx3 == ATK_InvalidID);

#if 0 // C API needs some additions to make this work...
  const char * name3 = ATK_datagroup_get_group_name(parent, idx3);
  EXPECT_TRUE(name3 == NULL);
  EXPECT_TRUE(ATK_isNameValid(name3) == 0);
#endif

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// createViewAndBuffer()
// destroyViewAndBuffer()
// hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_viewbuffer)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * group = ATK_datagroup_create_group(root, "parent");

  ATK_dataview * view = ATK_datagroup_create_view_and_buffer_simple(group, "view");
  EXPECT_TRUE( ATK_datagroup_get_parent(group) == root );
  EXPECT_TRUE( ATK_dataview_has_buffer(view) );

  EXPECT_TRUE( ATK_datagroup_has_view(group, "view") );

  ATK_datagroup_destroy_view_and_buffer(group, "view");

  EXPECT_FALSE( ATK_datagroup_has_view(group, "view") );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_group)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * group = ATK_datagroup_create_group(root, "group");
  EXPECT_TRUE( ATK_datagroup_get_parent(group) == root );

  EXPECT_TRUE( ATK_datagroup_has_group(root, "group") );


  ATK_datagroup_destroy_group(root, "group");
  EXPECT_FALSE( ATK_datagroup_has_group(root, "group") );

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * flds = ATK_datagroup_create_group(root, "fields");
  ATK_datagroup_create_view_and_buffer_simple(flds, "a");

  EXPECT_TRUE(ATK_datagroup_has_view(flds, "a"));

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * flds = ATK_datagroup_create_group(root, "fields");

  ATK_dataview * tmpview;
  tmpview = ATK_datagroup_create_view_and_buffer_simple(flds, "i0");
  ATK_dataview_allocate(tmpview, ATK_C_INT_T, 1);
  tmpview = ATK_datagroup_create_view_and_buffer_simple(flds, "f0");
  ATK_dataview_allocate(tmpview, ATK_C_FLOAT_T, 1);
  tmpview = ATK_datagroup_create_view_and_buffer_simple(flds, "d0");
  ATK_dataview_allocate(tmpview, ATK_C_DOUBLE_T, 1);

#ifdef XXX
  (*ATK_datagroup_get_view(flds, "i0")->getNode().as_int_ptr())   = 1;
  (*ATK_datagroup_get_view(flds, "f0")->getNode().as_float_ptr()) = 100.0;
  (*ATK_datagroup_get_view(flds, "d0")->getNode().as_double_ptr()) = 3000.0;
#endif
  {
    ATK_dataview * tmpview = ATK_datagroup_get_view(flds, "i0");
    ATK_databuffer * tmpbuf = ATK_dataview_get_buffer(tmpview);
    int * v = (int *) ATK_databuffer_get_data(tmpbuf);
    *v = 1;
  }
  {
    ATK_dataview * tmpview = ATK_datagroup_get_view(flds, "f0");
    ATK_databuffer * tmpbuf = ATK_dataview_get_buffer(tmpview);
    float * v = (float *) ATK_databuffer_get_data(tmpbuf);
    *v = 100.0;
  }
  {
    ATK_dataview * tmpview = ATK_datagroup_get_view(flds, "d0");
    ATK_databuffer * tmpbuf = ATK_dataview_get_buffer(tmpview);
    double * v = (double *) ATK_databuffer_get_data(tmpbuf);
    *v = 3000.0;
  }

  EXPECT_TRUE(ATK_datagroup_has_view(flds, "i0"));
  EXPECT_TRUE(ATK_datagroup_has_view(flds, "f0"));
  EXPECT_TRUE(ATK_datagroup_has_view(flds, "d0"));

  // test moving a view form feds7 to sub
  // flds->createGroup("sub")->moveView(flds->getView("d0"));
  ATK_datagroup * sub = ATK_datagroup_create_group(flds, "sub");
  ATK_datagroup_move_view(sub, ATK_datagroup_get_view(flds, "d0"));
  ATK_datagroup_print(flds);
  EXPECT_FALSE(ATK_datagroup_has_view(flds, "d0"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "sub"));
  EXPECT_TRUE(ATK_datagroup_has_view(sub, "d0"));

  // check the data value
  // double *d0_data =  flds->getGroup("sub")
  //                        ->getView("d0")
  //                        ->getNode().as_double_ptr();
  double * d0_data;
  {
    ATK_dataview * tmpview = ATK_datagroup_get_view(sub, "d0");
    ATK_databuffer * tmpbuf = ATK_dataview_get_buffer(tmpview);
    d0_data = (double *) ATK_databuffer_get_data(tmpbuf);
  }
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds top sub
  ATK_datagroup_copy_view(sub, ATK_datagroup_get_view(flds, "i0"));

  ATK_datagroup_print(flds);

  EXPECT_TRUE(ATK_datagroup_has_view(flds, "i0"));
  EXPECT_TRUE(ATK_datagroup_has_view(sub, "i0"));

#ifdef XXX
  // we expect the actual data  pointers to be the same
  EXPECT_EQ(ATK_datagroup_get_view(flds, "i0")->getNode().data_pointer(),
            ATK_datagroup_get_group("sub")->get_view("i0")->getNode().data_pointer());
#endif

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * flds = ATK_datagroup_create_group(root, "fields");

  ATK_datagroup * ga = ATK_datagroup_create_group(flds, "a");
  ATK_datagroup * gb = ATK_datagroup_create_group(flds, "b");
  ATK_datagroup * gc = ATK_datagroup_create_group(flds, "c");

  ATK_dataview * tmpview;
  tmpview = ATK_datagroup_create_view_and_buffer_simple(ga, "i0");
  ATK_dataview_allocate(tmpview, ATK_C_INT_T, 1);
  tmpview = ATK_datagroup_create_view_and_buffer_simple(gb, "f0");
  ATK_dataview_allocate(tmpview, ATK_C_FLOAT_T, 1);
  tmpview = ATK_datagroup_create_view_and_buffer_simple(gc, "d0");
  ATK_dataview_allocate(tmpview, ATK_C_DOUBLE_T, 1);

#ifdef XXX
  (*ATK_datagroup_get_view(ga, "i0")->getNode().as_int_ptr())   = 1;
  (*ATK_datagroup_get_view(gb, "f0")->getNode().as_float_ptr()) = 100.0;
  (*ATK_datagroup_get_view(gc, "d0")->getNode().as_double_ptr()) = 3000.0;
#endif
  {
    tmpview = ATK_datagroup_get_view(ga, "i0");
    int * v = (int *) ATK_dataview_get_data_pointer(tmpview);
    EXPECT_TRUE(v != NULL);
    if (v != NULL)
    {
      *v = 1;
    }
  }
  {
    tmpview = ATK_datagroup_get_view(gb, "f0");
    float * v = (float *) ATK_dataview_get_data_pointer(tmpview);
    EXPECT_TRUE(v != NULL);
    if (v != NULL)
    {
      *v = 100.0;
    }
  }
  {
    tmpview = ATK_datagroup_get_view(gc, "d0");
    double * v = (double *) ATK_dataview_get_data_pointer(tmpview);
    EXPECT_TRUE(v != NULL);
    if (v != NULL)
    {
      *v = 3000.0;
    }
  }

  // check that all sub groups exist
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

  //move "b" to a child of "sub"
  ATK_datagroup * sub = ATK_datagroup_create_group(flds, "sub");
  ATK_datagroup_move_group(sub, gb);

  ATK_datagroup_print(flds);

  EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "sub"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

  ATK_datagroup * tmpgrp = ATK_datagroup_get_group(flds, "sub");
  EXPECT_EQ(ATK_datagroup_get_group(tmpgrp, "b"), gb);

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  ATK_datastore * const ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * const grp = ATK_datagroup_create_group(root, "grp");

  const char * viewName1 = "viewBuffer1";
  const char * viewName2 = "viewBuffer2";

  // XXX const
  ATK_dataview * view1 = ATK_datagroup_create_view_and_buffer_simple(grp, viewName1);
  ATK_dataview * view2 = ATK_datagroup_create_view_and_buffer_simple(grp, viewName2);

  EXPECT_TRUE(ATK_datagroup_has_view(grp, viewName1));
  EXPECT_EQ( ATK_datagroup_get_view(grp, viewName1), view1 );

  EXPECT_TRUE(ATK_datagroup_has_view(grp, viewName2));
  EXPECT_EQ( ATK_datagroup_get_view(grp, viewName2), view2 );

  ATK_databuffer * tmpbuf = ATK_dataview_get_buffer(view1);
  ATK_IndexType bufferId1 = ATK_databuffer_get_index(tmpbuf);

  ATK_datagroup_destroy_view_and_buffer(grp, viewName1);


  EXPECT_FALSE(ATK_datagroup_has_view(grp, viewName1));
  EXPECT_EQ(ATK_datastore_get_num_buffers(ds), 1u);

  ATK_databuffer * buffer1 = ATK_datastore_get_buffer(ds, bufferId1);
  bool buffValid = true;
  if( buffer1 == NULL )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_alloc_view_and_buffer)
{
  ATK_datastore * const ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * const grp = ATK_datagroup_create_group(root, "grp");

  const char * viewName1 = "viewBuffer1";
  const char * viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  ATK_dataview * const view1 = ATK_datagroup_create_view_and_buffer_from_type
                                 (grp, viewName1, ATK_C_INT_T, 10);

  // this one is the Schema & method
  ATK_dataview * const view2 = ATK_datagroup_create_view_and_buffer_from_type
                                 (grp, viewName2, ATK_C_DOUBLE_T, 10);

  EXPECT_TRUE(ATK_datagroup_has_view(grp, viewName1));
  EXPECT_EQ( ATK_datagroup_get_view(grp, viewName1), view1 );

  EXPECT_TRUE(ATK_datagroup_has_view(grp, viewName2));
  EXPECT_EQ( ATK_datagroup_get_view(grp, viewName2), view2 );

#ifdef XXX
  int * v1_vals = (int *) ATK_dataview_get_data(view1);
  double * v2_vals = (double *)  ATK_dataview_get_data(view1);

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
    v2_vals[i] = i * 3.1415;
  }
#endif

  EXPECT_EQ(ATK_dataview_get_number_of_elements(view1), 10u);
  EXPECT_EQ(ATK_dataview_get_number_of_elements(view2), 10u);
  EXPECT_EQ(ATK_dataview_get_total_bytes(view1), 10 * sizeof(int));
  EXPECT_EQ(ATK_dataview_get_total_bytes(view2), 10 * sizeof(double));

  ATK_datagroup_destroy_view_and_buffer(grp, viewName1);
  ATK_datagroup_destroy_view_and_buffer(grp, viewName2);

  ATK_datastore_delete(ds);
}

#ifdef XXX
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_view_of_buffer_with_schema)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  // use create + alloc convenience methods
  // this one is the DataType & method
  ATK_dataview * base =  ATK_datagroup_create_view_and_buffer_from_type(root, "base",
                                                                        ATK_C_INT_T, 10);
#ifdef XXX
  int * base_vals = (int *) ATK_dataview_get_data(base);
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

  ATK_databuffer * base_buff = ATK_dataview_get_buffer(base);
  // create two views into this buffer
  // view for the first 5 values
  ATK_datagroup_createView(root, "sub_a", base_buff, ATK_C_INT_T, 5);
  // view for the second 5 values
  //  (schema call path case)
  Schema s(DataType::uint32(5, 5*sizeof(int)));
  ATK_datagroup_createView(root, "sub_b", base_buff, s);

  int * sub_a_vals = (int *) ATK_dataview_get_data(ATK_datagroup_get_view(root, "sub_a"));
  int * sub_b_vals = (int *) ATK_dataview_get_data(ATK_datagroup_get_view(root, "sub_b"));

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
    EXPECT_EQ(sub_b_vals[i], 20);
  }

  ATK_datastore_delete(ds);
}
#endif

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_simple)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * flds = ATK_datagroup_create_group(root, "fields");

  ATK_datagroup * ga = ATK_datagroup_create_group(flds, "a");

  ATK_dataview_allocate(ATK_datagroup_create_view_and_buffer_simple(ga, "i0"), ATK_C_INT_T, 1);

#ifdef XXX
  int * ival = (int *) ATK_dataview_get_data(ATK_datagroup_get_view(ga, "i0"));
  *ival = 1;
#endif

  EXPECT_TRUE(ATK_datagroup_has_group(root, "fields"));
  EXPECT_TRUE(ATK_datagroup_has_group(ATK_datagroup_get_group(root, "fields"), "a"));
  EXPECT_TRUE(ATK_datagroup_has_view(ATK_datagroup_get_group(ATK_datagroup_get_group(root, "fields"), "a"), "i0"));


  ATK_datagroup_save(root, "out_sidre_group_save_restore_simple","conduit");

  ATK_datastore_print(ds);

  ATK_datastore * ds2 = ATK_datastore_new();

  ATK_datagroup_load(ATK_datastore_get_root(ds2), "out_sidre_group_save_restore_simple","conduit");

  ATK_datastore_print(ds2);

  flds = ATK_datagroup_get_group(ATK_datastore_get_root(ds2), "fields");

  // check that all sub groups exist
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
#ifdef XXX
  EXPECT_EQ(ATK_datagroup_get_group(flds, "a")->get_view("i0")->getNode().as_int32(),1);
#endif

  ATK_datastore_print(ds2);

  ATK_datastore_delete(ds);
  ATK_datastore_delete(ds2);

}

#ifdef XXX
//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * flds = ATK_datagroup_create_group(root, "fields");

  ATK_datagroup * ga = ATK_datagroup_create_group(flds, "a");
  ATK_datagroup * gb = ATK_datagroup_create_group(flds, "b");
  ATK_datagroup * gc = ATK_datagroup_create_group(flds, "c");

  ATK_dataview * tmpview;
  tmpview = ATK_datagroup_create_view_and_buffer_simple(ga, "i0");
  ATK_dataview_allocate(tmpview, ATK_C_INT_T, 1);
  int * ival = (int *) ATK_dataview_get_data(tmpview);
  *ival = 1;

  tmpview = ATK_datagroup_create_view_and_buffer_simple(gb, "f0");
  ATK_dataview_allocate(tmpview, ATK_C_FLOAT_T, 1);
  float * fval = (float *) ATK_dataview_get_data(tmpview);
  *fval = 100.0;

  tmpview = ATK_datagroup_create_view_and_buffer_simple(gc, "d0");
  ATK_dataview_allocate(tmpview, ATK_C_DOUBLE_T, 1);
  double * dval = (double *) ATK_dataview_get_data(tmpview);
  *dval = 3000.0;

  // check that all sub groups exist
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

  ATK_datagroup_save(root, "out_sidre_group_save_restore_complex","conduit");

  ATK_datastore_print(ds);

  ATK_datastore * ds2 = ATK_datastore_new();
  root = ATK_datastore_get_root(ds2);

  ATK_datagroup_load(root, "out_sidre_group_save_restore_complex","conduit");

  flds = ATK_datagroup_get_group(root, "fields");
  // check that all sub groups exist
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
  EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

#ifdef XXX
  EXPECT_EQ(ATK_datagroup_get_group(flds, "a")->get_view("i0")->getNode().as_int(),1);
  EXPECT_NEAR(ATK_datagroup_get_group(flds, "b")->get_view("f0")->getNode().as_float(),100.0,  1e-12);
  EXPECT_NEAR(ATK_datagroup_get_group(flds, "c")->get_view("d0")->getNode().as_double(),3000.0, 1e-12);
#endif

  ATK_datastore_print(ds2);

  ATK_datastore_delete(ds);
  ATK_datastore_delete(ds2);

}
#endif

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
