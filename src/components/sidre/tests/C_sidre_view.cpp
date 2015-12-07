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

//------------------------------------------------------------------------------

TEST(C_sidre_view,create_views)
{
  SIDRE_datastore * ds   = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_dataview * dv_0 = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "field0",
                                                      SIDRE_INT_ID, 1);
  SIDRE_dataview * dv_1 = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "field1",
                                                      SIDRE_INT_ID, 1);

  SIDRE_databuffer * db_0 = SIDRE_dataview_get_buffer(dv_0);
  SIDRE_databuffer * db_1 = SIDRE_dataview_get_buffer(dv_1);

  EXPECT_EQ(SIDRE_databuffer_get_index(db_0), 0);
  EXPECT_EQ(SIDRE_databuffer_get_index(db_1), 1);
  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_buffer_from_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_dataview * dv = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "u0",
                                                      SIDRE_INT_ID, 10);

  EXPECT_EQ(SIDRE_dataview_get_type_id(dv), SIDRE_INT_ID);
  int * data_ptr = (int *) SIDRE_dataview_get_data_pointer(dv);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  SIDRE_dataview_print(dv);

  EXPECT_EQ(SIDRE_dataview_get_total_bytes(dv), sizeof(int) * 10);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_buffer_from_view_conduit_value)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  SIDRE_dataview * dv = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "u0", 
                                                      SIDRE_INT_ID, 10);
  int * data_ptr = (int *) SIDRE_dataview_get_data_pointer(dv);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i*i;
  }

  SIDRE_dataview_print(dv);

  EXPECT_EQ(SIDRE_dataview_get_total_bytes(dv), sizeof(int) * 10);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_array_multi_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_databuffer * dbuff = SIDRE_datastore_create_buffer(ds);

  SIDRE_databuffer_declare(dbuff, SIDRE_INT_ID, 10);
  SIDRE_databuffer_allocate_existing(dbuff);
  int * data_ptr = (int *) SIDRE_databuffer_get_data(dbuff);

  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = i;
  }

  SIDRE_databuffer_print(dbuff);

  EXPECT_EQ(SIDRE_databuffer_get_total_bytes(dbuff), sizeof(int) * 10);


  SIDRE_dataview * dv_e = 
     SIDRE_datagroup_create_view_into_buffer(root, "even", dbuff);
  SIDRE_dataview * dv_o = 
     SIDRE_datagroup_create_view_into_buffer(root, "odd", dbuff);
  EXPECT_TRUE(dv_e != NULL);
  EXPECT_TRUE(dv_o != NULL);
  EXPECT_EQ(SIDRE_databuffer_get_num_views(dbuff), 2u);

  SIDRE_dataview_apply_nelems_offset_stride(dv_e, 5, 0, 2);
  SIDRE_dataview_apply_nelems_offset_stride(dv_o, 5, 1, 2);

  SIDRE_dataview_print(dv_e);
  SIDRE_dataview_print(dv_o);

// Note: This is a big hack since the dataview get pointer method is broken
//       and the conduit support for this sort of thing doesn't exist for C code?
  int* dv_e_ptr = (int *) SIDRE_dataview_get_data_pointer(dv_e);
  int* dv_o_ptr = (int *) SIDRE_dataview_get_data_pointer(dv_o); dv_o_ptr++;
  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(dv_e_ptr[2*i], 2*i);
    EXPECT_EQ(dv_o_ptr[2*i], 2*i+1);
  }

  // Run similar test to above with different view apply method
  SIDRE_dataview * dv_e1 =
     SIDRE_datagroup_create_view_into_buffer(root, "even1", dbuff);
  SIDRE_dataview * dv_o1 =
     SIDRE_datagroup_create_view_into_buffer(root, "odd1", dbuff);
  EXPECT_TRUE(dv_e1 != NULL);
  EXPECT_TRUE(dv_o1 != NULL);
  EXPECT_EQ(SIDRE_databuffer_get_num_views(dbuff), 4u);

  SIDRE_dataview_apply_type_nelems_offset_stride(dv_e1, SIDRE_INT_ID, 5, 0, 2);
  SIDRE_dataview_apply_type_nelems_offset_stride(dv_o1, SIDRE_INT_ID, 5, 1, 2);

  SIDRE_dataview_print(dv_e1);
  SIDRE_dataview_print(dv_o1);

// Note: This is a big hack since the dataview get pointer method is broken
//       and the conduit support for this sort of thing doesn't exist for C code?
  int* dv_e1_ptr = (int *) SIDRE_dataview_get_data_pointer(dv_e1);
  int* dv_o1_ptr = (int *) SIDRE_dataview_get_data_pointer(dv_o1); dv_o1_ptr++;
  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(dv_e1_ptr[2*i], 2*i);
    EXPECT_EQ(dv_o1_ptr[2*i], 2*i+1);
  }


  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_array_depth_view)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_databuffer * dbuff = SIDRE_datastore_create_buffer(ds);

  const size_t depth_nelems = 10;

  SIDRE_databuffer_declare(dbuff, SIDRE_INT_ID, 4 * depth_nelems);
  SIDRE_databuffer_allocate_existing(dbuff);
  int * data_ptr = (int *) SIDRE_databuffer_get_data(dbuff);

  for(size_t i=0 ; i < 4 * depth_nelems ; i++)
  {
    data_ptr[i] = i / depth_nelems;
  }

  SIDRE_databuffer_print(dbuff);

  EXPECT_EQ(SIDRE_databuffer_get_num_elements(dbuff), 4 * depth_nelems);

  // create 4 "depth" views and apply offsets into buffer
  SIDRE_dataview* views[4];
  const char* view_names[4] = { "depth_0", "depth_1", "depth_2", "depth_3" };

  for (int id = 0; id < 4; ++id)
  {
     views[id] = SIDRE_datagroup_create_view_into_buffer(root, view_names[id], dbuff);
     SIDRE_dataview_apply_nelems_offset(views[id], depth_nelems, id*depth_nelems);
  }
  EXPECT_EQ(SIDRE_databuffer_get_num_views(dbuff), 4u);

  // print depth views...
  for (int id = 0; id < 4; ++id)
  {
     SIDRE_dataview_print(views[id]);
  }

  // check values in depth views...
  for (int id = 0; id < 4; ++id)
  {
// Note: This is a big hack since the dataview get pointer method is broken
//       and the conduit support for this sort of thing doesn't exist for C code?
     int* dv_ptr = (int *) SIDRE_dataview_get_data_pointer(views[id]); 
          dv_ptr += id * depth_nelems; 
     for (size_t i = 0; i < depth_nelems; ++i)
     {
        EXPECT_EQ(dv_ptr[i], id);
     }
  }

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

#if 0
// Similar to previous test, using other view creation methods
TEST(C_sidre_view,int_array_depth_view_2)
{
  SIDRE_datastore * ds = SIDRE_datastore_new();
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  SIDRE_databuffer * dbuff = SIDRE_datastore_create_buffer(ds);

  const size_t depth_nelems = 10;

  SIDRE_databuffer_declare(dbuff, SIDRE_INT_ID, 4 * depth_nelems);
  SIDRE_databuffer_allocate_existing(dbuff);
  int * data_ptr = (int *) SIDRE_databuffer_get_data(dbuff);

  for(size_t i=0 ; i < 4 * depth_nelems ; i++)
  {
    data_ptr[i] = i / depth_nelems;
  }

  SIDRE_databuffer_print(dbuff);

  EXPECT_EQ(SIDRE_databuffer_get_num_elements(dbuff), 4 * depth_nelems);

  // create 4 "depth" views and apply offsets into buffer
  SIDRE_dataview* views[4];
  const char* view_names[4] = { "depth_0", "depth_1", "depth_2", "depth_3" };

  for (int id = 0; id < 2; ++id)
  {
     views[id] = 
        SIDRE_datagroup_create_view_into_buffer_nelems_offset(root, 
           view_names[id], dbuff, depth_nelems, id*depth_nelems);
  }
  for (int id = 2; id < 4; ++id)
  {
     views[id] = 
        SIDRE_datagroup_create_view_into_buffer_type_nelems_offset(root, 
           view_names[id], dbuff, SIDRE_INT_ID, depth_nelems, id*depth_nelems);
  }
  EXPECT_EQ(SIDRE_databuffer_get_num_views(dbuff), 4u);

  // print depth views...
  for (int id = 0; id < 4; ++id)
  {
     SIDRE_dataview_print(views[id]);
  }

  // check values in depth views...
  for (int id = 0; id < 4; ++id)
  {
// Note: This is a big hack since the dataview get pointer method is broken
//       and the conduit support for this sort of thing doesn't exist for C code?
     int* dv_ptr = (int *) SIDRE_dataview_get_data_pointer(views[id]); 
          dv_ptr += id * depth_nelems; 
     for (size_t i = 0; i < depth_nelems; ++i)
     {
        EXPECT_EQ(dv_ptr[i], id);
     }
  }

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);

}
#endif

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_array_multi_view_resize)
{
  ///
  /// This example creates a 4 * 10 buffer of ints,
  /// and 4 views that point the 4 sections of 10 ints
  ///
  /// We then create a new buffer to support 4*12 ints
  /// and 4 views that point into them
  ///
  /// after this we use the old buffers to copy the values
  /// into the new views
  ///

  // create our main data store
  SIDRE_datastore * ds = SIDRE_datastore_new();
  // get access to our root data Group
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  // create a group to hold the "old" or data we want to copy
  SIDRE_datagroup * r_old = SIDRE_datagroup_create_group(root, "r_old");
  // create a view to hold the base buffer
  SIDRE_dataview * base_old = 
     SIDRE_datagroup_create_view_and_allocate_from_type(r_old, "base_data",
                                                      SIDRE_INT_ID, 40);

  int * data_ptr = (int *) SIDRE_dataview_get_data_pointer(base_old);


  // init the buff with values that align with the
  // 4 subsections.
  for(int i=0 ; i<10 ; i++)
  {
    data_ptr[i] = 1;
  }
  for(int i=10 ; i<20 ; i++)
  {
    data_ptr[i] = 2;
  }
  for(int i=20 ; i<30 ; i++)
  {
    data_ptr[i] = 3;
  }
  for(int i=30 ; i<40 ; i++)
  {
    data_ptr[i] = 4;
  }

#ifdef XXX
  /// setup our 4 views
  SIDRE_databuffer * buff_old = SIDRE_dataview_get_buffer(base_old);
  buff_old->getNode().print();
  SIDRE_dataview * r0_old = SIDRE_dataview_create_view(r_old, "r0",buff_old);
  SIDRE_dataview * r1_old = SIDRE_dataview_create_view(r_old, "r1",buff_old);
  SIDRE_dataview * r2_old = SIDRE_dataview_create_view(r_old, "r2",buff_old);
  SIDRE_dataview * r3_old = SIDRE_dataview_create_view(r_old, "r3",buff_old);

  // each view is offset by 10 * the # of bytes in a uint32
  // uint32(num_elems, offset)
  index_t offset =0;
  r0_old->apply(r0_old, DataType::uint32(10,offset));

  offset += sizeof(int) * 10;
  r1_old->apply(r1_old, DataType::uint32(10,offset));

  offset += sizeof(int) * 10;
  r2_old->apply(r2_old, DataType::uint32(10,offset));

  offset += sizeof(int) * 10;
  r3_old->apply(r3_old, DataType::uint32(10,offset));

  /// check that our views actually point to the expected data
  //
  uint32 * r0_ptr = r0_old->getNode().as_uint32_ptr();
  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r0_ptr[i], 1u);
    // check pointer relation
    EXPECT_EQ(&r0_ptr[i], &data_ptr[i]);
  }

  uint32 * r3_ptr = r3_old->getNode().as_uint32_ptr();
  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r3_ptr[i], 4u);
    // check pointer relation
    EXPECT_EQ(&r3_ptr[i], &data_ptr[i+30]);
  }

  // create a group to hold the "old" or data we want to copy into
  SIDRE_datagroup * r_new = SIDRE_datagroup_create_group(root, "r_new");
  // create a view to hold the base buffer
  SIDRE_dataview * base_new = SIDRE_datagroup_create_view_and_buffer_simple(r_new, "base_data");

  // alloc our buffer
  // create a buffer to hold larger subarrays
  base_new->allocate_from_type(base_new, DataType::uint32(4 * 12));
  int * base_new_data = (int *) SIDRE_databuffer_det_data(base_new);
  for (int i = 0 ; i < 4 * 12 ; ++i)
  {
    base_new_data[i] = 0;
  }

  SIDRE_databuffer * buff_new = SIDRE_dataview_get_buffer(base_new);
  buff_new->getNode().print();

  // create the 4 sub views of this array
  SIDRE_dataview * r0_new = SIDRE_datagroup_create_view(r_new, "r0",buff_new);
  SIDRE_dataview * r1_new = SIDRE_datagroup_create_view(r_new, "r1",buff_new);
  SIDRE_dataview * r2_new = SIDRE_datagroup_create_view(r_new, "r2",buff_new);
  SIDRE_dataview * r3_new = SIDRE_datagroup_create_view(r_new, "r3",buff_new);

  // apply views to r0,r1,r2,r3
  // each view is offset by 12 * the # of bytes in a uint32

  // uint32(num_elems, offset)
  offset =0;
  r0_new->apply(DataType::uint32(12,offset));

  offset += sizeof(int) * 12;
  r1_new->apply(DataType::uint32(12,offset));

  offset += sizeof(int) * 12;
  r2_new->apply(DataType::uint32(12,offset));

  offset += sizeof(int) * 12;
  r3_new->apply(DataType::uint32(12,offset));

  /// update r2 as an example first
  buff_new->getNode().print();
  r2_new->getNode().print();

  /// copy the subset of value
  r2_new->getNode().update(r2_old->getNode());
  r2_new->getNode().print();
  buff_new->getNode().print();


  /// check pointer values
  int * r2_new_ptr = (int *) SIDRE_dataview_get_data_pointer(r2_new);

  for(int i=0 ; i<10 ; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 3);
  }

  for(int i=10 ; i<12 ; i++)
  {
    EXPECT_EQ(r2_new_ptr[i], 0);     // assumes zero-ed alloc
  }


  /// update the other views
  r0_new->getNode().update(r0_old->getNode());
  r1_new->getNode().update(r1_old->getNode());
  r3_new->getNode().update(r3_old->getNode());

  buff_new->getNode().print();
#endif

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_view,int_array_realloc)
{
  ///
  /// info
  ///

  // create our main data store
  SIDRE_datastore * ds = SIDRE_datastore_new();
  // get access to our root data Group
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);

  // create a view to hold the base buffer
  SIDRE_dataview * a1 = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "a1", 
                                                      SIDRE_FLOAT_ID, 5);
  SIDRE_dataview * a2 = 
     SIDRE_datagroup_create_view_and_allocate_from_type(root, "a2", 
                                                      SIDRE_INT_ID, 5);

  float * a1_ptr = (float *) SIDRE_dataview_get_data_pointer(a1);
  int * a2_ptr = (int *)  SIDRE_dataview_get_data_pointer(a2);

  for(int i=0 ; i<5 ; i++)
  {
    a1_ptr[i] =  5.0;
    a2_ptr[i] = -5;
  }

  EXPECT_EQ(SIDRE_dataview_get_total_bytes(a1), sizeof(float)*5);
  EXPECT_EQ(SIDRE_dataview_get_total_bytes(a2), sizeof(int)*5);


  SIDRE_dataview_reallocate(a1, 10);
  SIDRE_dataview_reallocate(a2, 15);

  a1_ptr = (float *) SIDRE_dataview_get_data_pointer(a1);
  a2_ptr = (int *) SIDRE_dataview_get_data_pointer(a2);

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(a1_ptr[i],5.0);
    EXPECT_EQ(a2_ptr[i],-5);
  }

  for(int i=5 ; i<10 ; i++)
  {
    a1_ptr[i] = 10.0;
    a2_ptr[i] = -10;
  }

  for(int i=10 ; i<15 ; i++)
  {
    a2_ptr[i] = -15;
  }

  EXPECT_EQ(SIDRE_dataview_get_total_bytes(a1), sizeof(float)*10);
  EXPECT_EQ(SIDRE_dataview_get_total_bytes(a2), sizeof(int)*15);


  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);

}

//------------------------------------------------------------------------------

TEST(C_sidre_view,simple_opaque)
{
  // create our main data store
  SIDRE_datastore * ds = SIDRE_datastore_new();
  // get access to our root data Group
  SIDRE_datagroup * root = SIDRE_datastore_get_root(ds);
  int * src_data = (int *) malloc(sizeof(int));

  src_data[0] = 42;

  void * src_ptr = (void *)src_data;

  SIDRE_dataview * opq_view = SIDRE_datagroup_create_opaque_view(root, "my_opaque",src_ptr);

  // we shouldn't have any buffers
  EXPECT_EQ(SIDRE_datastore_get_num_buffers(ds), 0u);

  EXPECT_TRUE(SIDRE_dataview_is_opaque(opq_view));

  void * opq_ptr = SIDRE_dataview_get_opaque(opq_view);

  int * out_data = (int *)opq_ptr;
  EXPECT_EQ(opq_ptr,src_ptr);
  EXPECT_EQ(out_data[0],42);

  SIDRE_datastore_print(ds);
  SIDRE_datastore_delete(ds);
  free(src_data);
}
