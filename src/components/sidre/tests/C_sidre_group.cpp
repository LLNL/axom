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
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = ATK_datagroup_create_group(root, "test");
 
    //    EXPECT_TRUE(group->getName() == std::string("test") );
    EXPECT_TRUE(strcmp(ATK_datagroup_get_name(group), "test") == 0);

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_parent)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *parent = ATK_datagroup_create_group(root, "parent");
    ATK_datagroup *child = ATK_datagroup_create_group(parent, "child");
 
    EXPECT_TRUE( ATK_datagroup_get_parent(child) == parent );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_datastore)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = ATK_datagroup_create_group(root, "parent");
 
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
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
 
    ATK_datagroup *parent = ATK_datagroup_create_group(root, "parent");
    ATK_datagroup *child = ATK_datagroup_create_group(parent, "child");
    EXPECT_TRUE( ATK_datagroup_get_parent(child) == parent );

    EXPECT_TRUE( ATK_datagroup_has_group(parent, "child") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,has_view)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = ATK_datagroup_create_group(root, "parent");

    ATK_dataview *view = ATK_datagroup_create_view_and_buffer(parent, "view");

    EXPECT_TRUE( ATK_dataview_get_owning_group(view) == parent );

    EXPECT_TRUE( ATK_datagroup_has_view(parent, "view") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
#if 0
TEST(C_sidre_group,get_view_name_index)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = ATK_datagroup_create_group(root, "parent");
    ATK_dataview *view1 = ATK_datagroup_create_view_and_buffer(parent, "view1");
    ATK_dataview *view2 = ATK_datagroup_create_view_and_buffer(parent, "view2");

    EXPECT_EQ(ATK_datagroup_get_num_views(parent), 2u);

    ATK_IndexType idx1 = ATK_datagroup_get_view_index(parent, "view1");
    ATK_IndexType idx2 = ATK_datagroup_get_view_index(parent, "view2");

    char *name1 = ATK_datagroup_get_view_name(parent, idx1);
    char *name2 = ATK_datagroup_get_view_name(parent, idx2);
   
    EXPECT_TRUE(strcmp(name1, "view1") == 0);
    EXPECT_TRUE(strcmp(ATK_dataview_get_name(view1), name1) == 0);

    EXPECT_TRUE(strcmp(name2, "view2") == 0);
    EXPECT_TRUE(strcmp(ATK_dataview_get_name(view2), name2) == 0);

#if 0 // Leave out for now until we resolve error/warning/assert macro usage
    ATK_IndexType idx3 = ATK_datagroup_get_view_index(parent, "view3");
    char *name3 = ATK_datagroup_get_view_name(parent, idx3);

    EXPECT_EQ(idx3, InvalidID);
    EXPECT_TRUE(name3.empty());
#endif

    ATK_datastore_delete(ds);
}
#endif
//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
#if 0
TEST(C_sidre_group,get_group_name_index)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = ATK_datagroup_create_group(root, "parent");
    ATK_datagroup *group1 = ATK_datagroup_create_group(parent, "group1");
    ATK_datagroup *group2 = ATK_datagroup_create_group(parent, "group2");

    EXPECT_EQ(ATK_datagroup_get_num_groups(parent), 2u);

    ATK_IndexType idx1 = ATK_datagroup_get_group_index(parent, "group1");
    ATK_IndexType idx2 = ATK_datagroup_get_group_index(parent, "group2");

    const char *name1 = ATK_datagroup_get_group_name(parent, idx1);
    const char *name2 = ATK_datagroup_get_group_name(parent, idx2);

    EXPECT_TRUE(strcmp(name1, "group1") == 0);
    EXPECT_TRUE(strcmp(ATK_datagroup_get_name(group1), name1) == 0);

    EXPECT_EQ(name2, std::string("group2"));
    EXPECT_EQ(ATK_datagroup_get_name(group1), name2);

#if 0 // Leave out for now until we resolve error/warning/assert macro usage
    ATK_IndexType idx3 = ATK_datagroup_get_group_index(parent, "group3");
    std::string name3(ATK_datagroup_get_group_name(parent, idx3));

    EXPECT_EQ(idx3, InvalidID);
    EXPECT_TRUE(name3.empty());
#endif

    ATK_datastore_delete(ds);
}
#endif

//------------------------------------------------------------------------------
// createViewAndBuffer()
// destroyViewAndBuffer()
// hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_viewbuffer)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = ATK_datagroup_create_group(root, "parent");

    ATK_dataview *view = ATK_datagroup_create_view_and_buffer(group, "view");
    EXPECT_TRUE( ATK_datagroup_get_parent(group) == root );
    EXPECT_TRUE( ATK_dataview_has_buffer(view) );
 
    EXPECT_TRUE( ATK_datagroup_has_view(group, "view") );

    ATK_datagroup_destroy_view_and_buffer(group, "view");

    EXPECT_FALSE( ATK_datagroup_has_view(group, "view") );

    ATK_datastore_delete(ds);
}
#if 0

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_group)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = ATK_datagroup_create_group(root, "group");
    EXPECT_TRUE( ATK_datagroup_get_parent(group, ) == root );
    
    EXPECT_TRUE( ATK_datagroup_has_group(root, "group") );


    ATK_datagroup_destroy_group(root, "group");
    EXPECT_FALSE( ATK_datagroup_has_group(root, "group") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *flds = ATK_datagroup_create_group(root, "fields");
    ATK_datagroup_create_view_and_buffer(flds, "a");

    EXPECT_TRUE(ATK_datagroup_has_view(flds, "a"));

    ATK_datastore_delete(ds);
}
//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *flds = ATK_datagroup_create_group(root, "fields");

    ATK_datagroup_create_view_and_buffer(flds, "i0")->allocate(DataType::int32());
    ATK_datagroup_create_view_and_buffer(flds, "f0")->allocate(DataType::float32());
    ATK_datagroup_create_view_and_buffer(flds, "d0")->allocate(DataType::float64());

    (*ATK_datagroup_get_view(flds, "i0")->getNode().as_int32_ptr())   = 1;
    (*ATK_datagroup_get_view(flds, "f0")->getNode().as_float32_ptr()) = 100.0;
    (*ATK_datagroup_get_view(flds, "d0")->getNode().as_float64_ptr()) = 3000.0;

    EXPECT_TRUE(ATK_datagroup_has_view(flds, "i0"));
    EXPECT_TRUE(ATK_datagroup_has_view(flds, "f0"));
    EXPECT_TRUE(ATK_datagroup_has_view(flds, "d0"));

    // test moving a view form feds7 to sub
    //    ATK_datagroup_createGroup(flds, "sub")->moveView(ATK_datagroup_get_view(flds, "d0"));
    ATK_datagroup *sub = ATK_datagroup_create_group(flds, "sub");
    AT_datagroup_move_view();
    ATK_datagroup_print();
    EXPECT_FALSE(ATK_datagroup_has_view(flds, "d0"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "sub"));
    EXPECT_TRUE(ATK_datagroup_get_group(flds, "sub")->has_view("d0"));

    // check the data value
    float64 *d0_data =  ATK_datagroup_get_group(flds, "sub")
                            ->get_view("d0")
                            ->getNode().as_float64_ptr();
    EXPECT_NEAR(d0_data[0],3000.0,1e-12);
    
    // test copying a view from flds top sub
    ATK_datagroup_get_group(flds, "sub")->copyView(ATK_datagroup_get_view("i0"));

    ATK_datagroup_print();
    
    EXPECT_TRUE(ATK_datagroup_has_view(flds, "i0"));    
    EXPECT_TRUE(ATK_datagroup_get_group(flds, "sub")->has_view("i0"));

    // we expect the actual data  pointers to be the same
    EXPECT_EQ(ATK_datagroup_get_view(flds, "i0")->getNode().data_pointer(),
              ATK_datagroup_get_group("sub")->get_view("i0")->getNode().data_pointer());

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *flds = ATK_datagroup_create_group(root, "fields");

    ATK_datagroup *ga = ATK_datagroup_create_group(flds, "a");
    ATK_datagroup *gb = ATK_datagroup_create_group(flds, "b");
    ATK_datagroup *gc = ATK_datagroup_create_group(flds, "c");

    ATK_datagroup_create_view_and_buffer(ga, "i0")->allocate(DataType::int32());
    ATK_datagroup_create_view_and_buffer(gb, "f0")->allocate(DataType::float32());
    ATK_datagroup_create_view_and_buffer(gc, "d0")->allocate(DataType::float64());

    (*ATK_datagroup_get_view(ga, "i0")->getNode().as_int32_ptr())   = 1;
    (*ATK_datagroup_get_view(gb, "f0")->getNode().as_float32_ptr()) = 100.0;
    (*ATK_datagroup_get_view(gc, "d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

    //move "b" to a child of "sub"
    //    ATK_datagroup_createGroup(flds, "sub")->moveGroup(gb);
    ATK_datagroup *sub = ATK_datagroup_create_group(flds, "sub");
    ATK_datagroup_move_group(sub, gb);

    ATK_datagroup_print();
    
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "sub"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

    EXPECT_EQ(ATK_datagroup_get_group(flds, "sub")->get_group("b"),gb);

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  ATK_datastore * const ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  ATK_datagroup * const grp = ATK_datagroup_create_group(root, "grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  DataView const * const view1 = grp->create_view_and_buffer(grp, viewName1);
  DataView const * const view2 = grp->create_view_and_buffer(grp, viewName2);

  EXPECT_TRUE(grp->has_view(grp, viewName1));
  EXPECT_EQ( grp->get_view(grp, viewName1), view1 );

  EXPECT_TRUE(grp->has_view(grp, viewName2));
  EXPECT_EQ( grp->get_view(grp, viewName2), view2 );

  ATK_IndexType const bufferId1 = view1->getBuffer(view1)->getUID();

  grp->destroy_view_and_buffer(grp, viewName1);


  EXPECT_FALSE(grp->has_view(grp, viewName1));
  EXPECT_EQ(ds->getNumBuffers(ds), 1u);

  DataBuffer const * const buffer1 = ds->getBuffer(ds, bufferId1);
  bool buffValid = true;
  if( buffer1 == ATK_NULLPTR )
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
  ATK_datagroup *root = ATK_datastore_get_root(ds);
  ATK_datagroup * const grp = ATK_datagroup_create_group(root, "grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  ATK_dataview * const view1 = grp->create_view_and_buffer(grp, viewName1,
                                                          DataType::uint32(10));
  // this one is the Schema & method
  Schema s;
  s.set(DataType::float64(10));
  ATk_dataview * const view2 = grp->create_view_and_buffer(grp, viewName2,
                                                          s);

  EXPECT_TRUE(grp->has_view(grp, viewName1));
  EXPECT_EQ( grp->get_view(grp, viewName1), view1 );

  EXPECT_TRUE(grp->has_view(grp, viewName2));
  EXPECT_EQ( grp->get_view(grp, viewName2), view2 );

  
  uint32  *v1_vals = view1->getNode().as_uint32_ptr();
  float64 *v2_vals = view2->getNode().as_float64_ptr();
  
  for(int i=0;i<10;i++)
  {
      v1_vals[i] = i;
      v2_vals[i] = i * 3.1415;
  }


  EXPECT_EQ(view1->getSchema().total_bytes(), 10 * sizeof(uint32));
  EXPECT_EQ(view2->getSchema().total_bytes(), 10 * sizeof(float64));
    
  grp->destroy_view_and_buffer(grp, viewName1);
  grp->destroy_view_and_buffer(grp, viewName2);

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_view_of_buffer_with_schema)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  // use create + alloc convenience methods
  // this one is the DataType & method
  ATK_dataview *  base =  ATK_datagroup_create_view_and_buffer(root, "base",
                                                DataType::uint32(10)); 
  uint32 *base_vals = base->getNode().as_uint32_ptr();
  for(int i=0;i<10;i++)
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

  DataBuffer *base_buff = base->getBuffer(base);
  // create two views into this buffer
  // view for the first 5 values
  ATK_datagroup_createView(root, "sub_a",base_buff,DataType::uint32(5));
  // view for the second 5 values
  //  (schema call path case)
  Schema s(DataType::uint32(5,5*sizeof(uint32)));
  ATK_datagroup_createView(root, "sub_b",base_buff,s);

  uint32 *sub_a_vals = ATK_datagroup_get_view(root, "sub_a")->getNode().as_uint32_ptr();
  uint32 *sub_b_vals = ATK_datagroup_get_view(root, "sub_b")->getNode().as_uint32_ptr();

  for(int i=0;i<5;i++)
  {
     EXPECT_EQ(sub_a_vals[i], 10u);
     EXPECT_EQ(sub_b_vals[i], 20u);
  }

  ATK_datastore_delete(ds);
}




//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_simple)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *flds = ATK_datagroup_create_group(root, "fields");

    ATK_datagroup *ga = ATK_datagroup_create_group(flds, "a");

    ATK_datagroup_create_view_and_buffer(ga, "i0")->allocate(DataType::int32());

    (*ATK_datagroup_get_view(ga, "i0")->getNode().as_int32_ptr())   = 1;

    EXPECT_TRUE(ATK_datastore_get_root(ds)->has_group("fields"));
    EXPECT_TRUE(ATK_datastore_get_root(ds)->get_group("fields")->has_group("a"));
    EXPECT_TRUE(ATK_datastore_get_root(ds)->get_group("fields")->get_group("a")->has_view("i0"));
        
        
    ATK_datastore_get_root(ds)->save("out_sidre_group_save_restore_simple","conduit");

    ds->print();
    
    ATK_datastore *ds2 = ATK_datastore_new();

    ds2->getRoot(ds2)->load("out_sidre_group_save_restore_simple","conduit");
    
    ds2->print();

    flds = ds2->getRoot(ds2)->get_group("fields");
    // check that all sub groups exist
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
    EXPECT_EQ(ATK_datagroup_get_group(flds, "a")->get_view("i0")->getNode().as_int32(),1);
    
    ds2->print();
    
    ATK_datastore_delete(ds);
    ATK_datastore_delete(ds2);
    
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *flds = ATK_datagroup_create_group(root, "fields");

    ATK_datagroup *ga = ATK_datagroup_create_group(flds, "a");
    ATK_datagroup *gb = ATK_datagroup_create_group(flds, "b");
    ATK_datagroup *gc = ATK_datagroup_create_group(flds, "c");

    ATK_datagroup_create_view_and_buffer(ga, "i0")->allocate(DataType::int32());
    ATK_datagroup_create_view_and_buffer(gb, "f0")->allocate(DataType::float32());
    ATK_datagroup_create_view_and_buffer(gc, "d0")->allocate(DataType::float64());

    (*ATK_datagroup_get_view(ga, "i0")->getNode().as_int32_ptr())   = 1;
    (*ATK_datagroup_get_view(gb, "f0")->getNode().as_float32_ptr()) = 100.0;
    (*ATK_datagroup_get_view(gc, "d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

    ATK_datastore_get_root(ds)->save("out_sidre_group_save_restore_complex","conduit");

    ds->print();

    ATK_datastore *ds2 = ATK_datastore_new();


    ds2->getRoot()->load("out_sidre_group_save_restore_complex","conduit");

    flds = ds2->getRoot(ds2)->get_group("fields");
    // check that all sub groups exist
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "a"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "b"));
    EXPECT_TRUE(ATK_datagroup_has_group(flds, "c"));

    EXPECT_EQ(ATK_datagroup_get_group(flds, "a")->get_view("i0")->getNode().as_int32(),1);
    EXPECT_NEAR(ATK_datagroup_get_group(flds, "b")->get_view("f0")->getNode().as_float32(),100.0,  1e-12);
    EXPECT_NEAR(ATK_datagroup_get_group(flds, "c")->get_view("d0")->getNode().as_float64(),3000.0, 1e-12);

    ds2->print();

    ATK_datastore_delete(ds);
    ATK_datastore_delete(ds2);

}


#endif
