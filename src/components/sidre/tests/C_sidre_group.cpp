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
#if 0

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
 
    ATK_datagroup *parent = root->createGroup("parent");
    ATK_datagroup *child = parent->createGroup("child");
    EXPECT_TRUE( child->getParent() == parent );

    EXPECT_TRUE( parent->hasGroup("child") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify hasView()
//------------------------------------------------------------------------------
TEST(C_sidre_group,has_view)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = root->createGroup("parent");
    DataView *view = parent->createViewAndBuffer("view");

    EXPECT_TRUE( view->getOwningGroup() == parent );

    EXPECT_TRUE( parent->hasView("view") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_view_name_index)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = root->createGroup("parent");
    DataView *view1 = parent->createViewAndBuffer("view1");
    DataView *view2 = parent->createViewAndBuffer("view2");

    EXPECT_EQ(parent->getNumViews(), 2u);

    IDType idx1 = parent->getViewIndex("view1");
    IDType idx2 = parent->getViewIndex("view2");

    std::string name1(parent->getViewName(idx1));
    std::string name2(parent->getViewName(idx2));
   
    EXPECT_EQ(name1, std::string("view1")); 
    EXPECT_EQ(view1->getName(), name1); 

    EXPECT_EQ(name2, std::string("view2")); 
    EXPECT_EQ(view2->getName(), name2); 

#if 0 // Leave out for now until we resolve error/warning/assert macro usage
    IDType idx3 = parent->getViewIndex("view3");
    std::string name3(parent->getViewName(idx3));

    EXPECT_EQ(idx3, InvalidID);
    EXPECT_TRUE(name3.empty());
#endif

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(C_sidre_group,get_group_name_index)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);

    ATK_datagroup *parent = root->createGroup("parent");
    ATK_datagroup *group1 = parent->createGroup("group1");
    ATK_datagroup *group2 = parent->createGroup("group2");

    EXPECT_EQ(parent->getNumGroups(), 2u);

    IDType idx1 = parent->getGroupIndex("group1");
    IDType idx2 = parent->getGroupIndex("group2");

    std::string name1(parent->getGroupName(idx1));
    std::string name2(parent->getGroupName(idx2));

    EXPECT_EQ(name1, std::string("group1"));
    EXPECT_EQ(group1->getName(), name1);

    EXPECT_EQ(name2, std::string("group2"));
    EXPECT_EQ(group2->getName(), name2);

#if 0 // Leave out for now until we resolve error/warning/assert macro usage
    IDType idx3 = parent->getGroupIndex("group3");
    std::string name3(parent->getGroupName(idx3));

    EXPECT_EQ(idx3, InvalidID);
    EXPECT_TRUE(name3.empty());
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
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = root->createGroup("parent");

    DataView *view = group->createViewAndBuffer("view");
    EXPECT_TRUE( group->getParent() == root );
    EXPECT_TRUE( view->hasBuffer() );
 
    EXPECT_TRUE( group->hasView("view") );

    group->destroyViewAndBuffer("view");

    EXPECT_FALSE( group->hasView("view") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_has_group)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *root = ATK_datastore_get_root(ds);
    ATK_datagroup *group = root->createGroup("group");
    EXPECT_TRUE( group->getParent() == root );
    
    EXPECT_TRUE( root->hasGroup("group") );


    root->destroyGroup("group");
    EXPECT_FALSE( root->hasGroup("group") );

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *flds = ATK_datastore_get_root(ds)->createGroup("fields");
    flds->createViewAndBuffer("a");

    EXPECT_TRUE(flds->hasView("a"));

    ATK_datastore_delete(ds);
}
//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *flds = ATK_datastore_get_root(ds)->createGroup("fields");

    flds->createViewAndBuffer("i0")->allocate(DataType::int32());
    flds->createViewAndBuffer("f0")->allocate(DataType::float32());
    flds->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*flds->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*flds->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*flds->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    EXPECT_TRUE(flds->hasView("i0"));
    EXPECT_TRUE(flds->hasView("f0"));
    EXPECT_TRUE(flds->hasView("d0"));

    // test moving a view form feds7 to sub
    flds->createGroup("sub")->moveView(flds->getView("d0"));
    flds->print();
    EXPECT_FALSE(flds->hasView("d0"));
    EXPECT_TRUE(flds->hasGroup("sub"));
    EXPECT_TRUE(flds->getGroup("sub")->hasView("d0"));

    // check the data value
    float64 *d0_data =  flds->getGroup("sub")
                            ->getView("d0")
                            ->getNode().as_float64_ptr();
    EXPECT_NEAR(d0_data[0],3000.0,1e-12);
    
    // test copying a view from flds top sub
    flds->getGroup("sub")->copyView(flds->getView("i0"));

    flds->print();
    
    EXPECT_TRUE(flds->hasView("i0"));    
    EXPECT_TRUE(flds->getGroup("sub")->hasView("i0"));

    // we expect the actual data  pointers to be the same
    EXPECT_EQ(flds->getView("i0")->getNode().data_pointer(),
              flds->getGroup("sub")->getView("i0")->getNode().data_pointer());

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *flds = ATK_datastore_get_root(ds)->createGroup("fields");

    ATK_datagroup *ga = flds->createGroup("a");
    ATK_datagroup *gb = flds->createGroup("b");
    ATK_datagroup *gc = flds->createGroup("c");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());
    gb->createViewAndBuffer("f0")->allocate(DataType::float32());
    gc->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*gb->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*gc->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    //move "b" to a child of "sub"
    flds->createGroup("sub")->moveGroup(gb);

    flds->print();
    
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("sub"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("sub")->getGroup("b"),gb);

    ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  ATK_datastore * const ds = ATK_datastore_new();
  ATK_datagroup * const grp = ATK_datastore_get_root(ds)->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  DataView const * const view1 = grp->createViewAndBuffer(viewName1);
  DataView const * const view2 = grp->createViewAndBuffer(viewName2);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  IDType const bufferId1 = view1->getBuffer()->getUID();

  grp->destroyViewAndBuffer(viewName1);


  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
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
  ATK_datagroup * const grp = ATK_datastore_get_root(ds)->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  DataView * const view1 = grp->createViewAndBuffer(viewName1,
                                                          DataType::uint32(10));
  // this one is the Schema & method
  Schema s;
  s.set(DataType::float64(10));
  DataView * const view2 = grp->createViewAndBuffer(viewName2,
                                                          s);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  
  uint32  *v1_vals = view1->getNode().as_uint32_ptr();
  float64 *v2_vals = view2->getNode().as_float64_ptr();
  
  for(int i=0;i<10;i++)
  {
      v1_vals[i] = i;
      v2_vals[i] = i * 3.1415;
  }


  EXPECT_EQ(view1->getSchema().total_bytes(), 10 * sizeof(uint32));
  EXPECT_EQ(view2->getSchema().total_bytes(), 10 * sizeof(float64));
    
  grp->destroyViewAndBuffer(viewName1);
  grp->destroyViewAndBuffer(viewName2);

  ATK_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,create_view_of_buffer_with_schema)
{
  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);
  // use create + alloc convenience methods
  // this one is the DataType & method
  DataView *  base =  root->createViewAndBuffer("base",
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

  DataBuffer *base_buff = base->getBuffer();
  // create two views into this buffer
  // view for the first 5 values
  root->createView("sub_a",base_buff,DataType::uint32(5));
  // view for the second 5 values
  //  (schema call path case)
  Schema s(DataType::uint32(5,5*sizeof(uint32)));
  root->createView("sub_b",base_buff,s);

  uint32 *sub_a_vals = root->getView("sub_a")->getNode().as_uint32_ptr();
  uint32 *sub_b_vals = root->getView("sub_b")->getNode().as_uint32_ptr();

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
    ATK_datagroup *flds = ATK_datastore_get_root(ds)->createGroup("fields");

    ATK_datagroup *ga = flds->createGroup("a");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;

    EXPECT_TRUE(ATK_datastore_get_root(ds)->hasGroup("fields"));
    EXPECT_TRUE(ATK_datastore_get_root(ds)->getGroup("fields")->hasGroup("a"));
    EXPECT_TRUE(ATK_datastore_get_root(ds)->getGroup("fields")->getGroup("a")->hasView("i0"));
        
        
    ATK_datastore_get_root(ds)->save("out_sidre_group_save_restore_simple","conduit");

    ds->print();
    
    ATK_datastore *ds2 = ATK_datastore_new();

    ds2->getRoot()->load("out_sidre_group_save_restore_simple","conduit");
    
    ds2->print();

    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    
    ds2->print();
    
    ATK_datastore_delete(ds);    
    ATK_datastore_delete(ds2);    
    
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
    ATK_datastore *ds = ATK_datastore_new();
    ATK_datagroup *flds = ATK_datastore_get_root(ds)->createGroup("fields");

    ATK_datagroup *ga = flds->createGroup("a");
    ATK_datagroup *gb = flds->createGroup("b");
    ATK_datagroup *gc = flds->createGroup("c");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());
    gb->createViewAndBuffer("f0")->allocate(DataType::float32());
    gc->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*gb->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*gc->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    ATK_datastore_get_root(ds)->save("out_sidre_group_save_restore_complex","conduit");

    ds->print();

    ATK_datastore *ds2 = ATK_datastore_new();


    ds2->getRoot()->load("out_sidre_group_save_restore_complex","conduit");

    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    EXPECT_NEAR(flds->getGroup("b")->getView("f0")->getNode().as_float32(),100.0,  1e-12);
    EXPECT_NEAR(flds->getGroup("c")->getView("d0")->getNode().as_float64(),3000.0, 1e-12);

    ds2->print();

    ATK_datastore_delete(ds);
    ATK_datastore_delete(ds2);

}


#endif
