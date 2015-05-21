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

#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
using asctoolkit::common::IDType;

using namespace conduit;

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_name)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("test");
 
    EXPECT_TRUE(group->getName() == std::string("test") );

    delete ds;
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(sidre_group,get_parent)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *parent = root->createGroup("parent");
    DataGroup *child = parent->createGroup("child");
 
    EXPECT_TRUE( child->getParent() == parent );

    delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group,get_datastore)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("parent");
 
    EXPECT_TRUE( group->getDataStore() == ds );

    DataStore const * const_ds = group->getDataStore();
    EXPECT_TRUE( const_ds == ds );

    delete ds;
}

//------------------------------------------------------------------------------
// hasGroup()
//------------------------------------------------------------------------------
TEST(sidre_group,has_child)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
 
    DataGroup *parent = root->createGroup("parent");
    DataGroup *child = parent->createGroup("child");
    EXPECT_TRUE( child->getParent() == parent );

    EXPECT_TRUE( parent->hasGroup("child") );

    delete ds;
}

//------------------------------------------------------------------------------
// createViewAndBuffer()
// destroyViewAndBuffer()
// hasView()
//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_has_viewbuffer)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("parent");

    DataView *view = group->createViewAndBuffer("view");
    EXPECT_TRUE( group->getParent() == root );
    EXPECT_TRUE( view->hasBuffer() );
 
    EXPECT_TRUE( group->hasView("view") );

    group->destroyViewAndBuffer("view");

    EXPECT_FALSE( group->hasView("view") );

    delete ds;
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_has_group)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("group");
    EXPECT_TRUE( group->getParent() == root );
    
    EXPECT_TRUE( root->hasGroup("group") );


    root->destroyGroup("group");
    EXPECT_FALSE( root->hasGroup("group") );

    delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,group_name_collisions)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");
    flds->createViewAndBuffer("a");

    EXPECT_TRUE(flds->hasView("a"));

    delete ds;
}
//------------------------------------------------------------------------------
TEST(sidre_group,view_copy_move)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

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

    delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,groups_move_copy)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");
    DataGroup *gb = flds->createGroup("b");
    DataGroup *gc = flds->createGroup("c");

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

    delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_view_and_buffer)
{
  DataStore * const ds = new DataStore();
  DataGroup * const grp = ds->getRoot()->createGroup("grp");

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
  EXPECT_EQ(ds->getNumberOfBuffers(),1u);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
  bool buffValid = true;
  if( buffer1 == ATK_NULLPTR )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);

  delete ds;
}


//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_alloc_view_and_buffer)
{
  DataStore * const ds = new DataStore();
  DataGroup * const grp = ds->getRoot()->createGroup("grp");

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

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,create_view_of_buffer_with_schema)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
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
     EXPECT_EQ(sub_a_vals[i],10);
     EXPECT_EQ(sub_b_vals[i],20);
  }

  delete ds;
}




//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_simple)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;

    EXPECT_TRUE(ds->getRoot()->hasGroup("fields"));
    EXPECT_TRUE(ds->getRoot()->getGroup("fields")->hasGroup("a"));
    EXPECT_TRUE(ds->getRoot()->getGroup("fields")->getGroup("a")->hasView("i0"));
        
        
    ds->getRoot()->save("out_ds_group_save_restore_simple","conduit");
    
    ds->print();
    
    DataStore *ds2 = new DataStore();

    ds2->getRoot()->load("out_ds_group_save_restore_simple","conduit");
    
    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    
    ds2->print();
    
    delete ds;    
    delete ds2;
    
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_complex)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");
    DataGroup *gb = flds->createGroup("b");
    DataGroup *gc = flds->createGroup("c");

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

    ds->getRoot()->save("out_ds_group_save_restore_complex","conduit");

    ds->print();

    DataStore *ds2 = new DataStore();


    ds2->getRoot()->load("out_ds_group_save_restore_complex","conduit");

    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    EXPECT_NEAR(flds->getGroup("b")->getView("f0")->getNode().as_float32(),100.0,  1e-12);
    EXPECT_NEAR(flds->getGroup("c")->getView("d0")->getNode().as_float64(),3000.0, 1e-12);

    ds2->print();

    delete ds;
    delete ds2;

}


