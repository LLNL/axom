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
using asctoolkit::sidre::IndexType;
using asctoolkit::sidre::InvalidIndex;
using asctoolkit::sidre::isNameValid;
using asctoolkit::sidre::indexIsValid;
using asctoolkit::sidre::DataType;

using namespace conduit;


// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_name)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * group = root->createGroup("test");

  EXPECT_TRUE(group->getName() == std::string("test") );

  DataGroup * group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == ATK_NULLPTR);

  delete ds;
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(sidre_group,get_parent)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * parent = root->createGroup("parent");
  DataGroup * child = parent->createGroup("child");

  EXPECT_TRUE( child->getParent() == parent );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group,get_datastore)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * group = root->createGroup("parent");

  EXPECT_TRUE( group->getDataStore() == ds );

  DataStore const * const_ds = group->getDataStore();
  EXPECT_TRUE( const_ds == ds );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify hasGroup()
//------------------------------------------------------------------------------
TEST(sidre_group,has_group)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataGroup * child = parent->createGroup("child");
  EXPECT_TRUE( child->getParent() == parent );

  EXPECT_TRUE( parent->hasGroup("child") );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify hasView()
//------------------------------------------------------------------------------
TEST(sidre_group,has_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataView * view = parent->createViewAndBuffer("view");

  EXPECT_TRUE( view->getOwningGroup() == parent );

  EXPECT_TRUE( parent->hasView("view") );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view_name_index)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataView * view1 = parent->createViewAndBuffer("view1");
  DataView * view2 = parent->createViewAndBuffer("view2");

  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getViewIndex("view1");
  IndexType idx2 = parent->getViewIndex("view2");

  const std::string& name1 = parent->getViewName(idx1);
  const std::string& name2 = parent->getViewName(idx2);

  EXPECT_EQ(name1, std::string("view1"));
  EXPECT_EQ(view1->getName(), name1);

  EXPECT_EQ(name2, std::string("view2"));
  EXPECT_EQ(view2->getName(), name2);

  IndexType idx3 = parent->getViewIndex("view3");
  EXPECT_TRUE(idx3 == InvalidIndex);

  const std::string& name3 = parent->getViewName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(isNameValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(sidre_group,get_group_name_index)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataGroup * group1 = parent->createGroup("group1");
  DataGroup * group2 = parent->createGroup("group2");

  EXPECT_EQ(parent->getNumGroups(), 2u);

  IndexType idx1 = parent->getGroupIndex("group1");
  IndexType idx2 = parent->getGroupIndex("group2");

  const std::string& name1 = parent->getGroupName(idx1);
  const std::string& name2 = parent->getGroupName(idx2);

  EXPECT_EQ(name1, std::string("group1"));
  EXPECT_EQ(group1->getName(), name1);

  EXPECT_EQ(name2, std::string("group2"));
  EXPECT_EQ(group2->getName(), name2);

  IndexType idx3 = parent->getGroupIndex("group3");
  EXPECT_TRUE(idx3 == InvalidIndex);

  const std::string& name3 = parent->getGroupName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(isNameValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// createViewAndBuffer()
// destroyViewAndBuffer()
// hasView()
//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_has_viewbuffer)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * group = root->createGroup("parent");

  DataView * view = group->createViewAndBuffer("view");
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
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * group = root->createGroup("group");
  EXPECT_TRUE( group->getParent() == root );

  EXPECT_TRUE( root->hasGroup("group") );


  root->destroyGroup("group");
  EXPECT_FALSE( root->hasGroup("group") );

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,group_name_collisions)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");
  flds->createViewAndBuffer("a");

  EXPECT_TRUE(flds->hasView("a"));

  delete ds;
}
//------------------------------------------------------------------------------
TEST(sidre_group,view_copy_move)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  flds->createViewAndBuffer("i0")->allocate(DataType::c_int());
  flds->createViewAndBuffer("f0")->allocate(DataType::c_float());
  flds->createViewAndBuffer("d0")->allocate(DataType::c_double());

  // TODO - change to setScalar when ready
  (*flds->getView("i0")->getNode().as_int_ptr())   = 1;
  (*flds->getView("f0")->getNode().as_float_ptr()) = 100.0;
  (*flds->getView("d0")->getNode().as_double_ptr()) = 3000.0;

  EXPECT_TRUE(flds->hasView("i0"));
  EXPECT_TRUE(flds->hasView("f0"));
  EXPECT_TRUE(flds->hasView("d0"));

  // test moving a view from flds to sub
  flds->createGroup("sub")->moveView(flds->getView("d0"));
  flds->print();
  EXPECT_FALSE(flds->hasView("d0"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->getGroup("sub")->hasView("d0"));

  // check the data value
  double * d0_data =  flds->getGroup("sub")
                     ->getView("d0")
                     ->getValue();
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds to sub
  flds->getGroup("sub")->copyView(flds->getView("i0"));

  flds->print();

  EXPECT_TRUE(flds->hasView("i0"));
  EXPECT_TRUE(flds->getGroup("sub")->hasView("i0"));

  // we expect the actual data  pointers to be the same
  EXPECT_EQ(flds->getView("i0")->getDataPointer(),
            flds->getGroup("sub")->getView("i0")->getDataPointer());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,groups_move_copy)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds->createGroup("b");
  DataGroup * gc = flds->createGroup("c");

  ga->createViewAndBuffer("i0")->allocate(DataType::c_int());
  gb->createViewAndBuffer("f0")->allocate(DataType::c_float());
  gc->createViewAndBuffer("d0")->allocate(DataType::c_double());

  // TODO - change to setScalar when ready
  (*ga->getView("i0")->getNode().as_int_ptr())   = 1;
  (*gb->getView("f0")->getNode().as_float_ptr()) = 100.0;
  (*gc->getView("d0")->getNode().as_double_ptr()) = 3000.0;

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

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndBuffer(viewName1);


  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE( buffer1 == ATK_NULLPTR );

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
                                                    DataType::c_int(10));
  // this one is the Schema & method
  Schema s;
  s.set(DataType::c_double(10));
  DataView * const view2 = grp->createViewAndBuffer(viewName2,
                                                    s);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );


  int * v1_vals = view1->getValue();
  double * v2_vals = view2->getValue();

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
    v2_vals[i] = i * 3.1415;
  }


  EXPECT_EQ(view1->getNumberOfElements(), 10u);
  EXPECT_EQ(view2->getNumberOfElements(), 10u);
  EXPECT_EQ(view1->getTotalBytes(), 10 * sizeof(int));
  EXPECT_EQ(view2->getTotalBytes(), 10 * sizeof(double));

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
  DataView * base =  root->createViewAndBuffer("base",
                                               DataType::c_int(10));
  int * base_vals = base->getValue();
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

  DataBuffer * base_buff = base->getBuffer();
  // create two views into this buffer
  // view for the first 5 values
  root->createView("sub_a", base_buff, DataType::c_int(5));
  // view for the second 5 values
  //  (schema call path case)
  Schema s(DataType::c_int(5,5*sizeof(int)));
  root->createView("sub_b",base_buff,s);

  int * sub_a_vals = root->getView("sub_a")->getValue();
  int * sub_b_vals = root->getView("sub_b")->getValue();

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
    EXPECT_EQ(sub_b_vals[i], 20);
  }

  delete ds;
}




//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_simple)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  DataGroup * ga = flds->createGroup("a");

  ga->createViewAndBuffer("i0")->allocate(DataType::c_int());

  // TODO - change to setScalar when ready
  (*ga->getView("i0")->getNode().as_int_ptr())   = 1;

  EXPECT_TRUE(ds->getRoot()->hasGroup("fields"));
  EXPECT_TRUE(ds->getRoot()->getGroup("fields")->hasGroup("a"));
  EXPECT_TRUE(ds->getRoot()->getGroup("fields")->getGroup("a")->hasView("i0"));


  ds->getRoot()->save("out_sidre_group_save_restore_simple","conduit");

  ds->print();

  DataStore * ds2 = new DataStore();

  ds2->getRoot()->load("out_sidre_group_save_restore_simple","conduit");

  ds2->print();

  flds = ds2->getRoot()->getGroup("fields");
  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  int testvalue = flds->getGroup("a")->getView("i0")->getValue();
  EXPECT_EQ(testvalue,1);

  ds2->print();

  delete ds;
  delete ds2;

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_complex)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds->createGroup("b");
  DataGroup * gc = flds->createGroup("c");

  ga->createViewAndBuffer("i0")->allocate(DataType::c_int());
  gb->createViewAndBuffer("f0")->allocate(DataType::c_float());
  gc->createViewAndBuffer("d0")->allocate(DataType::c_double());

  // TODO - change to setScalar when ready
  (*ga->getView("i0")->getNode().as_int_ptr())   = 1;
  (*gb->getView("f0")->getNode().as_float_ptr()) = 100.0;
  (*gc->getView("d0")->getNode().as_double_ptr()) = 3000.0;

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  ds->print();

  ds->getRoot()->save("out_sidre_group_save_restore_complex","conduit");

  DataStore * ds2 = new DataStore();


  ds2->getRoot()->load("out_sidre_group_save_restore_complex","conduit");

  flds = ds2->getRoot()->getGroup("fields");
  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  // TODO - change to setScalar when ready
  EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int(),1);
  EXPECT_NEAR(flds->getGroup("b")->getView("f0")->getNode().as_float(),100.0,  1e-12);
  EXPECT_NEAR(flds->getGroup("c")->getView("d0")->getNode().as_double(),3000.0, 1e-12);

  ds2->print();

  delete ds;
  delete ds2;

}

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
