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
#include <cstring>

#include "sidre/sidre.hpp"

using asctoolkit::sidre::SidreLength;
using asctoolkit::sidre::TypeID;
using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
using asctoolkit::sidre::IndexType;
using asctoolkit::sidre::InvalidIndex;
using asctoolkit::sidre::nameIsValid;
using asctoolkit::sidre::indexIsValid;
using asctoolkit::sidre::DataType;

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
}

//------------------------------------------------------------------------------
// getNameWithPath()
//------------------------------------------------------------------------------
TEST(sidre_group,get_name_with_path)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * group =
    root->createGroup("test1")->createGroup("test2")->createGroup("test3");
  DataGroup * group2 = root->getGroup("test1/test2/test3");

  EXPECT_EQ(group, group2);

  // Now verify that code will not create missing groups.
  // TODO - improve error handling so this isn't fatal.
//  DataGroup * group3 = root->createGroup("testa")->createGroup("testb")->createGroup("testc");
//  DataGroup * group_bad = root->getGroup("testa/BAD/testc");

//  (void)group3;

//  EXPECT_EQ(group_bad, root->getGroup("testa") );

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
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(sidre_group,get_group)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataGroup * child = parent->createGroup("child");
  EXPECT_TRUE( child->getParent() == parent );

  EXPECT_TRUE( parent->getGroup("child") == child );
  // check error condition
  EXPECT_TRUE( parent->getGroup("non-existant group") == ATK_NULLPTR );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataView * view = parent->createView("view");

  EXPECT_TRUE( parent->getView("view") == view );

  // check error condition
  EXPECT_TRUE( parent->getView("non-existant view") == ATK_NULLPTR );

  delete ds;
}
//------------------------------------------------------------------------------
// Verify getViewWithPath()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view_with_path)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataView * view =
    root->createGroup("group1")->createGroup("group2")->createView("view1");
  DataView * view2 = root->getView("group1/group2/view1");

  EXPECT_EQ( view, view2 );

  delete ds;
}


//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view_names_and_indicies)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataView * view1 = parent->createView("view1");
  DataView * view2 = parent->createView("view2");

  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getViewIndex("view1");
  IndexType idx2 = parent->getViewIndex("view2");

  const std::string& name1 = parent->getViewName(idx1);
  const std::string& name2 = parent->getViewName(idx2);

  EXPECT_EQ(name1, std::string("view1"));
  EXPECT_EQ(view1->getName(), name1);

  EXPECT_EQ(name2, std::string("view2"));
  EXPECT_EQ(view2->getName(), name2);

  // check error conditions
  IndexType idx3 = parent->getViewIndex("view3");
  EXPECT_TRUE(idx3 == InvalidIndex);

  const std::string& name3 = parent->getViewName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(nameIsValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_view_index)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  DataGroup * parent = root->createGroup("parent");
  DataView * view1 = parent->createView("view1");
  DataView * view2 = parent->createView("view2");

  DataGroup * emptyGroup = root->createGroup("emptyGroup");

  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getFirstValidViewIndex();
  IndexType idx2 = parent->getNextValidViewIndex(idx1);

  const std::string& name1 = parent->getViewName(idx1);
  const std::string& name2 = parent->getViewName(idx2);

  EXPECT_EQ(name1, std::string("view1"));
  EXPECT_EQ(view1->getName(), name1);

  EXPECT_EQ(name2, std::string("view2"));
  EXPECT_EQ(view2->getName(), name2);

  // check error conditions
  IndexType badidx1 = emptyGroup->getFirstValidViewIndex();
  IndexType badidx2 = emptyGroup->getNextValidViewIndex(badidx1);

  EXPECT_TRUE(badidx1 == InvalidIndex);
  EXPECT_TRUE(badidx2 == InvalidIndex);

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

  // check error conditions
  IndexType idx3 = parent->getGroupIndex("group3");
  EXPECT_TRUE(idx3 == InvalidIndex);

  const std::string& name3 = parent->getGroupName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(nameIsValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// createView()
// createViewAndAllocate()
// destroyView()
// destroyViewAndData()
// hasView()
//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_has_view)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  DataGroup * group = root->createGroup("parent");

  DataView * view = group->createView("view");
  EXPECT_TRUE( group->getParent() == root );
  EXPECT_FALSE( view->hasBuffer() );

  EXPECT_TRUE( group->hasView("view") );
  // try creating view again, should be a no-op.
  EXPECT_TRUE( group->createView("view") == ATK_NULLPTR );

  group->destroyView("view");
  // destroy already destroyed group.  Should be a no-op, not a failure
  group->destroyView("view");

  EXPECT_FALSE( group->hasView("view") );

  // try api call that specifies specific type and length
  group->createViewAndAllocate( "viewWithLength1",
                                asctoolkit::sidre::FLOAT_ID, 50 );

  // error condition check - try again with duplicate name, should be a no-op
  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLength1",
                                             asctoolkit::sidre::FLOAT64_ID,
                                             50 ) == ATK_NULLPTR );
  group->destroyViewAndData("viewWithLength1");
  EXPECT_FALSE( group->hasView("viewWithLength1") );

  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLengthBadLen",
                                             asctoolkit::sidre::FLOAT64_ID,
                                             -1 ) == ATK_NULLPTR );

  // try api call that specifies data type in another way
  group->createViewAndAllocate( "viewWithLength2", DataType::float64(50) );
  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLength2",
                                             DataType::float64(
                                               50) ) == ATK_NULLPTR );
  // destroy this view using index
  group->destroyViewAndData( group->getFirstValidViewIndex() );

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

  // should be a no-op, not a failure
  root->destroyGroup("group");

  DataGroup * group2 = root->createGroup("group2");
  // shut up compiler about unused variable
  (void)group2;
  root->destroyGroup( root->getFirstValidGroupIndex() );

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,group_name_collisions)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");
  flds->createView("a");

  EXPECT_TRUE(flds->hasView("a"));

  // attempt to create duplicate group name

  DataGroup * badGroup = ds->getRoot()->createGroup("fields");
  EXPECT_TRUE( badGroup == ATK_NULLPTR );

  // check error condition
  // attempt to create duplicate view name.
  EXPECT_TRUE(flds->createView("a") == ATK_NULLPTR);

  delete ds;
}
//------------------------------------------------------------------------------
#if 0
TEST(sidre_group,view_copy_move)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  flds->createViewAndAllocate("i0", DataType::c_int());
  flds->createViewAndAllocate("f0", DataType::c_float());
  flds->createViewAndAllocate("d0", DataType::c_double());

  flds->getView("i0")->setScalar(1);
  flds->getView("f0")->setScalar(100.0);
  flds->getView("d0")->setScalar(3000.0);

  EXPECT_TRUE(flds->hasView("i0"));
  EXPECT_TRUE(flds->hasView("f0"));
  EXPECT_TRUE(flds->hasView("d0"));

  // test moving a view from flds to sub
  flds->createGroup("sub")->moveView(flds->getView("d0"));
  // flds->print();
  EXPECT_FALSE(flds->hasView("d0"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->getGroup("sub")->hasView("d0"));

  // check the data value
  double * d0_data =  flds->getGroup("sub")->getView("d0")->getData();
  EXPECT_NEAR(d0_data[0],3000.0,1e-12);

  // test copying a view from flds to sub
  flds->getGroup("sub")->copyView(flds->getView("i0"));

  // flds->print();

  EXPECT_TRUE(flds->hasView("i0"));
  EXPECT_TRUE(flds->getGroup("sub")->hasView("i0"));

  // we expect the data pointers to be the same
  int * i0_ptr = flds->getView("i0")->getData();
  int * sub_io0_ptr = flds->getGroup("sub")->getView("i0")->getData();
  EXPECT_EQ(i0_ptr, sub_io0_ptr);

  delete ds;
}
#endif
//------------------------------------------------------------------------------
#if 0
TEST(sidre_group,groups_move_copy)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds->createGroup("b");
  DataGroup * gc = flds->createGroup("c");

  ga->createViewAndAllocate("i0", DataType::c_int());
  gb->createViewAndAllocate("f0", DataType::c_float());
  gc->createViewAndAllocate("d0", DataType::c_double());

  ga->getView("i0")->setScalar(1);
  gb->getView("f0")->setScalar(100.0);
  gc->getView("d0")->setScalar(3000.0);

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  // move "b" to a child of "sub"
  flds->createGroup("sub")->moveGroup(gb);

  // flds->print();

  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->hasGroup("c"));

  EXPECT_EQ(flds->getGroup("sub")->getGroup("b"),gb);

  delete ds;
}
#endif
//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_view_and_buffer2)
{
  DataStore * const ds = new DataStore();
  DataGroup * const grp = ds->getRoot()->createGroup("grp");

  std::string viewName1("viewBuffer1");
  std::string viewName2("viewBuffer2");

  DataView * view1 = grp->createViewAndAllocate(viewName1,
                                                asctoolkit::sidre::INT_ID, 1);
  DataView * view2 = grp->createViewAndAllocate(viewName2,
                                                asctoolkit::sidre::INT_ID, 1);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndData(viewName1);

  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE( buffer1 == ATK_NULLPTR );

  DataView const * const view3 = grp->createView("viewBuffer3");
  grp->destroyViewsAndData();
  // should be no-op
  grp->destroyViewsAndData();
  // shut up compiler about unused variable
  (void)view3;

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
  DataView * const view1 = grp->createViewAndAllocate(viewName1,
                                                      DataType::c_int(10));

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  int * v1_vals = view1->getData();

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
  }

  EXPECT_EQ(view1->getNumElements(), 10u);
  EXPECT_EQ(view1->getTotalBytes(),
            static_cast<asctoolkit::sidre::SidreLength>(10 * sizeof(int)));

  grp->destroyViewAndData(viewName1);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,create_view_of_buffer_with_schema)
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();
  // use create + alloc convenience methods
  // this one is the DataType & method
  DataView * base =  root->createViewAndAllocate("base", DataType::c_int(10));
  int * base_vals = base->getData();
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
  //
  // view for the first 5 values
  root->createView("sub_a", base_buff)->apply(DataType::c_int(5));

  int * sub_a_vals = root->getView("sub_a")->getData();

  for(int i=0 ; i<5 ; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_empty)
{
  DataStore * ds = new DataStore();

  ds->save("sidre_empty_datastore", "conduit");

  delete ds;

  ds = new DataStore();
  ds->load("sidre_empty_datastore", "conduit");

  EXPECT_TRUE(ds->getNumBuffers() == 0 );
  EXPECT_TRUE(ds->getRoot()->getNumGroups() == 0 );
  EXPECT_TRUE(ds->getRoot()->getNumViews() == 0 );

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_api)
{
  DataStore * ds1 = new DataStore();

  ds1->getRoot()->createViewScalar<int>("i0", 1);

  // All three of these should be produce identical files.
  // (I've diff'ed the output files manually.)

  // Call from datastore, no group, defaults to root group
  ds1->save("sidre_save_api1", "conduit");
  // Call from group
  ds1->getRoot()->save("sidre_save_api2", "conduit");
  // Call from datastore, pass in group to save
  ds1->save("sidre_save_api3", "conduit", ds1->getRoot());

  // Text output ( for debugging )
  ds1->save("sidre_save_text_output", "text", ds1->getRoot());

  DataStore * ds2 = new DataStore();
  DataStore * ds3 = new DataStore();
  DataStore * ds4 = new DataStore();

  ds2->load("sidre_save_api1", "conduit");
  ds3->load("sidre_save_api2", "conduit", ds3->getRoot() );
  ds4->getRoot()->load("sidre_save_api3", "conduit");

  EXPECT_TRUE( ds1->getRoot()->isEquivalentTo(ds2->getRoot()) );
  EXPECT_TRUE( ds2->getRoot()->isEquivalentTo(ds3->getRoot()) );
  EXPECT_TRUE( ds3->getRoot()->isEquivalentTo(ds4->getRoot()) );

  delete ds1;
  delete ds2;
  delete ds3;
  delete ds4;

  // Why don't these pass??? Need to ask Noah about this...
  // Trying to make sure sub trees are same here.
#if 0
  DataStore * ds_new = new DataStore();
  DataGroup * api1 = ds_new->getRoot()->createGroup("api1");
  DataGroup * api2 = ds_new->getRoot()->createGroup("api2");
  DataGroup * api3 = ds_new->getRoot()->createGroup("api3");

  api1->load("sidre_save_api1", "conduit");
  api2->load("sidre_save_api2", "conduit");
  api3->load("sidre_save_api3", "conduit");

  EXPECT_TRUE( api1->isEquivalentTo( api2) );
  EXPECT_TRUE( api2->isEquivalentTo( api3) );
#endif

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_scalars_and_strings)
{
  DataStore * ds1 = new DataStore();

  ds1->getRoot()->createViewScalar<int>("i0", 1);
  ds1->getRoot()->createViewScalar<float>("f0", 1.0);
  ds1->getRoot()->createViewScalar<double>("d0", 10.0);

  ds1->save("sidre_save_scalars_and_strings", "conduit");

  DataStore * ds2 = new DataStore();

  ds2->load("sidre_save_scalars_and_strings", "conduit");

  EXPECT_TRUE( ds1->getRoot()->isEquivalentTo( ds2->getRoot()) );

  delete ds1;
  delete ds2;
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_external_data)
{
  int foo[100];
  for (int i =0; i < 100; ++i)
  {
    foo[i] = i;
  }

  DataStore * ds = new DataStore();

  DataView * view = ds->getRoot()->createView("external_array", &foo[0] );
  view->apply( asctoolkit::sidre::INT_ID, 100 );

  ds->save("sidre_save_external", "conduit");

  // Text output ( for debugging )
  ds->save("sidre_save_text_output2", "text");

  delete ds;

  // Now load back in.
  ds = new DataStore();
  ds->load("sidre_save_external", "conduit");

  // All this code should change after we re-write how we handle restoring external data.
  // Right now, the external data is coming back in the view's node and we have to do extra work to
  // restore it to the user's pointer.
  view = ds->getRoot()->getView("external_array");

  int* new_data_pointer = new int[100];

  EXPECT_TRUE( view->getTotalBytes() == sizeof( int[100] ) );

  std::memcpy(&new_data_pointer[0], view->getVoidPtr(), view->getTotalBytes() );

  // Will set view back to EMPTY and reset node.  Will leave description alone.
  view->setExternalDataPtr( ATK_NULLPTR );

  view->setExternalDataPtr( new_data_pointer );

  for (int i = 0; i < 100; ++i)
  {
    EXPECT_TRUE( static_cast<int*>( view->getVoidPtr() )[i] == i );
  }

  delete[] new_data_pointer;
  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_complex)
{
  DataStore * ds = new DataStore();
  DataGroup * flds = ds->getRoot()->createGroup("fields");

  DataGroup * ga = flds->createGroup("a");
  DataGroup * gb = flds->createGroup("b");
  DataGroup * gc = flds->createGroup("c");

  ga->createViewScalar<int>("i0", 100.0);
  ga->createViewScalar<double>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int100", asctoolkit::sidre::INT_ID, 100);
  int* data_ptr = gc->getView("int100")->getArray();
  for (int i =0; i < 100; ++i)
  {
    data_ptr[i] = i;
  }

  ds->save("sidre_mixed_types","conduit");

  DataStore * ds2 = new DataStore();

  ds2->load("sidre_mixed_types","conduit");

  EXPECT_TRUE( ds->getRoot()->isEquivalentTo(ds2->getRoot()) );

  delete ds;

  flds = ds2->getRoot()->getGroup("fields");

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  EXPECT_EQ(flds->getGroup("a")->getView("i0")->getData<int>(),100.0);
  EXPECT_NEAR(flds->getGroup("a")->getView("d0")->getData<double>(),3000.0, 1e-12);

  int* new_data_ptr = flds->getGroup("c")->getView("int100")->getArray();
  for (int i = 0; i < 100; ++i)
  {
    EXPECT_TRUE( new_data_ptr[i] == i);
  }

  // TODO - Figure out the right way to get the string value our of conduit node!!
  //char * char_ptr = flds->getGroup("b")->getView("s0")->getString();
  //EXPECT_TRUE( std::string(char_ptr) == "foo" );



  //ds2->print();

  delete ds2;
  }

//------------------------------------------------------------------------------
// isEquivalentTo()
//------------------------------------------------------------------------------
TEST(sidre_group,is_equivalent_to)
{
  DataStore * ds = new DataStore();

  //These are the parents for two separate subtrees of the root group.
  //Everything below them will be created identically.
  DataGroup * parent1 = ds->getRoot()->createGroup("parent1");
  DataGroup * parent2 = ds->getRoot()->createGroup("parent2");

  //The flds1 and flds2 groups will be compared for equivalence
  DataGroup * flds1 = parent1->createGroup("fields");
  DataGroup * flds2 = parent2->createGroup("fields");

  DataGroup * ga1 = flds1->createGroup("a");
  DataGroup * gb1 = flds1->createGroup("b");
  DataGroup * gc1 = flds1->createGroup("c");
  DataGroup * ga2 = flds2->createGroup("a");
  DataGroup * gb2 = flds2->createGroup("b");
  DataGroup * gc2 = flds2->createGroup("c");

  ga1->createViewScalar("i0", 1 );
  gb1->createViewScalar("f0", 100.0f );
  gc1->createViewScalar("d0", 3000.00);
  ga2->createViewScalar("i0", 1);
  gb2->createViewScalar("f0", 100.0f);
  gc2->createViewScalar("d0", 3000.00);

  ga1->getView("i0")->setScalar(1);
  gb1->getView("f0")->setScalar( 100.0f );
  gc1->getView("d0")->setScalar(3000.00);
  ga2->getView("i0")->setScalar(1);
  gb2->getView("f0")->setScalar( 100.0f );
  gc2->getView("d0")->setScalar(3000.00);

  // Groups were created identically, so should be equivalent.
  EXPECT_TRUE(flds1->isEquivalentTo(flds2));

  // Add something extra to flds2, making them not equivalent.
  gc2->createViewAndAllocate("extra", DataType::c_double());

  EXPECT_FALSE(flds1->isEquivalentTo(flds2));

  delete ds;

}
