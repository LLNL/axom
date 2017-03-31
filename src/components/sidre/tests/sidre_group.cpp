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

using axom::sidre::SidreLength;
using axom::sidre::TypeID;
using axom::sidre::Buffer;
using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::DataView;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::nameIsValid;
using axom::sidre::indexIsValid;
using axom::sidre::DataType;
using axom::sidre::INT_ID;
using axom::sidre::FLOAT64_ID;

// Test protocols
int nprotocols = 3;
std::string const protocols[] = { "sidre_json", "sidre_hdf5", "json" };

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_name)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
  Group * group = root->createGroup("test");

  EXPECT_TRUE(group->getName() == std::string("test") );

  Group * group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == AXOM_NULLPTR);
}

//------------------------------------------------------------------------------
// getPath(), getPathName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_path_name)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
  EXPECT_EQ(root->getParent(), root);
  EXPECT_EQ(root->getName(), "");
  Group * group = root->createGroup("test/a/b/c");
  Group * grp2 = root->getGroup("test/a");
  Group * grp3 = root->getGroup("test");

  EXPECT_EQ(root->getName(), std::string(""));
  EXPECT_EQ(root->getPath(), std::string(""));
  EXPECT_EQ(root->getPathName(), std::string(""));

  EXPECT_EQ(grp2->getName(), std::string("a") );
  EXPECT_EQ(grp2->getPath(), std::string("test") );
  EXPECT_EQ(grp2->getPathName(), std::string("test/a") );

  EXPECT_EQ(grp3->getName(), std::string("test") );
  EXPECT_EQ(grp3->getPath(), std::string("") );
  EXPECT_EQ(grp3->getPathName(), std::string("test") );

  EXPECT_EQ(group->getName(), std::string("c") );
  EXPECT_EQ(group->getPath(), std::string("test/a/b") );
  EXPECT_EQ(group->getPathName(), std::string("test/a/b/c") );
}

//------------------------------------------------------------------------------
// createGroup(), getGroup(), hasGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group,group_with_path)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  // Test full path access when building incrementally
  Group * group =
    root->createGroup("test1")->createGroup("test2")->createGroup("test3");
  Group * group2 = root->getGroup("test1/test2/test3");

  EXPECT_TRUE(AXOM_NULLPTR != group2);
  EXPECT_EQ(group, group2);

  // Test incremental access when building full path
  Group * groupP = root->createGroup("testA/testB/testC");
  Group * groupP2 =
    root->getGroup("testA")->getGroup("testB")->getGroup("testC");

  EXPECT_TRUE(AXOM_NULLPTR != groupP2);
  EXPECT_EQ(groupP, groupP2);
  // test non-const getGroup() with path
  Group * groupPParent = root->getGroup("testA/testB");
  EXPECT_EQ(groupP->getParent(), groupPParent);
  EXPECT_EQ(groupP->getParent()->getName(), "testB");


  // Now verify that code will not create missing groups.

  root->createGroup("testa")->createGroup("testb")->createGroup("testc");
  Group * group_bada = root->getGroup("BAD/testb/testc");
  Group * group_badb = root->getGroup("testa/BAD/testc");
  Group * group_badc = root->getGroup("testa/testb/BAD");

  EXPECT_EQ(group_bada, static_cast<void *>(AXOM_NULLPTR) );
  EXPECT_EQ(group_badb, static_cast<void *>(AXOM_NULLPTR) );
  EXPECT_EQ(group_badc, static_cast<void *>(AXOM_NULLPTR) );

  // Test hasGroup with paths.

  EXPECT_FALSE(root->hasGroup("BAD/testb/testc"));
  EXPECT_FALSE(root->hasGroup("testa/BAD/testc"));
  EXPECT_FALSE(root->hasGroup("testa/testb/BAD"));

  EXPECT_TRUE(root->hasGroup("test1"));
  EXPECT_TRUE(root->hasGroup("test1/test2"));
  EXPECT_TRUE(root->hasGroup("test1/test2/test3"));
  Group * group_testa = root->getGroup("testa");
  EXPECT_TRUE(group_testa->hasGroup("testb"));
  EXPECT_TRUE(group_testa->hasGroup("testb/testc"));
  EXPECT_FALSE(group_testa->hasGroup("testb/BAD"));
  EXPECT_FALSE(group_testa->hasGroup("testb/testc/BAD"));

  unsigned int testbnumgroups = group_testa->getGroup("testb")->getNumGroups();
  Group * group_cdup = group_testa->createGroup("testb/testc");

  EXPECT_EQ(group_cdup, static_cast<void *>(AXOM_NULLPTR));
  EXPECT_EQ(group_testa->getGroup("testb")->getNumGroups(), testbnumgroups);

  delete ds;

}

//------------------------------------------------------------------------------
// createGroup(), destroyGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group,destroy_group_with_path)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  // Test full path access when building incrementally
  Group * group = root->createGroup("test1/test2/test3");
  (void)group;

  const std::size_t exp_no_groups = 0;
  const std::size_t exp_one_group = 1;

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("test1")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("test1/test2")->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1/test2/test3")->getNumGroups());

  root->destroyGroup("test1/test2");

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1")->getNumGroups());
  EXPECT_FALSE(root->hasGroup("test1/test2/test3"));
  EXPECT_FALSE(root->hasGroup("test1/test2"));

  root->destroyGroup("test1/BAD");

  EXPECT_EQ(exp_one_group, root->getNumGroups());
  EXPECT_EQ(exp_no_groups, root->getGroup("test1")->getNumGroups());

  delete ds;
}


//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(sidre_group,get_parent)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
  Group * parent = root->createGroup("parent");
  Group * child = parent->createGroup("child");

  EXPECT_TRUE( child->getParent() == parent );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group,get_datastore)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
  Group * group = root->createGroup("parent");

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
  Group * root = ds->getRoot();

  Group * parent = root->createGroup("parent");
  Group * child = parent->createGroup("child");
  EXPECT_TRUE( child->getParent() == parent );

  EXPECT_TRUE( parent->getGroup("child") == child );
  // check error condition
  EXPECT_TRUE( parent->getGroup("non-existant group") == AXOM_NULLPTR );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  Group * parent = root->createGroup("parent");
  DataView * view = parent->createView("view");

  EXPECT_TRUE( parent->getView("view") == view );

  // check error condition
  EXPECT_TRUE( parent->getView("non-existant view") == AXOM_NULLPTR );

  delete ds;
}

//------------------------------------------------------------------------------
// createView, hasView(), getView(), destroyView() with path strings
//------------------------------------------------------------------------------
TEST(sidre_group,view_with_path)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  // Test with full path access when building incrementally
  DataView * view =
    root->createGroup("group1")->createGroup("group2")->createView("view1");
  DataView * view2 = root->getView("group1/group2/view1");

  EXPECT_TRUE(AXOM_NULLPTR != view2);
  EXPECT_EQ( view, view2 );


  // Test incremental access when building with full path
  DataView * viewP = root->createView("groupA/groupB/viewA");
  DataView * viewP2 =
    root->getGroup("groupA")->getGroup("groupB")->getView("viewA");

  EXPECT_TRUE(AXOM_NULLPTR != viewP2);
  EXPECT_EQ( viewP, viewP2 );

  // Now verify that bad paths just return null, and don't create missing groups
  DataView * v_bad1 = root->getView("BAD/groupB/viewA");
  DataView * v_bad2 = root->getView("groupA/BAD/viewA");
  DataView * v_bad3 = root->getView("groupA/groupB/BAD");

  EXPECT_EQ(v_bad1, static_cast<void *>(AXOM_NULLPTR));
  EXPECT_EQ(v_bad2, static_cast<void *>(AXOM_NULLPTR));
  EXPECT_EQ(v_bad3, static_cast<void *>(AXOM_NULLPTR));

  const std::size_t exp_no_groups = 0;
  const std::size_t exp_one_group = 1;
  const std::size_t exp_two_group = 2;


  EXPECT_EQ(exp_two_group, root->getNumGroups());
  EXPECT_TRUE(root->hasGroup("group1"));
  EXPECT_TRUE(root->hasGroup("groupA"));
  EXPECT_EQ(exp_one_group, root->getGroup("group1")->getNumGroups());
  EXPECT_TRUE(root->hasGroup("group1/group2"));
  EXPECT_EQ(exp_no_groups, root->getGroup("group1/group2")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("group1/group2")->getNumViews());
  EXPECT_EQ(root->getGroup("group1/group2")->getView("view1"), view);
  EXPECT_EQ(root->getGroup("group1")->getView("group2/view1"), view);

  EXPECT_EQ(exp_one_group, root->getGroup("groupA")->getNumGroups());
  EXPECT_TRUE(root->hasGroup("groupA/groupB"));
  EXPECT_EQ(exp_no_groups, root->getGroup("groupA/groupB")->getNumGroups());
  EXPECT_EQ(exp_one_group, root->getGroup("groupA/groupB")->getNumViews());
  EXPECT_EQ(root->getGroup("groupA/groupB")->getView("viewA"), viewP);
  EXPECT_EQ(root->getGroup("groupA")->getView("groupB/viewA"), viewP);

  root->destroyView("group1/group2/view1");

  EXPECT_EQ(exp_no_groups, root->getGroup("group1/group2")->getNumViews());
  EXPECT_FALSE(root->getGroup("group1/group2")->hasView("view1"));
  EXPECT_EQ(root->getGroup("group1/group2")->getView("view1"),
            static_cast<void *>(AXOM_NULLPTR));
  EXPECT_FALSE(root->hasView("group1/group2/view1"));
  EXPECT_EQ(root->getView("group1/group2/view1"),
            static_cast<void *>(AXOM_NULLPTR));

  Group * groupA = root->getGroup("groupA");
  EXPECT_TRUE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), viewP);
  EXPECT_TRUE(root->hasView("groupA/groupB/viewA"));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), viewP);

  groupA->destroyView("groupB/viewA");

  EXPECT_EQ(exp_no_groups, groupA->getGroup("groupB")->getNumViews());
  EXPECT_FALSE(groupA->getGroup("groupB")->hasView("viewA"));
  EXPECT_EQ(groupA->getGroup("groupB")->getView("viewA"),
            static_cast<void *>(AXOM_NULLPTR));
  EXPECT_FALSE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), static_cast<void *>(AXOM_NULLPTR));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"),
            static_cast<void *>(AXOM_NULLPTR));

  delete ds;
}


//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view_names_and_indicies)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  Group * parent = root->createGroup("parent");
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
  Group * root = ds->getRoot();

  Group * parent = root->createGroup("parent");
  DataView * view1 = parent->createView("view1");
  DataView * view2 = parent->createView("view2");

  Group * emptyGroup = root->createGroup("emptyGroup");

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
  Group * root = ds->getRoot();

  Group * parent = root->createGroup("parent");
  Group * group1 = parent->createGroup("group1");
  Group * group2 = parent->createGroup("group2");

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
  Group * root = ds->getRoot();
  Group * group = root->createGroup("parent");

  DataView * view = group->createView("view");
  EXPECT_TRUE( group->getParent() == root );
  EXPECT_FALSE( view->hasBuffer() );
  EXPECT_TRUE( group->hasView("view") );
  IndexType iview = group->getViewIndex("view");
  EXPECT_EQ(0, iview);

  // try creating view again, should be a no-op.
  EXPECT_TRUE( group->createView("view") == AXOM_NULLPTR );

  // Create another view to make sure destroyView only destroys one view
  group->createView("viewfiller");
  EXPECT_EQ(2u, group->getNumViews());
  IndexType iviewfiller = group->getViewIndex("viewfiller");
  EXPECT_EQ(1, iviewfiller);

  group->destroyView("view");
  EXPECT_EQ(1u, group->getNumViews());
  // Check if index changed
  EXPECT_EQ(iviewfiller, group->getViewIndex("viewfiller"));

  // destroy already destroyed group.  Should be a no-op, not a failure
  group->destroyView("view");
  EXPECT_EQ(1u, group->getNumViews());
  EXPECT_FALSE( group->hasView("view") );

  // try api call that specifies specific type and length
  group->createViewAndAllocate( "viewWithLength1",INT_ID, 50 );
  IndexType iview2 = group->getViewIndex("viewWithLength1");
  EXPECT_EQ(iview, iview2);  // reuse slot

  // error condition check - try again with duplicate name, should be a no-op
  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLength1", FLOAT64_ID,
                                             50 ) == AXOM_NULLPTR );
  group->destroyViewAndData("viewWithLength1");
  EXPECT_FALSE( group->hasView("viewWithLength1") );

  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLengthBadLen", FLOAT64_ID,
                                             -1 ) == AXOM_NULLPTR );

  // try api call that specifies data type in another way
  group->createViewAndAllocate( "viewWithLength2", DataType::float64(50) );
  EXPECT_TRUE( group->createViewAndAllocate( "viewWithLength2",
                                             DataType::float64(
                                               50) ) == AXOM_NULLPTR );
  // destroy view and its buffer using index
  IndexType indx = group->getFirstValidViewIndex();
  IndexType bindx = group->getView( indx )->getBuffer()->getIndex();
  group->destroyViewAndData( indx );
  EXPECT_TRUE( ds->getBuffer(bindx) == NULL );

  // Destroy view but not the buffer
  view = group->createViewAndAllocate( "viewWithLength2", INT_ID, 50 );
  Buffer * buff = view->getBuffer();
  group->destroyView("viewWithLength2");
  EXPECT_TRUE( buff->isAllocated() );

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
  Group * root = ds->getRoot();
  Group * group = root->createGroup("group");
  EXPECT_TRUE( group->getParent() == root );

  EXPECT_TRUE( root->hasGroup("group") );

  root->destroyGroup("group");
  EXPECT_FALSE( root->hasGroup("group") );

  // should be a no-op, not a failure
  root->destroyGroup("group");

  Group * group2 = root->createGroup("group2");
  // shut up compiler about unused variable
  (void)group2;
  root->destroyGroup( root->getFirstValidGroupIndex() );

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,group_name_collisions)
{
  DataStore * ds = new DataStore();
  Group * flds = ds->getRoot()->createGroup("fields");
  flds->createView("a");

  EXPECT_TRUE(flds->hasChildView("a"));

  // attempt to create duplicate group name
  EXPECT_TRUE( ds->getRoot()->createGroup("fields") == AXOM_NULLPTR );

  // attempt to create duplicate view name.
  EXPECT_TRUE(flds->createView("a") == AXOM_NULLPTR);

  // attempt to create a group named the same as an existing view
  EXPECT_TRUE(flds->createGroup("a") == AXOM_NULLPTR);

  // attempt to create a view named the same as an existing group
  EXPECT_TRUE(ds->getRoot()->createView("fields") == AXOM_NULLPTR);

  ds->getRoot()->createGroup("here//is/path");
  ds->getRoot()->createGroup("éch≈o/Ωd");
  ds->getRoot()->createGroup("../group/..");

  IndexType idx = ds->getRoot()->getFirstValidGroupIndex();
  while (idx != InvalidIndex)
  {
    std::cout << ds->getRoot()->getGroup(idx)->getName() << std::endl;
    idx = ds->getRoot()->getNextValidGroupIndex(idx);
  }

  delete ds;
}
//------------------------------------------------------------------------------

TEST(sidre_group,view_copy_move)
{
  DataStore * ds = new DataStore();
  Group * flds = ds->getRoot()->createGroup("fields");
  int * buffdata;
  int extdata[10];

  DataView * views[6];
  std::string names[6];

  // Create view in different states
  views[0] = flds->createView("empty0");
  views[1] = flds->createView("empty1", INT_ID, 10);
  views[2] = flds->createViewAndAllocate("buffer", INT_ID, 10);
  views[3] = flds->createView("external", INT_ID, 10)->setExternalDataPtr(
    extdata);
  views[4] = flds->createViewScalar("scalar", 25);
  views[5] = flds->createViewString("string", "I am string");

  buffdata = flds->getView("buffer")->getData();
  for (int i=0 ; i < 10 ; ++i)
  {
    extdata[i] = i;
    buffdata[i] = i + 100;
  }

  for (int i = 0 ; i < 6 ; ++i)
  {
    names[i] = views[i]->getName();
    EXPECT_TRUE(flds->hasView(names[i]));
  }

  // test moving a view from flds to sub1
  Group * sub1 = flds->createGroup("sub1");

  // flds->print();

  for (int i = 0 ; i < 6 ; ++i)
  {
    sub1->moveView(views[i]);
    EXPECT_FALSE(flds->hasView(names[i]));
    EXPECT_TRUE(sub1->hasView(names[i]));

    // moving to same group is a no-op
    sub1->moveView(views[i]);
    EXPECT_TRUE(sub1->hasView(names[i]));
  }

  // flds->print();

  Group * sub2 = flds->createGroup("sub2");

  for (int i = 0 ; i < 6 ; ++i)
  {
    sub2->copyView(views[i]);
    EXPECT_TRUE(sub1->hasView(names[i]));
    EXPECT_TRUE(sub2->hasView(names[i]));
  }

  // Check copies
  DataView * view1 = sub1->getView("empty0");
  DataView * view2 = sub2->getView("empty0");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isEmpty());
  EXPECT_FALSE(view2->isDescribed());
  EXPECT_FALSE(view2->isAllocated());
  EXPECT_FALSE(view2->isApplied());

  view1 = sub1->getView("empty1");
  view2 = sub2->getView("empty1");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isEmpty());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_FALSE(view2->isAllocated());
  EXPECT_FALSE(view2->isApplied());

  view1 = sub1->getView("buffer");
  view2 = sub2->getView("buffer");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->hasBuffer());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_TRUE(view1->getBuffer() == view2->getBuffer());
  EXPECT_EQ(2, view1->getBuffer()->getNumViews());

  view1 = sub1->getView("external");
  view2 = sub2->getView("external");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isExternal());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_EQ(view1->getVoidPtr(), view2->getVoidPtr());

  view1 = sub1->getView("scalar");
  view2 = sub2->getView("scalar");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isScalar());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  EXPECT_EQ( view1->getData<int>(), view2->getData<int>());

  view1 = sub1->getView("string");
  view2 = sub2->getView("string");
  EXPECT_NE(view1, view2);
  EXPECT_TRUE(view2->isString());
  EXPECT_TRUE(view2->isDescribed());
  EXPECT_TRUE(view2->isAllocated());
  EXPECT_TRUE(view2->isApplied());
  const char * svalue = view1->getString();
  EXPECT_TRUE(strcmp("I am string", svalue) == 0);

  // flds->print();

  delete ds;
}

//------------------------------------------------------------------------------

TEST(sidre_group,groups_move_copy)
{
  DataStore * ds = new DataStore();
  Group * flds = ds->getRoot()->createGroup("fields");

  Group * ga = flds->createGroup("a");
  Group * gb = flds->createGroup("b");
  Group * gc = flds->createGroup("c");

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

//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_view_and_buffer2)
{
  DataStore * const ds = new DataStore();
  Group * const grp = ds->getRoot()->createGroup("grp");

  std::string viewName1("viewBuffer1");
  std::string viewName2("viewBuffer2");

  DataView * view1 = grp->createViewAndAllocate(viewName1, INT_ID, 1);
  DataView * view2 = grp->createViewAndAllocate(viewName2, INT_ID, 1);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndData(viewName1);

  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  Buffer const * const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE( buffer1 == AXOM_NULLPTR );

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
  Group * const grp = ds->getRoot()->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  DataView * const view1 = grp->createViewAndAllocate(viewName1,
                                                      DataType::c_int(10));

  EXPECT_TRUE(grp->hasChildView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  int * v1_vals = view1->getData();

  for(int i=0 ; i<10 ; i++)
  {
    v1_vals[i] = i;
  }

  EXPECT_EQ(view1->getNumElements(), 10u);
  EXPECT_EQ(view1->getTotalBytes(),
            static_cast<axom::sidre::SidreLength>(10 * sizeof(int)));

  grp->destroyViewAndData(viewName1);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,create_view_of_buffer_with_schema)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
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

  Buffer * base_buff = base->getBuffer();

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
TEST(sidre_group,save_restore_empty_datastore)
{
  const std::string file_path_base("sidre_empty_datastore_");
  DataStore * ds1 = new DataStore();

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    ds1->getRoot()->save(file_path, protocols[i]);
  }

  delete ds1;

  // Only restore conduit_hdf5
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    EXPECT_TRUE(ds2->getNumBuffers() == 0 );
    EXPECT_TRUE(root2->getNumGroups() == 0 );
    EXPECT_TRUE(root2->getNumViews() == 0 );

    delete ds2;
  }
}


//------------------------------------------------------------------------------
// make sure the hdf5 methods are consistent with the path based methods
//------------------------------------------------------------------------------
TEST(sidre_group,save_load_via_hdf5_ids)
{

  DataStore ds_save;
  // populate the datastore
  Group * root = ds_save.getRoot();
  root->createViewScalar<int>("i0", 1);
  root->createViewAndAllocate("vals", INT_ID, 5);
  // set values for the "vals" array
  int * vals_ptr =  root->getView("vals")->getData();
  for (int i = 0 ; i < 5 ; ++i)
  {
    vals_ptr[i] = i;
  }

  // save using the sidre_hdf5 protocol
  root->save("out_save_load_via_hdf5_ids.sidre_hdf5", "sidre_hdf5");

  // load via path based
  DataStore ds_load_generic;
  ds_load_generic.getRoot()->load("out_save_load_via_hdf5_ids.sidre_hdf5",
                                  "sidre_hdf5");

  // load via hdf5 id
  DataStore ds_load_hdf5;

  hid_t h5_id = H5Fopen("out_save_load_via_hdf5_ids.sidre_hdf5",
                        H5F_ACC_RDWR,
                        H5P_DEFAULT);
  EXPECT_TRUE(h5_id >= 0);

  // this implies protocol == "sidre_hdf5"
  ds_load_hdf5.getRoot()->load(h5_id);

  // ? Does isEquivalentTo check values?
  // check path based with source
  EXPECT_TRUE( ds_load_generic.getRoot()->isEquivalentTo(ds_save.getRoot()) );

  // check hdf5 based with source
  EXPECT_TRUE( ds_load_hdf5.getRoot()->isEquivalentTo(ds_save.getRoot()) );

  // check path based vs hdf5 based
  EXPECT_TRUE( ds_load_generic.getRoot()->isEquivalentTo(ds_load_hdf5.getRoot()) );

  // close hdf5 handle
  EXPECT_TRUE(H5Fclose(h5_id) >=0);
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_api)
{
  const std::string file_path_base("sidre_save_subtree_");
  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);

  // These should be produce identical files.

  // No group provided, defaults to root group
  root1->save("sidre_save_fulltree_conduit", "json");

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

#if 0
  DataStore * ds2 = new DataStore();
  Group * root2 = ds2->getRoot();
  root2->load("sidre_save_fulltree_conduit", "conduit");
  EXPECT_TRUE( ds2->getRoot()->isEquivalentTo(root1) );
  delete ds2;

  DataStore * ds3 = new DataStore();
  Group * root3 = ds3->getRoot();
  root3->load("sidre_save_subtree_conduit", "conduit", ds3->getRoot() );
  EXPECT_TRUE( ds3->getRoot()->isEquivalentTo(root1) );
  delete ds3;
#endif

  DataStore * ds4 = new DataStore();
  Group * root4 = ds4->getRoot();
  root4->load("sidre_save_subtree_sidre_hdf5", "sidre_hdf5");
  EXPECT_TRUE( ds4->getRoot()->isEquivalentTo(root1) );
  delete ds4;

  delete ds1;

  // Why don't these pass??? Need to ask Noah about this...
  // Trying to make sure sub trees are same here.
#if 0
  DataStore * ds_new = new DataStore();
  Group * tree1 = ds_new->getRoot()->createGroup("api1");
  Group * tree2 = ds_new->getRoot()->createGroup("api2");

  api1->load("sidre_save_subtree", "conduit");
  api2->load("sidre_save_subtree", "conduit");

  EXPECT_TRUE( api1->isEquivalentTo( api2) );
#endif

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_scalars_and_strings)
{
  const std::string file_path_base("sidre_save_scalars_and_strings_");
  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);
  root1->createViewScalar<float>("f0", 1.0);
  root1->createViewScalar<double>("d0", 10.0);
  root1->createViewString("s0", "I am a string");

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //	  continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }


  // Only restore conduit_hdf
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    EXPECT_TRUE( root1->isEquivalentTo( root2 ));

    int i0 = root2->getView("i0")->getScalar();
    float f0 = root2->getView("f0")->getScalar();
    double d0 = root2->getView("d0")->getScalar();
    const char * s0 = root2->getView("s0")->getString();

    EXPECT_EQ( 1, i0);
    EXPECT_EQ( 1.0, f0);
    EXPECT_EQ( 10.0, d0);
    EXPECT_EQ( std::string(s0), "I am a string");

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group,rename_group)
{
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();
  Group * child1 = root->createGroup("g_a");
  Group * child2 = root->createGroup("g_b");
  Group * child3 = root->createGroup("g_c");

  bool success = child1->rename("g_r");
  EXPECT_TRUE( success );
  EXPECT_TRUE( child1->getName() == "g_r" );
  EXPECT_TRUE( root->hasGroup("g_r") );
  EXPECT_FALSE( root->hasGroup("g_a") );

  success = child2->rename("fields/g_s");
  EXPECT_FALSE( success );
  EXPECT_TRUE( child2->getName() == "g_b" );

  success = child3->rename("g_b");
  EXPECT_FALSE( success );
  EXPECT_TRUE( child3->getName() == "g_c" );

}



//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_name_change)
{
  const std::string file_path_base("sidre_save_name_change_");
  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();
  Group * child1 = root1->createGroup("child1");

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  bool success = child1->getView("s0")->rename("s0_renamed");

  EXPECT_TRUE( success );
  EXPECT_FALSE( child1->hasView("s0") );
  EXPECT_TRUE( child1->hasView("s0_renamed") );

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //	  continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    child1->save(file_path, protocols[i]);
  }


  // Only restore conduit_hdf
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();
    Group * child2 = root2->createGroup("child2");

    EXPECT_EQ( child2->getName(), "child2" );

    child2->load(file_path, protocols[i]);

    EXPECT_EQ( child2->getName(), "child1" );

    EXPECT_TRUE( root1->isEquivalentTo( root2 ) );

    int i0 = child2->getView("i0")->getScalar();
    const char * s0 = child2->getView("s0_renamed")->getString();

    EXPECT_EQ( 1, i0 );
    EXPECT_EQ( std::string(s0), "I am a string" );

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_external_data)
{
  const std::string file_path_base("sidre_save_external_");

  int nfoo = 10;
  int foo1[nfoo], foo2[nfoo], * foo3, foo4[nfoo];
  int int2d1[nfoo*2], int2d2[nfoo*2];
  SidreLength shape[] = { nfoo, 2 };

  for (int i = 0 ; i < nfoo ; ++i)
  {
    foo1[i] = i;
    foo2[i] = 0;
    foo4[i] = i;
  }
  for (int i = 0 ; i < 2*nfoo ; ++i)
  {
    int2d1[i] = i;
    int2d2[i] = 0;
  }
  foo3 = NULL;

  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();

  root1->createView("external_array", INT_ID, nfoo, foo1);
  root1->createView("empty_array", INT_ID, nfoo, foo3);
  // XXX this falls into createView(name, type, ndims, shape)
  // root1->createView("empty_array", INT_ID, nfoo, NULL);
  root1->createView("external_undescribed")->setExternalDataPtr(foo4);
  root1->createView("int2d", INT_ID, 2, shape, int2d1);

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  delete ds1;

  // Now load back in.
  // Only restore conduit protocol_hdf5
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    SidreLength extents[7];
    int rank;

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    // load has set the type and size of the view.
    // Now set the external address before calling loadExternal.
    DataView * view1 = root2->getView("external_array");
    EXPECT_TRUE(view1->isExternal());
    EXPECT_TRUE(view1->isDescribed());
    EXPECT_EQ(view1->getTypeID(), INT_ID);
    EXPECT_EQ(view1->getNumElements(), nfoo);
    EXPECT_TRUE(view1->getVoidPtr() == AXOM_NULLPTR);
    view1->setExternalDataPtr(foo2);

    DataView * view2 = root2->getView("empty_array");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_TRUE(view2->getVoidPtr() == AXOM_NULLPTR);
    view2->setExternalDataPtr(foo3);

    DataView * view3 = root2->getView("external_undescribed");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_FALSE(view3->isDescribed());
    EXPECT_TRUE(view3->getVoidPtr() == AXOM_NULLPTR);
    // Set "external_array" and "external_undescribed" to the same external array
    // since it was created that way.  However, "external_undescribed" was not
    // written to the dump since it is undescribed.
    view3->setExternalDataPtr(foo2);

    DataView * view4 = root2->getView("int2d");
    EXPECT_FALSE(view4->isEmpty());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_TRUE(view4->getVoidPtr() == AXOM_NULLPTR);
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), nfoo*2);
    EXPECT_EQ(view4->getNumDimensions(), 2);
    rank = view4->getShape(7, extents);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(extents[0] == nfoo && extents[1] == 2);
    view4->setExternalDataPtr(int2d2);

    // Read external data into views
    root2->loadExternalData(file_path);

    // Make sure addresses have not changed
    EXPECT_TRUE(view1->getVoidPtr() == static_cast<void *>(foo2));
    EXPECT_TRUE(view2->getVoidPtr() == static_cast<void *>(foo3));    // AXOM_NULLPTR
    EXPECT_TRUE(view3->getVoidPtr() == static_cast<void *>(foo2));
    EXPECT_TRUE(view4->getVoidPtr() == static_cast<void *>(int2d2));

    for (int j = 0 ; i < nfoo ; ++i)
    {
      EXPECT_TRUE( foo1[j] == foo2[j] );
    }
    for (int j = 0 ; i < 2*nfoo ; ++i)
    {
      EXPECT_TRUE( int2d1[j] == int2d2[j] );
    }

    delete ds2;
  }
}

//------------------------------------------------------------------------------

// Check the association between views and buffers to make sure it is what we expect.
// This checks more than isEquivalentTo.

static void save_restore_buffer_association(const std::string & msg,
                                            DataStore * ds)
{
  const SidreLength len = 10;

  SCOPED_TRACE(msg);

  // Make sure all buffers were created
  if (ds->getNumBuffers() != 4u)
  {
    EXPECT_EQ(ds->getNumBuffers(), 4u);
    return;
  }

  Group * root = ds->getRoot();

  // Get all views and their buffers
  DataView * view1 = root->getView("undescribed_attached_buffer");
  ASSERT_TRUE(view1->hasBuffer());
  Buffer * buff1a = view1->getBuffer();
  ASSERT_FALSE(buff1a->isDescribed());
  ASSERT_FALSE(buff1a->isAllocated());

  DataView * view2 = root->getView("unallocated_attached_buffer");
  ASSERT_TRUE(view2->hasBuffer());
  Buffer * buff2a = view2->getBuffer();
  ASSERT_TRUE(buff2a->isDescribed());
  ASSERT_FALSE(buff2a->isAllocated());

  DataView * view3 = root->getView("undescribed_view_described_buffer");
  ASSERT_TRUE(view3->hasBuffer());
  Buffer * buff3a = view3->getBuffer();
  ASSERT_TRUE(buff3a->isDescribed());
  ASSERT_TRUE(buff3a->isAllocated());

  DataView * view4 = root->getView("describe_view_described_buffer");
  ASSERT_TRUE(view4->hasBuffer());
  Buffer * buff3b = view4->getBuffer();

  DataView * view5 = root->getView("even");
  ASSERT_TRUE(view5->hasBuffer());
  Buffer * buff3c = view5->getBuffer();

  DataView * view6 = root->getView("odd");
  ASSERT_TRUE(view6->hasBuffer());
  Buffer * buff3d = view6->getBuffer();

  DataView * view7 = root->getView("empty_described");
  ASSERT_FALSE(view7->hasBuffer());
  ASSERT_TRUE(view7->isEmpty());

  DataView * view8 = root->getView("allocated");
  ASSERT_TRUE(view8->hasBuffer());
  Buffer * buff4a = view8->getBuffer();
  ASSERT_TRUE(buff4a->isDescribed());
  ASSERT_TRUE(buff4a->isAllocated());

  // Check shared buffer  3a == 3b == 3c == 3d
  ASSERT_EQ(buff3a, buff3b);
  ASSERT_EQ(buff3b, buff3c);
  ASSERT_EQ(buff3c, buff3d);

  // Check unique buffers
  ASSERT_NE(buff1a, buff2a);
  ASSERT_NE(buff1a, buff3a);
  ASSERT_NE(buff1a, buff4a);

  ASSERT_NE(buff2a, buff3a);
  ASSERT_NE(buff2a, buff4a);

  ASSERT_NE(buff3a, buff4a);

  // Check contents of buffers
  ASSERT_EQ(buff3a->getNumElements(), len);
  int * idata = buff3a->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 100);
  }

  ASSERT_EQ(buff4a->getNumElements(), len);
  idata = buff4a->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 200);
  }

}

TEST(sidre_group,save_restore_buffer)
{
  const std::string file_path_base("sidre_save_buffer_");
  const SidreLength len = 10;

  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();
  Buffer * buff1 = ds1->createBuffer();
  Buffer * buff2 = ds1->createBuffer(INT_ID, len);
  Buffer * buff3 = ds1->createBuffer(INT_ID, len)->allocate();

  int * idata = buff3->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii + 100;
  }

  root1->createView("undescribed_attached_buffer", buff1);
  root1->createView("unallocated_attached_buffer", buff2);

  // These views share a buffer
  root1->createView("undescribed_view_described_buffer", buff3);
  root1->createView("describe_view_described_buffer", INT_ID, len, buff3);
  root1->createView("even", buff3)->apply(INT_ID, 5, 0, 2);
  root1->createView("odd", buff3)->apply(INT_ID, 5, 1, 2);

  root1->createView("empty_described", INT_ID, len);
  DataView * view = root1->createViewAndAllocate("allocated", INT_ID, len);

  idata = view->getData();
  for (int ii = 0 ; ii < len ; ++ii)
  {
    idata[ii] = ii + 200;
  }

  save_restore_buffer_association("original datastore", ds1);

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  // Now load back in.
  // Only restore conduit protocol_hdf5
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    bool isequivalent = root1->isEquivalentTo(root2);
    EXPECT_TRUE( isequivalent );
    if (isequivalent)
    {
      save_restore_buffer_association("loaded datastore", ds2);
    }

    delete ds2;
  }

  delete ds1;

}

//------------------------------------------------------------------------------

TEST(sidre_group,save_restore_other)
{
  const std::string file_path_base("sidre_save_other_");
  const int ndata = 10;
  SidreLength shape1[] = {ndata, 2};
  DataStore * ds1 = new DataStore();
  Group * root1 = ds1->getRoot();

  root1->createView("empty_view");
  root1->createView("empty_described", INT_ID, ndata);
  root1->createView("empty_shape", INT_ID, 2, shape1);

  root1->createViewAndAllocate("buffer_shape", INT_ID, 2, shape1);

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  delete ds1;

  // Now load back in.
  // Only restore conduit protocol_hdf5
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    SidreLength shape2[7];
    int rank;

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    DataView * view1 = root2->getView("empty_view");
    EXPECT_TRUE(view1->isEmpty());
    EXPECT_FALSE(view1->isDescribed());

    DataView * view2 = root2->getView("empty_described");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_EQ(view2->getNumElements(), ndata);

    DataView * view3 = root2->getView("empty_shape");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_TRUE(view3->isDescribed());
    EXPECT_EQ(view3->getTypeID(), INT_ID);
    EXPECT_EQ(view3->getNumElements(), ndata*2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view3->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    DataView * view4 = root2->getView("buffer_shape");
    EXPECT_TRUE(view4->hasBuffer());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), ndata*2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view4->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    delete ds2;
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_complex)
{
  const std::string file_path_base("sidre_mixed_types_");
  DataStore * ds1 = new DataStore();
  Group * flds = ds1->getRoot()->createGroup("fields");

  Group * ga = flds->createGroup("a");
  Group * gb = flds->createGroup("b");
  Group * gc = flds->createGroup("c");
  int ndata = 10;

  ga->createViewScalar<int>("i0", 100.0);
  ga->createViewScalar<double>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", INT_ID, ndata);
  int * data_ptr = gc->getView("int10")->getArray();
  for (int i = 0 ; i < ndata ; ++i)
  {
    data_ptr[i] = i;
  }

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    ds1->getRoot()->save(file_path, protocols[i]);
  }

  // Only restore conduit_hdf5 protocol
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();

    ds2->getRoot()->load(file_path, protocols[i]);

    EXPECT_TRUE( ds1->getRoot()->isEquivalentTo(ds2->getRoot()) );

    flds = ds2->getRoot()->getGroup("fields");

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getData<int>(),100.0);
    EXPECT_NEAR(flds->getGroup("a")->getView(
                  "d0")->getData<double>(),3000.0, 1e-12);

    int * new_data_ptr = flds->getGroup("c")->getView("int10")->getArray();
    for (int i = 0 ; i < ndata ; ++i)
    {
      EXPECT_TRUE( new_data_ptr[i] == i);
    }

    const char * char_ptr = flds->getView("b/s0")->getString();
    EXPECT_TRUE( std::string(char_ptr) == "foo" );

    //ds2->print();

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
// isEquivalentTo()
//------------------------------------------------------------------------------
TEST(sidre_group,is_equivalent_to)
{
  DataStore * ds = new DataStore();

  //These are the parents for two separate subtrees of the root group.
  //Everything below them will be created identically.
  Group * parent1 = ds->getRoot()->createGroup("parent1");
  Group * parent2 = ds->getRoot()->createGroup("parent2");

  //The flds1 and flds2 groups will be compared for equivalence
  Group * flds1 = parent1->createGroup("fields");
  Group * flds2 = parent2->createGroup("fields");

  Group * ga1 = flds1->createGroup("a");
  Group * gb1 = flds1->createGroup("b");
  Group * gc1 = flds1->createGroup("c");

  Group * gc2 = flds2->createGroup("c");    // Note: flds2 groups added in different order
  Group * gb2 = flds2->createGroup("b");
  Group * ga2 = flds2->createGroup("a");

  ga1->createViewScalar("i0", 1 );
  gb1->createViewScalar("f0", 100.0f );
  gc1->createViewScalar("d0", 3000.00);
  gc1->createViewScalar("d1", 6000.00);
  gc1->createViewScalar("d2", 9000.00);

  ga2->createViewScalar("i0", 1);
  gb2->createViewScalar("f0", 100.0f);
  gc2->createViewScalar("d2", 9000.00);         // Note: views of gc2 added in different order
  gc2->createViewScalar("d1", 6000.00);
  gc2->createViewScalar("d0", 3000.00);


  // Groups were created identically, so should be equivalent.
  EXPECT_TRUE(flds1->isEquivalentTo(flds2));

  // Add something extra to flds2, making them not equivalent.
  gc2->createViewAndAllocate("extra", DataType::c_double());

  EXPECT_FALSE(flds1->isEquivalentTo(flds2));

  delete ds;

}



//------------------------------------------------------------------------------
TEST(sidre_group,save_load_all_protocols)
{
  const std::string file_path_base("sidre_save_load_all_protocols.");
  DataStore ds;

  Group * flds = ds.getRoot()->createGroup("fields");

  Group * ga = flds->createGroup("a");
  Group * gb = flds->createGroup("b");
  Group * gc = flds->createGroup("c");
  int ndata = 10;

  // prep a tree that can exactly restored by all
  // i/o protocols.
  // Specially, use int64 and float64 b/c the
  // json i/o case uses those types for parsed integers
  // and floating point numbers.

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64 * data_ptr = gc->getView("int10")->getArray();
  for (int i = 0 ; i < ndata ; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  // show the source tree
  SLIC_INFO("Source tree");
  ds.print();

  //
  // test all protocols
  //
  std::vector<std::string> protocols;
  protocols.push_back("sidre_hdf5");
  protocols.push_back("sidre_conduit_json");
  protocols.push_back("sidre_json");

  protocols.push_back("conduit_hdf5");
  protocols.push_back("conduit_bin");
  protocols.push_back("conduit_json");
  protocols.push_back("json");

  for (size_t i = 0 ; i < protocols.size() ; ++i)
  {
    SLIC_INFO("Testing protocol: " << protocols[i]);
    const std::string file_path = file_path_base + protocols[i];
    // save using current protocol
    ds.getRoot()->save(file_path, protocols[i]);

    DataStore ds_load;
    ds_load.getRoot()->load(file_path, protocols[i]);

    SLIC_INFO("Tree from protocol: " <<  protocols[i]);
    // show the result
    ds_load.print();

    Group * ds_load_root = ds_load.getRoot();
    // check that the sidre hierarchy is equiv
    EXPECT_TRUE( ds.getRoot()->isEquivalentTo(ds_load_root));

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView(
                "fields/a/i0")->getData<conduit::int64>(),100);
    EXPECT_NEAR(ds_load_root->getView(
                  "fields/a/d0")->getData<conduit::float64>(),3000.00,1e-12);
    EXPECT_EQ(ds_load_root->getView("fields/b/s0")->getString(),
              std::string("foo"));

    conduit::int64 * load_data_ptr =
      ds_load_root->getView("fields/c/int10")->getData();
    for(int j=0 ; j< ndata ; j++)
    {
      EXPECT_EQ(data_ptr[j],load_data_ptr[j]);
    }

  }
}
