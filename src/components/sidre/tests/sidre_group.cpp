/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "axom/config.hpp"  // for AXOM_USE_HDF5
#include "sidre/sidre.hpp"

#include "gtest/gtest.h"

#include <cstring>
#include <vector>

using axom::sidre::SidreLength;
using axom::sidre::TypeID;
using axom::sidre::Buffer;
using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::View;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::nameIsValid;
using axom::sidre::indexIsValid;
using axom::sidre::DataType;
using axom::sidre::INT_ID;
using axom::sidre::FLOAT64_ID;

namespace
{
// Test protocols
#ifdef AXOM_USE_HDF5
int nprotocols = 3;
std::string const protocols[] = { "sidre_json", "sidre_hdf5", "json" };
#else
int nprotocols = 2;
std::string const protocols[] = { "sidre_json", "json" };
#endif

// Function to return a vector of available sidre protocols
std::vector<std::string> getAvailableSidreProtocols()
{
  std::vector<std::string> protocols;

#ifdef AXOM_USE_HDF5
  protocols.push_back("sidre_hdf5");
#endif
  protocols.push_back("sidre_json");
  protocols.push_back("sidre_conduit_json");

#ifdef AXOM_USE_HDF5
  protocols.push_back("conduit_hdf5");
#endif
  protocols.push_back("conduit_bin");
  protocols.push_back("conduit_json");
  protocols.push_back("json");

  return protocols;
}



} // end anonymous namespace

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("test");

  EXPECT_TRUE(root->isRoot());
  EXPECT_FALSE(group->isRoot());

  EXPECT_TRUE(group->getName() == std::string("test") );

  Group* group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == AXOM_NULLPTR);
}

//------------------------------------------------------------------------------
// getPath(), getPathName()
//------------------------------------------------------------------------------
TEST(sidre_group,get_path_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  EXPECT_EQ(root->getParent(), root);
  EXPECT_EQ(root->getName(), "");
  Group* group = root->createGroup("test/a/b/c");
  Group* grp2 = root->getGroup("test/a");
  Group* grp3 = root->getGroup("test");

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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test full path access when building incrementally
  Group* group =
    root->createGroup("test1")->createGroup("test2")->createGroup("test3");
  Group* group2 = root->getGroup("test1/test2/test3");

  EXPECT_TRUE(AXOM_NULLPTR != group2);
  EXPECT_EQ(group, group2);

  // Test incremental access when building full path
  Group* groupP = root->createGroup("testA/testB/testC");
  Group* groupP2 =
    root->getGroup("testA")->getGroup("testB")->getGroup("testC");

  EXPECT_TRUE(AXOM_NULLPTR != groupP2);
  EXPECT_EQ(groupP, groupP2);
  // test non-const getGroup() with path
  Group* groupPParent = root->getGroup("testA/testB");
  EXPECT_EQ(groupP->getParent(), groupPParent);
  EXPECT_EQ(groupP->getParent()->getName(), "testB");


  // Now verify that code will not create missing groups.

  root->createGroup("testa")->createGroup("testb")->createGroup("testc");
  Group* group_bada = root->getGroup("BAD/testb/testc");
  Group* group_badb = root->getGroup("testa/BAD/testc");
  Group* group_badc = root->getGroup("testa/testb/BAD");

  EXPECT_EQ(group_bada, static_cast<void*>(AXOM_NULLPTR) );
  EXPECT_EQ(group_badb, static_cast<void*>(AXOM_NULLPTR) );
  EXPECT_EQ(group_badc, static_cast<void*>(AXOM_NULLPTR) );

  // Test hasGroup with paths.

  EXPECT_FALSE(root->hasGroup("BAD/testb/testc"));
  EXPECT_FALSE(root->hasGroup("testa/BAD/testc"));
  EXPECT_FALSE(root->hasGroup("testa/testb/BAD"));

  EXPECT_TRUE(root->hasGroup("test1"));
  EXPECT_TRUE(root->hasGroup("test1/test2"));
  EXPECT_TRUE(root->hasGroup("test1/test2/test3"));
  Group* group_testa = root->getGroup("testa");
  EXPECT_TRUE(group_testa->hasGroup("testb"));
  EXPECT_TRUE(group_testa->hasGroup("testb/testc"));
  EXPECT_FALSE(group_testa->hasGroup("testb/BAD"));
  EXPECT_FALSE(group_testa->hasGroup("testb/testc/BAD"));

  EXPECT_EQ(root->getNumGroups(), 3u);
  EXPECT_TRUE(root->hasGroup(0));
  EXPECT_TRUE(root->hasGroup(1));
  EXPECT_TRUE(root->hasGroup(2));
  EXPECT_FALSE(root->hasGroup(3));
  EXPECT_FALSE(root->hasGroup(InvalidIndex));

  unsigned int testbnumgroups = group_testa->getGroup("testb")->getNumGroups();
  Group* group_cdup = group_testa->createGroup("testb/testc");

  EXPECT_EQ(group_cdup, static_cast<void*>(AXOM_NULLPTR));
  EXPECT_EQ(group_testa->getGroup("testb")->getNumGroups(), testbnumgroups);

  delete ds;

}

//------------------------------------------------------------------------------
// createGroup(), destroyGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group,destroy_group_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test full path access when building incrementally
  Group* group = root->createGroup("test1/test2/test3");
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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");

  EXPECT_TRUE( child->getParent() == parent );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group,get_datastore)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  EXPECT_TRUE( group->getDataStore() == ds );

  DataStore const* const_ds = group->getDataStore();
  EXPECT_TRUE( const_ds == ds );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(sidre_group,get_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");
  EXPECT_TRUE( child->getParent() == parent );

  EXPECT_TRUE( parent->getGroup("child") == child );
  EXPECT_TRUE( parent->getGroup(0) == child );

  // check error condition
  EXPECT_TRUE( parent->getGroup("non-existant group") == AXOM_NULLPTR );

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view = parent->createView("view");

  EXPECT_TRUE( parent->getView("view") == view );
  EXPECT_TRUE( parent->getView(0) == view );

  // check error condition
  EXPECT_TRUE( parent->getView("non-existant view") == AXOM_NULLPTR );

  delete ds;
}

//------------------------------------------------------------------------------
// createView, hasView(), getView(), destroyView() with path strings
//------------------------------------------------------------------------------
TEST(sidre_group,view_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test with full path access when building incrementally
  View* view =
    root->createGroup("group1")->createGroup("group2")->createView("view1");
  View* view2 = root->getView("group1/group2/view1");

  EXPECT_TRUE(AXOM_NULLPTR != view2);
  EXPECT_EQ( view, view2 );


  // Test incremental access when building with full path
  View* viewP = root->createView("groupA/groupB/viewA");
  View* viewP2 =
    root->getGroup("groupA")->getGroup("groupB")->getView("viewA");

  EXPECT_TRUE(AXOM_NULLPTR != viewP2);
  EXPECT_EQ( viewP, viewP2 );

  // Now verify that bad paths just return null, and don't create missing groups
  View* v_bad1 = root->getView("BAD/groupB/viewA");
  View* v_bad2 = root->getView("groupA/BAD/viewA");
  View* v_bad3 = root->getView("groupA/groupB/BAD");

  EXPECT_EQ(v_bad1, static_cast<void*>(AXOM_NULLPTR));
  EXPECT_EQ(v_bad2, static_cast<void*>(AXOM_NULLPTR));
  EXPECT_EQ(v_bad3, static_cast<void*>(AXOM_NULLPTR));

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
            static_cast<void*>(AXOM_NULLPTR));
  EXPECT_FALSE(root->hasView("group1/group2/view1"));
  EXPECT_EQ(root->getView("group1/group2/view1"),
            static_cast<void*>(AXOM_NULLPTR));

  Group* groupA = root->getGroup("groupA");
  EXPECT_TRUE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), viewP);
  EXPECT_TRUE(root->hasView("groupA/groupB/viewA"));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), viewP);

  groupA->destroyView("groupB/viewA");

  EXPECT_EQ(exp_no_groups, groupA->getGroup("groupB")->getNumViews());
  EXPECT_FALSE(groupA->getGroup("groupB")->hasView("viewA"));
  EXPECT_EQ(groupA->getGroup("groupB")->getView("viewA"),
            static_cast<void*>(AXOM_NULLPTR));
  EXPECT_FALSE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), static_cast<void*>(AXOM_NULLPTR));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"),
            static_cast<void*>(AXOM_NULLPTR));

  delete ds;
}


//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group,get_view_names_and_indicies)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view1 = parent->createView("view1");
  View* view2 = parent->createView("view2");

  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getViewIndex("view1");
  IndexType idx2 = parent->getViewIndex("view2");
  EXPECT_EQ(idx1, view1->getIndex());
  EXPECT_EQ(idx2, view2->getIndex());

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
// Verify getFirstValidGroupIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_group_index)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* group1 = parent->createGroup("group1");
  Group* group2 = parent->createGroup("group2");
  EXPECT_EQ(parent->getNumGroups(), 2u);

  IndexType idx1 = parent->getFirstValidGroupIndex();
  IndexType idx2 = parent->getNextValidGroupIndex(idx1);
  IndexType idx3 = parent->getNextValidGroupIndex(idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(InvalidIndex, idx3);

  Group* group1out = parent->getGroup(idx1);
  Group* group2out = parent->getGroup(idx2);
  EXPECT_EQ(group1, group1out);
  EXPECT_EQ(group2, group2out);

  // check error conditions
  Group* emptyGroup = root->createGroup("emptyGroup");
  IndexType badidx1 = emptyGroup->getFirstValidGroupIndex();
  IndexType badidx2 = emptyGroup->getNextValidGroupIndex(badidx1);

  EXPECT_TRUE(badidx1 == InvalidIndex);
  EXPECT_TRUE(badidx2 == InvalidIndex);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group,get_first_and_next_view_index)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view1 = parent->createView("view1");
  View* view2 = parent->createView("view2");
  EXPECT_EQ(parent->getNumViews(), 2u);

  IndexType idx1 = parent->getFirstValidViewIndex();
  IndexType idx2 = parent->getNextValidViewIndex(idx1);
  IndexType idx3 = parent->getNextValidViewIndex(idx2);
  EXPECT_EQ(0, idx1);
  EXPECT_EQ(1, idx2);
  EXPECT_EQ(InvalidIndex, idx3);

  View* view1out = parent->getView(idx1);
  View* view2out = parent->getView(idx2);
  EXPECT_EQ(view1, view1out);
  EXPECT_EQ(view2, view2out);

  // check error conditions
  Group* emptyGroup = root->createGroup("emptyGroup");
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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* group1 = parent->createGroup("group1");
  Group* group2 = parent->createGroup("group2");

  EXPECT_EQ(parent->getNumGroups(), 2u);

  IndexType idx1 = parent->getGroupIndex("group1");
  IndexType idx2 = parent->getGroupIndex("group2");
  EXPECT_EQ(idx1, group1->getIndex());
  EXPECT_EQ(idx2, group2->getIndex());

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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  View* view = group->createView("view");
  EXPECT_TRUE( group->getParent() == root );
  EXPECT_FALSE( view->hasBuffer() );
  EXPECT_TRUE( group->hasView("view") );
  IndexType iview = group->getViewIndex("view");
  EXPECT_EQ(0, iview);
  iview = view->getIndex();
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
  Buffer* buff = view->getBuffer();
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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("group");
  EXPECT_TRUE( group->getParent() == root );

  EXPECT_TRUE( root->hasGroup("group") );

  root->destroyGroup("group");
  EXPECT_FALSE( root->hasGroup("group") );

  // should be a no-op, not a failure
  root->destroyGroup("group");

  Group* group2 = root->createGroup("group2");
  // shut up compiler about unused variable
  (void)group2;
  root->destroyGroup( root->getFirstValidGroupIndex() );

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,group_name_collisions)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");
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
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");
  int* buffdata;
  int extdata[10];

  View* views[6];
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
  Group* sub1 = flds->createGroup("sub1");

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

  Group* sub2 = flds->createGroup("sub2");

  for (int i = 0 ; i < 6 ; ++i)
  {
    sub2->copyView(views[i]);
    EXPECT_TRUE(sub1->hasView(names[i]));
    EXPECT_TRUE(sub2->hasView(names[i]));
  }

  // Check copies
  View* view1 = sub1->getView("empty0");
  View* view2 = sub2->getView("empty0");
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
  const char* svalue = view1->getString();
  EXPECT_TRUE(strcmp("I am string", svalue) == 0);

  // flds->print();

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,groups_move_copy)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");

  const double f0value = 100.0;
  const double val = 101.0;

  Group* bschild = gb->createGroup("childOfB");

  ga->createView("i0")->setScalar(1);
  gb->createView("f0")->setScalar(f0value);
  gc->createView("d0")->setScalar(3000.0);
  bschild->createView("val")->setScalar(val);

  int buffercount = ds->getNumBuffers();

  // check that all sub groups exist
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("b"));
  EXPECT_TRUE(flds->hasGroup("c"));

  // move "b" to a child of "sub"
  EXPECT_EQ(1, gb->getIndex());
  EXPECT_EQ(flds, gb->getParent());
  Group* gsub = flds->createGroup("sub");
  Group* gb0 = gsub->moveGroup(gb);

  // gb0 is an alias to gb
  EXPECT_EQ(gb, gb0);
  EXPECT_EQ(0, gb->getIndex());
  EXPECT_EQ(gsub, gb->getParent());
  EXPECT_EQ(gb0->getNumGroups(), 1);
  EXPECT_EQ(gb0->getGroup("childOfB"), bschild);
  EXPECT_EQ(bschild->getNumGroups(), 0);
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  EXPECT_EQ(gb0->getNumViews(), 1);
  EXPECT_TRUE(gb0->hasChildView("f0"));
  if (gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  EXPECT_EQ(bschild->getNumViews(), 1);
  EXPECT_TRUE(bschild->hasChildView("val"));
  if (gb0->hasChildView("val"))
  {
    EXPECT_EQ((double)gb0->getView("val")->getScalar(), val);
  }

  EXPECT_EQ(flds->getNumGroups(), 3);
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->hasGroup("c"));

  EXPECT_EQ(flds->getGroup("sub")->getGroup("b"),gb);

  // verify that we can copy a group into an empty group
  Group* containCopy = ds->getRoot()->createGroup("containCopy");
  Group* theCopy = containCopy->copyGroup(flds);
  EXPECT_TRUE(theCopy->isEquivalentTo(flds));
  EXPECT_EQ(containCopy->getNumGroups(), 1);
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  // verify that we can copy a group, when there is no name clash
  Group* anotherCopy = ds->getRoot()->createGroup("anotherCopy");
  anotherCopy->createGroup("futureSiblingGroup");
  Group* theOtherCopy = anotherCopy->copyGroup(flds);
  EXPECT_EQ(anotherCopy->getNumGroups(), 2);
  EXPECT_TRUE(theOtherCopy->isEquivalentTo(flds));
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  // verify that we cannot copy a group when there is a name clash
  Group* otherB = containCopy->createGroup("b");
  otherB->createView("f1")->setScalar(42.0);
  otherB->createGroup("Q");
  Group* triedCopy = gsub->copyGroup(otherB);
  EXPECT_EQ(triedCopy, static_cast<void*>(AXOM_NULLPTR));
  EXPECT_EQ(gsub->getNumGroups(), 1);
  EXPECT_TRUE(gsub->hasChildGroup("b"));
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  EXPECT_EQ(gb0->getNumGroups(), 1);
  EXPECT_EQ(gb0->getGroup("childOfB"), bschild);
  EXPECT_EQ(gb0->getNumViews(), 1);
  EXPECT_TRUE(gb0->hasChildView("f0"));
  if (gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group,create_destroy_view_and_buffer2)
{
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string viewName1("viewBuffer1");
  std::string viewName2("viewBuffer2");

  View* view1 = grp->createViewAndAllocate(viewName1, INT_ID, 1);
  View* view2 = grp->createViewAndAllocate(viewName2, INT_ID, 1);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndData(viewName1);

  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  Buffer const* const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE( buffer1 == AXOM_NULLPTR );

  View const* const view3 = grp->createView("viewBuffer3");
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
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  View* const view1 = grp->createViewAndAllocate(viewName1,
                                                 DataType::c_int(10));

  EXPECT_TRUE(grp->hasChildView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  int* v1_vals = view1->getData();

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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  // use create + alloc convenience methods
  // this one is the DataType & method
  View* base =  root->createViewAndAllocate("base", DataType::c_int(10));
  int* base_vals = base->getData();
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

  Buffer* base_buff = base->getBuffer();

  // create two views into this buffer
  //
  // view for the first 5 values
  root->createView("sub_a", base_buff)->apply(DataType::c_int(5));

  int* sub_a_vals = root->getView("sub_a")->getData();

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
  DataStore* ds1 = new DataStore();

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

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    EXPECT_TRUE(ds2->getNumBuffers() == 0 );
    EXPECT_TRUE(root2->getNumGroups() == 0 );
    EXPECT_TRUE(root2->getNumViews() == 0 );

    delete ds2;
  }
}

#ifdef AXOM_USE_HDF5
//------------------------------------------------------------------------------
// make sure the hdf5 methods are consistent with the path based methods
//------------------------------------------------------------------------------
TEST(sidre_group,save_load_via_hdf5_ids)
{

  DataStore ds_save;
  // populate the datastore
  Group* root = ds_save.getRoot();
  root->createViewScalar<int>("i0", 1);
  root->createViewAndAllocate("vals", INT_ID, 5);
  // set values for the "vals" array
  int* vals_ptr =  root->getView("vals")->getData();
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
#endif  // AXOM_USE_HDF5

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_api)
{
  const std::string file_path_base("sidre_save_subtree_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);

  // These should be produce identical files.

  // No group provided, defaults to root group
  root1->save("sidre_save_fulltree_conduit", "json");

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

//These are commented out because createViewScalar<int> creates a scalar View
//with a 32-bit int, but when conduit reads a raw integer from json file, it
//stores it as a 64-bit int, so the isEquivalentTo test fails.

#if 0
  DataStore* ds2 = new DataStore();
  Group* root2 = ds2->getRoot();
  root2->load("sidre_save_fulltree_conduit", "json");
  EXPECT_TRUE( ds2->getRoot()->isEquivalentTo(root1) );
  delete ds2;

  DataStore* ds3 = new DataStore();
  Group* root3 = ds3->getRoot();
  root3->load("sidre_save_subtree_sidre_json", "sidre_json");
  EXPECT_TRUE( ds3->getRoot()->isEquivalentTo(root1) );
  delete ds3;
#endif

#ifdef AXOM_USE_HDF5
  DataStore* ds4 = new DataStore();
  Group* root4 = ds4->getRoot();
  root4->load("sidre_save_subtree_sidre_hdf5", "sidre_hdf5");
  EXPECT_TRUE( ds4->getRoot()->isEquivalentTo(root1) );
  delete ds4;
#endif

#if 0
  DataStore* ds5 = new DataStore();
  Group* root5 = ds5->getRoot();
  root5->load("sidre_save_subtree_json", "json");
  EXPECT_TRUE( ds5->getRoot()->isEquivalentTo(root1) );
  delete ds5;
#endif
  delete ds1;

  // Test loading of same subtree to different parts of a Group tree
  DataStore* ds_new = new DataStore();
  Group* tree1 = ds_new->getRoot()->createGroup("api1");
  Group* tree2 = ds_new->getRoot()->createGroup("api2");

  Group* load1 = tree1->createGroup("subtree");
  Group* load2 = tree2->createGroup("subtree");

  load1->load("sidre_save_subtree_sidre_json", "sidre_json");
  load2->load("sidre_save_subtree_sidre_json", "sidre_json");

  EXPECT_TRUE( load1->isEquivalentTo( load2) );

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_scalars_and_strings)
{
  const std::string file_path_base("sidre_save_scalars_and_strings_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);
  root1->createViewScalar<float>("f0", 1.0);
  root1->createViewScalar<double>("d0", 10.0);
  root1->createViewString("s0", "I am a string");

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    EXPECT_TRUE( root1->isEquivalentTo( root2 ));

    int i0 = root2->getView("i0")->getScalar();
    float f0 = root2->getView("f0")->getScalar();
    double d0 = root2->getView("d0")->getScalar();
    const char* s0 = root2->getView("s0")->getString();

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
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* child1 = root->createGroup("g_a");
  Group* child2 = root->createGroup("g_b");
  Group* child3 = root->createGroup("g_c");

  // rename should not change the index
  EXPECT_EQ(0, child1->getIndex());
  bool success = child1->rename("g_r");
  EXPECT_TRUE( success );
  EXPECT_EQ( "g_r", child1->getName() );
  EXPECT_EQ(0, child1->getIndex());
  EXPECT_TRUE( root->hasGroup("g_r") );
  EXPECT_FALSE( root->hasGroup("g_a") );

  // try to rename to path
  success = child2->rename("fields/g_s");
  EXPECT_FALSE( success );
  EXPECT_EQ( "g_b", child2->getName() );

  // Try to rename to existing group name
  success = child3->rename("g_b");
  EXPECT_FALSE( success );
  EXPECT_EQ( "g_c", child3->getName() );

  // Rename root group
  EXPECT_EQ(InvalidIndex, root->getIndex());
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("", root->getName());
  root->rename("newroot");
  EXPECT_EQ(InvalidIndex, root->getIndex());
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("newroot", root->getName());

}

//------------------------------------------------------------------------------
TEST(sidre_group,save_restore_name_change)
{
  const std::string file_path_base("sidre_save_name_change_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Group* child1 = root1->createGroup("child1");

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  bool success = child1->getView("s0")->rename("s0_renamed");

  EXPECT_TRUE( success );
  EXPECT_FALSE( child1->hasView("s0") );
  EXPECT_TRUE( child1->hasView("s0_renamed") );

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    child1->save(file_path, protocols[i]);
  }


  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();
    Group* child2 = root2->createGroup("child2");

    EXPECT_EQ( child2->getName(), "child2" );

    child2->load(file_path, protocols[i]);

    EXPECT_EQ( child2->getName(), "child1" );

    EXPECT_TRUE( root1->isEquivalentTo( root2 ) );

    int i0 = child2->getView("i0")->getScalar();
    const char* s0 = child2->getView("s0_renamed")->getString();

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

  const int nfoo = 10;
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

  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

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
  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];
    SidreLength extents[7];
    int rank;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    // load has set the type and size of the view.
    // Now set the external address before calling loadExternal.
    View* view1 = root2->getView("external_array");
    EXPECT_TRUE(view1->isExternal());
    EXPECT_TRUE(view1->isDescribed());
    EXPECT_EQ(view1->getTypeID(), INT_ID);
    EXPECT_EQ(view1->getNumElements(), nfoo);
    EXPECT_TRUE(view1->getVoidPtr() == AXOM_NULLPTR);
    view1->setExternalDataPtr(foo2);

    View* view2 = root2->getView("empty_array");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_TRUE(view2->getVoidPtr() == AXOM_NULLPTR);
    view2->setExternalDataPtr(foo3);

    View* view3 = root2->getView("external_undescribed");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_FALSE(view3->isDescribed());
    EXPECT_TRUE(view3->getVoidPtr() == AXOM_NULLPTR);
    // Set "external_array" and "external_undescribed" to the same external
    // array
    // since it was created that way.  However, "external_undescribed" was not
    // written to the dump since it is undescribed.
    view3->setExternalDataPtr(foo2);

    View* view4 = root2->getView("int2d");
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
    EXPECT_TRUE(view1->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view2->getVoidPtr() == static_cast<void*>(foo3));     // AXOM_NULLPTR
    EXPECT_TRUE(view3->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view4->getVoidPtr() == static_cast<void*>(int2d2));

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

// Check the association between views and buffers to make sure it is what we
// expect.
// This checks more than isEquivalentTo.

static void save_restore_buffer_association(const std::string & msg,
                                            DataStore* ds)
{
  const SidreLength len = 10;

  SCOPED_TRACE(msg);

  // Make sure all buffers were created
  if (ds->getNumBuffers() != 4u)
  {
    EXPECT_EQ(ds->getNumBuffers(), 4u);
    return;
  }

  Group* root = ds->getRoot();

  // Get all views and their buffers
  View* view1 = root->getView("undescribed_attached_buffer");
  ASSERT_TRUE(view1->hasBuffer());
  Buffer* buff1a = view1->getBuffer();
  ASSERT_FALSE(buff1a->isDescribed());
  ASSERT_FALSE(buff1a->isAllocated());

  View* view2 = root->getView("unallocated_attached_buffer");
  ASSERT_TRUE(view2->hasBuffer());
  Buffer* buff2a = view2->getBuffer();
  ASSERT_TRUE(buff2a->isDescribed());
  ASSERT_FALSE(buff2a->isAllocated());

  View* view3 = root->getView("undescribed_view_described_buffer");
  ASSERT_TRUE(view3->hasBuffer());
  Buffer* buff3a = view3->getBuffer();
  ASSERT_TRUE(buff3a->isDescribed());
  ASSERT_TRUE(buff3a->isAllocated());

  View* view4 = root->getView("describe_view_described_buffer");
  ASSERT_TRUE(view4->hasBuffer());
  Buffer* buff3b = view4->getBuffer();

  View* view5 = root->getView("even");
  ASSERT_TRUE(view5->hasBuffer());
  Buffer* buff3c = view5->getBuffer();

  View* view6 = root->getView("odd");
  ASSERT_TRUE(view6->hasBuffer());
  Buffer* buff3d = view6->getBuffer();

  View* view7 = root->getView("empty_described");
  ASSERT_FALSE(view7->hasBuffer());
  ASSERT_TRUE(view7->isEmpty());

  View* view8 = root->getView("allocated");
  ASSERT_TRUE(view8->hasBuffer());
  Buffer* buff4a = view8->getBuffer();
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
  int* idata = buff3a->getData();
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

  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Buffer* buff1 = ds1->createBuffer();
  Buffer* buff2 = ds1->createBuffer(INT_ID, len);
  Buffer* buff3 = ds1->createBuffer(INT_ID, len)->allocate();

  int* idata = buff3->getData();
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
  View* view = root1->createViewAndAllocate("allocated", INT_ID, len);

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
  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

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
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

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
  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];
    SidreLength shape2[7];
    int rank;

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);

    View* view1 = root2->getView("empty_view");
    EXPECT_TRUE(view1->isEmpty());
    EXPECT_FALSE(view1->isDescribed());

    View* view2 = root2->getView("empty_described");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_EQ(view2->getNumElements(), ndata);

    View* view3 = root2->getView("empty_shape");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_TRUE(view3->isDescribed());
    EXPECT_EQ(view3->getTypeID(), INT_ID);
    EXPECT_EQ(view3->getNumElements(), ndata*2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view3->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    View* view4 = root2->getView("buffer_shape");
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
  DataStore* ds1 = new DataStore();
  Group* flds = ds1->getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
  int ndata = 10;

  ga->createViewScalar<int>("i0", 100.0);
  ga->createViewScalar<double>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", INT_ID, ndata);
  int* data_ptr = gc->getView("int10")->getArray();
  for (int i = 0 ; i < ndata ; ++i)
  {
    data_ptr[i] = i;
  }

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    ds1->getRoot()->save(file_path, protocols[i]);
  }

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();

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

    int* new_data_ptr = flds->getGroup("c")->getView("int10")->getArray();
    for (int i = 0 ; i < ndata ; ++i)
    {
      EXPECT_TRUE( new_data_ptr[i] == i);
    }

    const char* char_ptr = flds->getView("b/s0")->getString();
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
  DataStore* ds = new DataStore();

  //These are the parents for two separate subtrees of the root group.
  //Everything below them will be created identically.
  Group* parent1 = ds->getRoot()->createGroup("parent1");
  Group* parent2 = ds->getRoot()->createGroup("parent2");

  //The flds1 and flds2 groups will be compared for equivalence
  Group* flds1 = parent1->createGroup("fields");
  Group* flds2 = parent2->createGroup("fields");

  Group* ga1 = flds1->createGroup("a");
  Group* gb1 = flds1->createGroup("b");
  Group* gc1 = flds1->createGroup("c");

  Group* gc2 = flds2->createGroup("c");     // Note: flds2 groups added in
                                            // different order
  Group* gb2 = flds2->createGroup("b");
  Group* ga2 = flds2->createGroup("a");

  ga1->createViewScalar("i0", 1 );
  gb1->createViewScalar("f0", 100.0f );
  gc1->createViewScalar("d0", 3000.00);
  gc1->createViewScalar("d1", 6000.00);
  gc1->createViewScalar("d2", 9000.00);

  ga2->createViewScalar("i0", 1);
  gb2->createViewScalar("f0", 100.0f);
  gc2->createViewScalar("d2", 9000.00);         // Note: views of gc2 added in
                                                // different order
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

  Group* flds = ds.getRoot()->createGroup("fields");

  Group* ga = flds->createGroup("a");
  Group* gb = flds->createGroup("b");
  Group* gc = flds->createGroup("c");
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
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
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
  std::vector<std::string> protocols = getAvailableSidreProtocols();
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

    Group* ds_load_root = ds_load.getRoot();
    // check that the sidre hierarchy is equiv
    EXPECT_TRUE( ds.getRoot()->isEquivalentTo(ds_load_root));

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView(
                "fields/a/i0")->getData<conduit::int64>(),100);
    EXPECT_NEAR(ds_load_root->getView(
                  "fields/a/d0")->getData<conduit::float64>(),3000.00,1e-12);
    EXPECT_EQ(ds_load_root->getView("fields/b/s0")->getString(),
              std::string("foo"));

    conduit::int64* load_data_ptr =
      ds_load_root->getView("fields/c/int10")->getData();
    for(int j=0 ; j< ndata ; j++)
    {
      EXPECT_EQ(data_ptr[j],load_data_ptr[j]);
    }

  }
}


//------------------------------------------------------------------------------
TEST(sidre_group,save_load_preserve_contents)
{
  const std::string file_path_tree0("sidre_save_load_preserve_contents.tree0.");
  const std::string file_path_tree1("sidre_save_load_preserve_contents.tree1.");
  DataStore ds;

  Group* tree0 = ds.getRoot()->createGroup("tree0");

  Group* ga = tree0->createGroup("a");
  Group* gb = tree0->createGroup("b");
  Group* gc = tree0->createGroup("c");
  int ndata = 10;

  // Prepare a tree that can be exactly restored by all I/O protocols.
  // Specifically, use int64 and float64 because the json I/O case
  // uses those types for parsed integers and floating point numbers.

  ga->createViewScalar<conduit::int64>("i0", 100);
  ga->createViewScalar<conduit::float64>("d0", 3000.00);
  gb->createViewString("s0", "foo");

  gc->createViewAndAllocate("int10", DataType::int64(ndata));
  conduit::int64* data_ptr = gc->getView("int10")->getArray();
  for (int i = 0 ; i < ndata ; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  std::vector<std::string> protocols = getAvailableSidreProtocols();
  for(size_t i = 0 ; i < protocols.size() ; ++i)
  {
    std::string& protocol = protocols[i];

    std::string file_path0 = file_path_tree0 + protocol;
    tree0->save(file_path0, protocol);

    Group* tree1 = tree0->createGroup("tree1");

    Group* gx = tree1->createGroup("x");
    Group* gy = tree1->createGroup("y");
    Group* gz = tree1->createGroup("z");

    gx->createViewAndAllocate("int20", DataType::int64(ndata*2));
    conduit::int64* data_ptr20 = gx->getView("int20")->getArray();
    for (int i = 0 ; i < ndata*2 ; ++i)
    {
      data_ptr20[i] = (conduit::int64)(i*2);
    }
    gy->createViewScalar<conduit::int64>("i0", 400);
    gz->createViewScalar<conduit::float64>("d0", 17.00);

    std::string file_path1 = file_path_tree1 + protocol;

    tree1->save(file_path1, protocol);

    // show the source tree
    SLIC_INFO("Source tree");
    ds.print();

    DataStore ds_load;
    Group* loadtree0 = ds_load.getRoot()->createGroup("tree0");
    loadtree0->load(file_path0, protocol);
    loadtree0->load(file_path1, protocol, true);

    SLIC_INFO("Tree from protocol: " << protocol);
    // show the result
    ds_load.print();

    Group* ds_load_root = ds_load.getRoot();

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView(
                "tree1/a/i0")->getData<conduit::int64>(),100);
    EXPECT_NEAR(ds_load_root->getView(
                  "tree1/a/d0")->getData<conduit::float64>(),3000.00,1e-12);
    EXPECT_EQ(ds_load_root->getView("tree1/b/s0")->getString(),
              std::string("foo"));
    EXPECT_EQ(ds_load_root->getView(
                "tree1/y/i0")->getData<conduit::int64>(),400);
    EXPECT_NEAR(ds_load_root->getView(
                  "tree1/z/d0")->getData<conduit::float64>(),17.00,1e-12);

    conduit::int64* load_data_ptr =
      ds_load_root->getView("tree1/c/int10")->getData();
    for(int j=0 ; j< ndata ; j++)
    {
      EXPECT_EQ(data_ptr[j],load_data_ptr[j]);
    }
    load_data_ptr = ds_load_root->getView("tree1/x/int20")->getData();
    for(int j=0 ; j< ndata*2 ; j++)
    {
      EXPECT_EQ(data_ptr20[j],load_data_ptr[j]);
    }

    // Destroy the group so the name can be reused by the next protocol
    tree0->destroyGroup("tree1");
  }

}
