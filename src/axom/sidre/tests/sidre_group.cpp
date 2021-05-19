// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"  // for AXOM_USE_HDF5
#include "axom/sidre/core/sidre.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

#include "gtest/gtest.h"

#include <cstring>
#include <vector>

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::FLOAT64_ID;
using axom::sidre::Group;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::INT_ID;
using axom::sidre::InvalidIndex;
using axom::sidre::nameIsValid;
using axom::sidre::TypeID;
using axom::sidre::View;

namespace
{
// Test protocols
#ifdef AXOM_USE_HDF5
int nprotocols = 3;
std::string const protocols[] = {"sidre_json", "sidre_hdf5", "json"};
#else
int nprotocols = 2;
std::string const protocols[] = {"sidre_json", "json"};
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

}  // end anonymous namespace

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(sidre_group, get_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("test");

  EXPECT_TRUE(root->isRoot());
  EXPECT_FALSE(group->isRoot());

  EXPECT_TRUE(group->getName() == std::string("test"));

  Group* group2 = root->getGroup("foo");
  EXPECT_TRUE(group2 == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// getPath(), getPathName()
//------------------------------------------------------------------------------
TEST(sidre_group, get_path_name)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  EXPECT_EQ(root->getParent(), root);
  const Group* croot = ds->getRoot();
  EXPECT_EQ(croot, root);
  EXPECT_EQ(root->getName(), "");
  Group* group = root->createGroup("test/a/b/c");
  Group* grp2 = root->getGroup("test/a");
  Group* grp3 = root->getGroup("test");

  EXPECT_EQ(root->getName(), std::string(""));
  EXPECT_EQ(root->getPath(), std::string(""));
  EXPECT_EQ(root->getPathName(), std::string(""));

  EXPECT_EQ(grp2->getName(), std::string("a"));
  EXPECT_EQ(grp2->getPath(), std::string("test"));
  EXPECT_EQ(grp2->getPathName(), std::string("test/a"));

  EXPECT_EQ(grp3->getName(), std::string("test"));
  EXPECT_EQ(grp3->getPath(), std::string(""));
  EXPECT_EQ(grp3->getPathName(), std::string("test"));

  EXPECT_EQ(group->getName(), std::string("c"));
  EXPECT_EQ(group->getPath(), std::string("test/a/b"));
  EXPECT_EQ(group->getPathName(), std::string("test/a/b/c"));

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup(), getGroup(), hasGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, group_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test full path access when building incrementally
  Group* group =
    root->createGroup("test1")->createGroup("test2")->createGroup("test3");
  Group* group2 = root->getGroup("test1/test2/test3");

  EXPECT_TRUE(nullptr != group2);
  EXPECT_EQ(group, group2);

  // Test incremental access when building full path
  Group* groupP = root->createGroup("testA/testB/testC");
  Group* groupP2 =
    root->getGroup("testA")->getGroup("testB")->getGroup("testC");

  EXPECT_TRUE(nullptr != groupP2);
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

  EXPECT_EQ(group_bada, static_cast<void*>(nullptr));
  EXPECT_EQ(group_badb, static_cast<void*>(nullptr));
  EXPECT_EQ(group_badc, static_cast<void*>(nullptr));

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

  EXPECT_EQ(group_cdup, static_cast<void*>(nullptr));
  EXPECT_EQ(group_testa->getGroup("testb")->getNumGroups(), testbnumgroups);

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup(), destroyGroup()  with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, destroy_group_with_path)
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
TEST(sidre_group, get_parent)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");

  EXPECT_TRUE(child->getParent() == parent);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(sidre_group, get_datastore)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  EXPECT_TRUE(group->getDataStore() == ds);

  DataStore const* const_ds = group->getDataStore();
  EXPECT_TRUE(const_ds == ds);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroup()
//------------------------------------------------------------------------------
TEST(sidre_group, get_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  Group* child = parent->createGroup("child");
  EXPECT_TRUE(child->getParent() == parent);

  EXPECT_TRUE(parent->getGroup("child") == child);
  EXPECT_TRUE(parent->getGroup(0) == child);

  // check error condition
  EXPECT_TRUE(parent->getGroup("non-existant group") == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getView()
//------------------------------------------------------------------------------
TEST(sidre_group, get_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  View* view = parent->createView("view");

  EXPECT_TRUE(parent->getView("view") == view);
  EXPECT_TRUE(parent->getView(0) == view);

  // check error condition
  EXPECT_TRUE(parent->getView("non-existant view") == nullptr);

  delete ds;
}

//------------------------------------------------------------------------------
// createView, hasView(), getView(), destroyView() with path strings
//------------------------------------------------------------------------------
TEST(sidre_group, view_with_path)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Test with full path access when building incrementally
  View* view =
    root->createGroup("group1")->createGroup("group2")->createView("view1");
  View* view2 = root->getView("group1/group2/view1");

  EXPECT_TRUE(nullptr != view2);
  EXPECT_EQ(view, view2);

  // Test incremental access when building with full path
  View* viewP = root->createView("groupA/groupB/viewA");
  View* viewP2 = root->getGroup("groupA")->getGroup("groupB")->getView("viewA");

  EXPECT_TRUE(nullptr != viewP2);
  EXPECT_EQ(viewP, viewP2);

  // Now verify that bad paths just return null, and don't create missing groups
  View* v_bad1 = root->getView("BAD/groupB/viewA");
  View* v_bad2 = root->getView("groupA/BAD/viewA");
  View* v_bad3 = root->getView("groupA/groupB/BAD");

  EXPECT_EQ(v_bad1, static_cast<void*>(nullptr));
  EXPECT_EQ(v_bad2, static_cast<void*>(nullptr));
  EXPECT_EQ(v_bad3, static_cast<void*>(nullptr));

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
            static_cast<void*>(nullptr));
  EXPECT_FALSE(root->hasView("group1/group2/view1"));
  EXPECT_EQ(root->getView("group1/group2/view1"), static_cast<void*>(nullptr));

  Group* groupA = root->getGroup("groupA");
  EXPECT_TRUE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), viewP);
  EXPECT_TRUE(root->hasView("groupA/groupB/viewA"));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), viewP);

  groupA->destroyView("groupB/viewA");

  EXPECT_EQ(exp_no_groups, groupA->getGroup("groupB")->getNumViews());
  EXPECT_FALSE(groupA->getGroup("groupB")->hasView("viewA"));
  EXPECT_EQ(groupA->getGroup("groupB")->getView("viewA"),
            static_cast<void*>(nullptr));
  EXPECT_FALSE(groupA->hasView("groupB/viewA"));
  EXPECT_EQ(groupA->getView("groupB/viewA"), static_cast<void*>(nullptr));
  EXPECT_EQ(root->getView("groupA/groupB/viewA"), static_cast<void*>(nullptr));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getViewName(), getViewIndex()
//------------------------------------------------------------------------------
TEST(sidre_group, get_view_names_and_indicies)
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
  EXPECT_FALSE(indexIsValid(idx3));

  const std::string& name3 = parent->getViewName(idx3);
  EXPECT_TRUE(name3.empty());
  EXPECT_FALSE(nameIsValid(name3));

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getFirstValidGroupIndex, getNextValidGroupIndex
//------------------------------------------------------------------------------
TEST(sidre_group, get_first_and_next_group_index)
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

  EXPECT_FALSE(indexIsValid(badidx1));
  EXPECT_FALSE(indexIsValid(badidx2));

  delete ds;
}

TEST(sidre_group, child_lists)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // parent is a Group in list format.
  Group* parent = root->createGroup("parent", true);

  // Create 10 unnamed Groups as children of parent.
  for(IndexType i = 0; i < 10; ++i)
  {
    Group* unnamed_group = parent->createUnnamedGroup();
    unnamed_group->createViewScalar("val", i);
  }

  // Create 15 unnamed Views as children of parent.
  for(IndexType i = 0; i < 15; ++i)
  {
    View* unnamed_view;
    if(i % 3 == 0)
    {
      unnamed_view = parent->createView("");
    }
    else if(i % 3 == 1)
    {
      unnamed_view = parent->createViewScalar("", i * i);
    }
    else
    {
      unnamed_view = parent->createViewString("", "foo");
    }
    if(!unnamed_view->isApplied())
    {
      unnamed_view->apply(INT_ID, i);
      unnamed_view->allocate(INT_ID, i);
      int* vdata = unnamed_view->getData();
      for(IndexType j = 0; j < i; ++j)
      {
        vdata[j] = j + 3;
      }
    }
  }

  // Create Group not in list format, show that it can't create
  // unnamed children.
  Group* not_list = root->createGroup("not_list", false);
  Group* dummy_group = not_list->createUnnamedGroup();
  View* dummy_view = not_list->createView("");
  EXPECT_TRUE(not_list->isUsingMap());
  EXPECT_EQ(not_list->getNumGroups(), 0);
  EXPECT_EQ(not_list->getNumViews(), 0);
  EXPECT_TRUE(dummy_group == nullptr);
  EXPECT_TRUE(dummy_view == nullptr);

  // Access data from unnamed Groups held by parent.
  std::set<int> scalars;
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    Group* unnamed_group = parent->getGroup(idx);
    View* val_view = unnamed_group->getView("val");
    IndexType val = val_view->getScalar();
    EXPECT_TRUE(val >= 0 && val < 10);

    scalars.insert(val);
  }

  EXPECT_TRUE(parent->isUsingList());
  EXPECT_TRUE(scalars.size() == 10);
  EXPECT_TRUE(parent->hasGroup(6));
  EXPECT_FALSE(parent->hasGroup(20));

  // Destroy five of the unnamed Groups held by parent.
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    if(idx % 2 == 1)
    {
      parent->destroyGroup(idx);
    }
  }

  // Add one more unnamed Group, so there should be six child Groups.
  (void)parent->createUnnamedGroup();
  EXPECT_EQ(parent->getNumGroups(), 6);

  // Access data from the unnamed Views.
  for(IndexType idx = parent->getFirstValidViewIndex(); indexIsValid(idx);
      idx = parent->getNextValidViewIndex(idx))
  {
    View* unnamed_view = parent->getView(idx);
    if(idx % 3 == 0)
    {
      EXPECT_EQ(unnamed_view->getTypeID(), INT_ID);
      IndexType num_elems = unnamed_view->getNumElements();
      EXPECT_EQ(num_elems, idx);
      int* vdata = unnamed_view->getData();
      for(IndexType j = 0; j < num_elems; ++j)
      {
        EXPECT_EQ(vdata[j], j + 3);
      }
    }
    else if(idx % 3 == 1)
    {
      EXPECT_TRUE(unnamed_view->isScalar());
      IndexType val = unnamed_view->getScalar();
      EXPECT_EQ(val, idx * idx);
    }
    else
    {
      EXPECT_TRUE(unnamed_view->isString());
      std::string vstr = unnamed_view->getString();
      EXPECT_EQ(vstr, std::string("foo"));
    }
  }

  root->destroyGroup("parent");

  delete ds;
}

//------------------------------------------------------------------------------
// Iterate Views with getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, iterate_groups)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  (void)parent->createView("view1");
  (void)parent->createView("view2");
  (void)parent->createView("view3");
  (void)parent->createGroup("g1");
  (void)parent->createView("view4");
  (void)parent->createView("view5");
  (void)parent->createView("view6");
  (void)parent->createGroup("g2");
  (void)parent->createGroup("g3");
  (void)parent->createView("view7");
  (void)parent->createView("view8");
  (void)parent->createView("view9");
  (void)parent->createGroup("g4");

  int groupcount = 0;
  for(IndexType idx = parent->getFirstValidGroupIndex(); indexIsValid(idx);
      idx = parent->getNextValidGroupIndex(idx))
  {
    groupcount += 1;
  }
  EXPECT_EQ(4, groupcount);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, get_first_and_next_view_index)
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

  EXPECT_FALSE(indexIsValid(badidx1));
  EXPECT_FALSE(indexIsValid(badidx2));

  delete ds;
}

//------------------------------------------------------------------------------
// Iterate Views with getFirstValidViewIndex, getNextValidViewIndex
//------------------------------------------------------------------------------
TEST(sidre_group, iterate_views)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* parent = root->createGroup("parent");
  (void)parent->createView("view1");
  (void)parent->createView("view2");
  (void)parent->createView("view3");
  (void)parent->createGroup("g1");
  (void)parent->createView("view4");
  (void)parent->createView("view5");
  (void)parent->createView("view6");
  (void)parent->createGroup("g2");
  (void)parent->createGroup("g3");
  (void)parent->createView("view7");
  (void)parent->createView("view8");
  (void)parent->createView("view9");
  (void)parent->createGroup("g4");

  int viewcount = 0;
  IndexType idx = parent->getFirstValidViewIndex();
  while(indexIsValid(idx))
  {
    viewcount += 1;

    idx = parent->getNextValidViewIndex(idx);
  }
  EXPECT_EQ(9, viewcount);

  delete ds;
}

//------------------------------------------------------------------------------
// Verify getGroupName(), getGroupIndex()
//------------------------------------------------------------------------------
TEST(sidre_group, get_group_name_index)
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
  EXPECT_FALSE(indexIsValid(idx3));

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
TEST(sidre_group, create_destroy_has_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("parent");

  View* view = group->createView("view");
  EXPECT_TRUE(group->getParent() == root);
  EXPECT_FALSE(view->hasBuffer());
  EXPECT_TRUE(group->hasView("view"));
  IndexType iview = group->getViewIndex("view");
  EXPECT_EQ(0, iview);
  iview = view->getIndex();
  EXPECT_EQ(0, iview);

  // try creating view again, should be a no-op.
  EXPECT_TRUE(group->createView("view") == nullptr);

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
  EXPECT_FALSE(group->hasView("view"));

  // try api call that specifies specific type and length
  group->createViewAndAllocate("viewWithLength1", INT_ID, 50);
  IndexType iview2 = group->getViewIndex("viewWithLength1");
  EXPECT_EQ(iview, iview2);  // reuse slot

  // error condition check - try again with duplicate name, should be a no-op
  EXPECT_TRUE(group->createViewAndAllocate("viewWithLength1", FLOAT64_ID, 50) ==
              nullptr);
  group->destroyViewAndData("viewWithLength1");
  EXPECT_FALSE(group->hasView("viewWithLength1"));

  EXPECT_TRUE(group->createViewAndAllocate("viewWithLengthBadLen",
                                           FLOAT64_ID,
                                           -1) == nullptr);

  // try api call that specifies data type in another way
  group->createViewAndAllocate("viewWithLength2", DataType::float64(50));
  EXPECT_TRUE(group->createViewAndAllocate("viewWithLength2",
                                           DataType::float64(50)) == nullptr);
  // destroy view and its buffer using index
  IndexType indx = group->getFirstValidViewIndex();
  IndexType bindx = group->getView(indx)->getBuffer()->getIndex();
  group->destroyViewAndData(indx);
  EXPECT_TRUE(ds->getBuffer(bindx) == NULL);

  // Destroy view but not the buffer
  view = group->createViewAndAllocate("viewWithLength2", INT_ID, 50);
  Buffer* buff = view->getBuffer();
  group->destroyView("viewWithLength2");
  EXPECT_TRUE(buff->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------
// createViewAndAllocate() with zero-sized array
//------------------------------------------------------------------------------
TEST(sidre_group, create_zero_sized_view)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  View* zeroSizedView = root->createViewAndAllocate("foo", INT_ID, 0);
  EXPECT_TRUE(zeroSizedView->isDescribed());
  EXPECT_TRUE(zeroSizedView->isAllocated());

  delete ds;
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_has_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* group = root->createGroup("group");
  EXPECT_TRUE(group->getParent() == root);

  EXPECT_TRUE(root->hasGroup("group"));

  root->destroyGroup("group");
  EXPECT_FALSE(root->hasGroup("group"));

  // should be a no-op, not a failure
  root->destroyGroup("group");

  Group* group2 = root->createGroup("group2");
  // shut up compiler about unused variable
  (void)group2;
  root->destroyGroup(root->getFirstValidGroupIndex());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, group_name_collisions)
{
  DataStore* ds = new DataStore();
  Group* flds = ds->getRoot()->createGroup("fields");
  flds->createView("a");

  EXPECT_TRUE(flds->hasChildView("a"));

  // attempt to create duplicate group name
  EXPECT_TRUE(ds->getRoot()->createGroup("fields") == nullptr);

  // attempt to create duplicate view name.
  EXPECT_TRUE(flds->createView("a") == nullptr);

  // attempt to create a group named the same as an existing view
  EXPECT_TRUE(flds->createGroup("a") == nullptr);

  // attempt to create a view named the same as an existing group
  EXPECT_TRUE(ds->getRoot()->createView("fields") == nullptr);

  ds->getRoot()->createGroup("here//is/path");
  ds->getRoot()->createGroup("éch≈o/Ωd");
  ds->getRoot()->createGroup("../group/..");

  IndexType idx = ds->getRoot()->getFirstValidGroupIndex();
  while(indexIsValid(idx))
  {
    std::cout << ds->getRoot()->getGroup(idx)->getName() << std::endl;
    idx = ds->getRoot()->getNextValidGroupIndex(idx);
  }

  delete ds;
}
//------------------------------------------------------------------------------

TEST(sidre_group, view_copy_move)
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
  views[3] =
    flds->createView("external", INT_ID, 10)->setExternalDataPtr(extdata);
  views[4] = flds->createViewScalar("scalar", 25);
  views[5] = flds->createViewString("string", "I am string");

  buffdata = flds->getView("buffer")->getData();
  for(int i = 0; i < 10; ++i)
  {
    extdata[i] = i;
    buffdata[i] = i + 100;
  }

  for(int i = 0; i < 6; ++i)
  {
    names[i] = views[i]->getName();
    EXPECT_TRUE(flds->hasView(names[i]));
  }

  // test moving a view from flds to sub1
  Group* sub1 = flds->createGroup("sub1");

  // flds->print();

  for(int i = 0; i < 6; ++i)
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

  for(int i = 0; i < 6; ++i)
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
  EXPECT_EQ(view1->getData<int>(), view2->getData<int>());

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
TEST(sidre_group, groups_move_copy)
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
  if(gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  EXPECT_EQ(bschild->getNumViews(), 1);
  EXPECT_TRUE(bschild->hasChildView("val"));
  if(gb0->hasChildView("val"))
  {
    EXPECT_EQ((double)gb0->getView("val")->getScalar(), val);
  }

  EXPECT_EQ(flds->getNumGroups(), 3);
  EXPECT_TRUE(flds->hasGroup("a"));
  EXPECT_TRUE(flds->hasGroup("sub"));
  EXPECT_TRUE(flds->hasGroup("c"));

  EXPECT_EQ(flds->getGroup("sub")->getGroup("b"), gb);

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
  EXPECT_EQ(triedCopy, static_cast<void*>(nullptr));
  EXPECT_EQ(gsub->getNumGroups(), 1);
  EXPECT_TRUE(gsub->hasChildGroup("b"));
  EXPECT_EQ(buffercount, ds->getNumBuffers());

  EXPECT_EQ(gb0->getNumGroups(), 1);
  EXPECT_EQ(gb0->getGroup("childOfB"), bschild);
  EXPECT_EQ(gb0->getNumViews(), 1);
  EXPECT_TRUE(gb0->hasChildView("f0"));
  if(gb0->hasChildView("f0"))
  {
    EXPECT_EQ((double)gb0->getView("f0")->getScalar(), f0value);
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_view_and_buffer2)
{
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string viewName1("viewBuffer1");
  std::string viewName2("viewBuffer2");

  View* view1 = grp->createViewAndAllocate(viewName1, INT_ID, 1);
  View* view2 = grp->createViewAndAllocate(viewName2, INT_ID, 1);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ(grp->getView(viewName1), view1);

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ(grp->getView(viewName2), view2);

  IndexType const bufferId1 = view1->getBuffer()->getIndex();

  grp->destroyViewAndData(viewName1);

  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumBuffers(), 1u);

  Buffer const* const buffer1 = ds->getBuffer(bufferId1);
  EXPECT_TRUE(buffer1 == nullptr);

  View const* const view3 = grp->createView("viewBuffer3");
  grp->destroyViewsAndData();
  // should be no-op
  grp->destroyViewsAndData();
  // shut up compiler about unused variable
  (void)view3;

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_destroy_alloc_view_and_buffer)
{
  DataStore* const ds = new DataStore();
  Group* const grp = ds->getRoot()->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  // use create + alloc convenience methods
  // this one is the DataType & method
  View* const view1 = grp->createViewAndAllocate(viewName1, DataType::c_int(10));

  EXPECT_TRUE(grp->hasChildView(viewName1));
  EXPECT_EQ(grp->getView(viewName1), view1);

  int* v1_vals = view1->getData();

  for(int i = 0; i < 10; i++)
  {
    v1_vals[i] = i;
  }

  EXPECT_EQ(view1->getNumElements(), 10u);
  EXPECT_EQ(view1->getTotalBytes(),
            static_cast<axom::sidre::IndexType>(10 * sizeof(int)));

  grp->destroyViewAndData(viewName1);

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, create_view_of_buffer_with_schema)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  // use create + alloc convenience methods
  // this one is the DataType & method
  View* base = root->createViewAndAllocate("base", DataType::c_int(10));
  int* base_vals = base->getData();
  for(int i = 0; i < 10; i++)
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

  for(int i = 0; i < 5; i++)
  {
    EXPECT_EQ(sub_a_vals[i], 10);
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_empty_datastore)
{
  const std::string file_path_base("sidre_empty_datastore_");
  DataStore* ds1 = new DataStore();

  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    ds1->getRoot()->save(file_path, protocols[i]);
  }

  delete ds1;

  // Only restore conduit_hdf5
  for(int i = 1; i < 2; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();
    Group* root2 = ds2->getRoot();
    root2->load(file_path, protocols[i]);

    EXPECT_TRUE(ds2->getNumBuffers() == 0);
    EXPECT_TRUE(root2->getNumGroups() == 0);
    EXPECT_TRUE(root2->getNumViews() == 0);

    delete ds2;
  }
}

#ifdef AXOM_USE_HDF5
//------------------------------------------------------------------------------
// make sure the hdf5 methods are consistent with the path based methods
//------------------------------------------------------------------------------
TEST(sidre_group, save_load_via_hdf5_ids)
{
  DataStore ds_save;
  // populate the datastore
  Group* root = ds_save.getRoot();
  root->createViewScalar<int>("i0", 1);
  root->createViewAndAllocate("vals", INT_ID, 5);
  // set values for the "vals" array
  int* vals_ptr = root->getView("vals")->getData();
  for(int i = 0; i < 5; ++i)
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

  hid_t h5_id =
    H5Fopen("out_save_load_via_hdf5_ids.sidre_hdf5", H5F_ACC_RDWR, H5P_DEFAULT);
  EXPECT_TRUE(h5_id >= 0);

  // this implies protocol == "sidre_hdf5"
  ds_load_hdf5.getRoot()->load(h5_id);

  // ? Does isEquivalentTo check values?
  // check path based with source
  EXPECT_TRUE(ds_load_generic.getRoot()->isEquivalentTo(ds_save.getRoot()));

  // check hdf5 based with source
  EXPECT_TRUE(ds_load_hdf5.getRoot()->isEquivalentTo(ds_save.getRoot()));

  // check path based vs hdf5 based
  EXPECT_TRUE(ds_load_generic.getRoot()->isEquivalentTo(ds_load_hdf5.getRoot()));

  // close hdf5 handle
  EXPECT_TRUE(H5Fclose(h5_id) >= 0);
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_root_restore_as_child)
{
  // We'll save the DataStore's root Group into a file then restore it
  // into a child group of another DataStore.

  // Create a DataStore; put some groups and views into it
  const std::string file_path_base("sidre_save_root_restore_as_child_");
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* child1 = root->createGroup("g_a");
  Group* child2 = root->createGroup("g_b");
  root->createGroup("g_c");  // We don't put anything into g_c

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  const int num_elems = 4;
  int* pa0 = child2->createViewAndAllocate("a0", INT_ID, num_elems)->getArray();
  double* pa1 =
    child2->createViewAndAllocate("a1", FLOAT64_ID, num_elems)->getArray();

  const double factor = 2.3;
  const double offset = -0.23;
  for(int i = 0; i < num_elems; ++i)
  {
    pa0[i] = i;
    pa1[i] = offset + i * factor;
  }

  // Save the DataStore's root
  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root->save(file_path, protocols[i]);
  }

  // Restore the original DataStore into a child group
  for(int i = 0; i < nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    DataStore* dscopy = new DataStore();
    Group* dsroot = dscopy->getRoot();
    const std::string file_path = file_path_base + protocols[i];

    const std::string group_base("group_");
    const std::string group_name = group_base + protocols[i];
    Group* cg = dsroot->createGroup(group_name);

    if(axom::utilities::filesystem::pathExists(file_path))
    {
      std::cout << "loading " << file_path << std::endl;
      cg->load(file_path, protocols[i]);

      EXPECT_TRUE(cg->isEquivalentTo(root, false));
      EXPECT_TRUE(root->isEquivalentTo(cg, false));
    }
    else
    {
      FAIL() << "file not present: " << file_path;
    }

    delete dscopy;
  }

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_child_restore_as_root)
{
  // We'll save a child Group into a file then restore it as the root of
  // a new DataStore.

  // Create a DataStore; put some groups and views into it
  const std::string file_path_base("sidre_save_child_restore_as_root_");
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  // child1 and all its descendents will get saved into files
  Group* child1 = root->createGroup("g_a");
  Group* child2 = child1->createGroup("g_b");
  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");
  // everything else won't be put into the files
  Group* child1a = root->createGroup("g_notSaved");
  child1a->createViewScalar<int>("i1", 42);
  root->createViewString("s1", "string view off the root, not saved");

  const int num_elems = 4;
  // included in files
  int* pa0 = child2->createViewAndAllocate("a0", INT_ID, num_elems)->getArray();
  double* pa1 =
    child2->createViewAndAllocate("a1", FLOAT64_ID, num_elems)->getArray();
  // not included
  int* pa2 = child1a->createViewAndAllocate("a2", INT_ID, num_elems)->getArray();
  int* pa3 = root->createViewAndAllocate("a3", INT_ID, num_elems)->getArray();

  const double factor = 2.3;
  const double offset = -0.23;
  for(int i = 0; i < num_elems; ++i)
  {
    pa0[i] = i;
    pa1[i] = offset + i * factor;
    pa2[i] = i + 2;
    pa3[i] = 4 - i;
  }

  // Save the Group in question (child1) into an archive
  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    child1->save(file_path, protocols[i]);
  }

  // Restore the saved child1 into a root group
  for(int i = 0; i < nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    DataStore* dscopy = new DataStore();
    const std::string file_path = file_path_base + protocols[i];
    if(axom::utilities::filesystem::pathExists(file_path))
    {
      dscopy->getRoot()->load(file_path, protocols[i]);

      EXPECT_TRUE(dscopy->getRoot()->isEquivalentTo(child1, false));
      EXPECT_TRUE(child1->isEquivalentTo(dscopy->getRoot(), false));
    }
    else
    {
      FAIL() << "file not present: " << file_path;
    }

    delete dscopy;
  }

  delete ds;
}
#endif  // AXOM_USE_HDF5

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_api)
{
  const std::string file_path_base("sidre_save_subtree_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);

  // These should be produce identical files.

  // No group provided, defaults to root group
  root1->save("sidre_save_fulltree_conduit", "json");

  for(int i = 0; i < nprotocols; ++i)
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
  EXPECT_TRUE(ds4->getRoot()->isEquivalentTo(root1));
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

  EXPECT_TRUE(load1->isEquivalentTo(load2));

  std::string newgroupname = "in case of blank";
  std::string groupname = newgroupname;
  bool loadSuccess = false;
  Group* load3 = load2->createGroupAndLoad(groupname,
                                           "sidre_save_subtree_sidre_json",
                                           "sidre_json",
                                           loadSuccess);

  EXPECT_NE(load3, (Group*)nullptr);
  EXPECT_TRUE(loadSuccess);
  EXPECT_EQ(groupname, "");
  if(load3 != nullptr)
  {
    EXPECT_EQ(newgroupname, load3->getName());
    EXPECT_TRUE(load1->isEquivalentTo(load3, false));
  }

  std::string anothergroupname = "another group";
  groupname = anothergroupname;
  loadSuccess = false;
  Group* load4 = load2->createGroupAndLoad(groupname,
                                           "sidre_save_subtree_sidre_json",
                                           "sidre_json",
                                           loadSuccess);

  EXPECT_NE(load4, (Group*)nullptr);
  EXPECT_TRUE(loadSuccess);
  if(load4 != nullptr)
  {
    EXPECT_EQ(anothergroupname, load4->getName());
    EXPECT_TRUE(load3->isEquivalentTo(load4, false));
  }

  groupname = anothergroupname;
  loadSuccess = false;
  Group* load5 = load2->createGroupAndLoad(groupname,
                                           "sidre_save_subtree_sidre_json",
                                           "sidre_json",
                                           loadSuccess);
  EXPECT_EQ(load5, (Group*)nullptr);
  EXPECT_FALSE(loadSuccess);

  delete ds_new;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_scalars_and_strings)
{
  const std::string file_path_base("sidre_save_scalars_and_strings_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createViewScalar<int>("i0", 1);
  root1->createViewScalar<float>("f0", 1.0);
  root1->createViewScalar<double>("d0", 10.0);
  root1->createViewString("s0", "I am a string");

  for(int i = 0; i < nprotocols; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  for(int i = 0; i < nprotocols; ++i)
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

    EXPECT_TRUE(root1->isEquivalentTo(root2));

    int i0 = root2->getView("i0")->getScalar();
    float f0 = root2->getView("f0")->getScalar();
    double d0 = root2->getView("d0")->getScalar();
    const char* s0 = root2->getView("s0")->getString();

    EXPECT_EQ(1, i0);
    EXPECT_EQ(1.0, f0);
    EXPECT_EQ(10.0, d0);
    EXPECT_EQ(std::string(s0), "I am a string");

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group, rename_group)
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();
  Group* child1 = root->createGroup("g_a");
  Group* child2 = root->createGroup("g_b");
  Group* child3 = root->createGroup("g_c");

  // rename should not change the index
  EXPECT_EQ(0, child1->getIndex());
  bool success = child1->rename("g_r");
  EXPECT_TRUE(success);
  EXPECT_EQ("g_r", child1->getName());
  EXPECT_EQ(0, child1->getIndex());
  EXPECT_TRUE(root->hasGroup("g_r"));
  EXPECT_FALSE(root->hasGroup("g_a"));

  // try to rename to path
  success = child2->rename("fields/g_s");
  EXPECT_FALSE(success);
  EXPECT_EQ("g_b", child2->getName());

  // Try to rename to existing group name
  success = child3->rename("g_b");
  EXPECT_FALSE(success);
  EXPECT_EQ("g_c", child3->getName());

  // Rename root group
  EXPECT_FALSE(indexIsValid(root->getIndex()));
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("", root->getName());
  root->rename("newroot");
  EXPECT_FALSE(indexIsValid(root->getIndex()));
  EXPECT_EQ(root, root->getParent());
  EXPECT_EQ("newroot", root->getName());

  delete ds;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_name_change)
{
  const std::string file_path_base("sidre_save_name_change_");
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Group* child1 = root1->createGroup("child1");

  child1->createViewScalar<int>("i0", 1);
  child1->createViewString("s0", "I am a string");

  bool success = child1->getView("s0")->rename("s0_renamed");

  EXPECT_TRUE(success);
  EXPECT_FALSE(child1->hasView("s0"));
  EXPECT_TRUE(child1->hasView("s0_renamed"));

  for(int i = 0; i < nprotocols; ++i)
  {
    //      if ( protocols[i] == "conduit_hdf5")
    //    continue;   // XXX - Does not work
    const std::string file_path = file_path_base + protocols[i];
    child1->save(file_path, protocols[i]);
  }

  std::string groupname;
  for(int i = 0; i < nprotocols; ++i)
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

    EXPECT_EQ(child2->getName(), "child2");

    child2->load(file_path, protocols[i], false, groupname);

    EXPECT_EQ(child2->getName(), "child2");

    child2->rename(groupname);

    EXPECT_TRUE(root1->isEquivalentTo(root2));

    int i0 = child2->getView("i0")->getScalar();
    const char* s0 = child2->getView("s0_renamed")->getString();

    EXPECT_EQ(1, i0);
    EXPECT_EQ(std::string(s0), "I am a string");

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_external_data)
{
  const std::string file_path_base("sidre_save_external_");

  const int nfoo = 10;
  int foo1[nfoo], foo2[nfoo], *foo3, foo4[nfoo];
  int int2d1[nfoo * 2], int2d2[nfoo * 2];
  IndexType shape[] = {nfoo, 2};

  for(int i = 0; i < nfoo; ++i)
  {
    foo1[i] = i;
    foo2[i] = 0;
    foo4[i] = i;
  }
  for(int i = 0; i < 2 * nfoo; ++i)
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

  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  delete ds1;

  // Now load back in.
  for(int i = 0; i < nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];
    IndexType extents[7];
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
    EXPECT_TRUE(view1->getVoidPtr() == nullptr);
    view1->setExternalDataPtr(foo2);

    View* view2 = root2->getView("empty_array");
    EXPECT_TRUE(view2->isEmpty());
    EXPECT_TRUE(view2->isDescribed());
    EXPECT_EQ(view2->getTypeID(), INT_ID);
    EXPECT_TRUE(view2->getVoidPtr() == nullptr);
    view2->setExternalDataPtr(foo3);

    View* view3 = root2->getView("external_undescribed");
    EXPECT_TRUE(view3->isEmpty());
    EXPECT_FALSE(view3->isDescribed());
    EXPECT_TRUE(view3->getVoidPtr() == nullptr);
    // Set "external_array" and "external_undescribed" to the same external
    // array
    // since it was created that way.  However, "external_undescribed" was not
    // written to the dump since it is undescribed.
    view3->setExternalDataPtr(foo2);

    View* view4 = root2->getView("int2d");
    EXPECT_FALSE(view4->isEmpty());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_TRUE(view4->getVoidPtr() == nullptr);
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), nfoo * 2);
    EXPECT_EQ(view4->getNumDimensions(), 2);
    rank = view4->getShape(7, extents);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(extents[0] == nfoo && extents[1] == 2);
    view4->setExternalDataPtr(int2d2);

    // Read external data into views
    root2->loadExternalData(file_path);

    // Make sure addresses have not changed
    EXPECT_TRUE(view1->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view2->getVoidPtr() == static_cast<void*>(foo3));  // nullptr
    EXPECT_TRUE(view3->getVoidPtr() == static_cast<void*>(foo2));
    EXPECT_TRUE(view4->getVoidPtr() == static_cast<void*>(int2d2));

    for(int j = 0; i < nfoo; ++i)
    {
      EXPECT_TRUE(foo1[j] == foo2[j]);
    }
    for(int j = 0; i < 2 * nfoo; ++i)
    {
      EXPECT_TRUE(int2d1[j] == int2d2[j]);
    }

    delete ds2;
  }
}

//------------------------------------------------------------------------------

// Check the association between views and buffers to make sure it is what we
// expect.
// This checks more than isEquivalentTo.

static void save_restore_buffer_association(const std::string& msg, DataStore* ds)
{
  const IndexType len = 10;

  SCOPED_TRACE(msg);

  // Make sure all buffers were created
  ASSERT_EQ(ds->getNumBuffers(), 4);

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
  for(int ii = 0; ii < len; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 100);
  }

  ASSERT_EQ(buff4a->getNumElements(), len);
  idata = buff4a->getData();
  for(int ii = 0; ii < len; ++ii)
  {
    ASSERT_EQ(idata[ii], ii + 200);
  }
}

TEST(sidre_group, save_restore_buffer)
{
  const std::string file_path_base("sidre_save_buffer_");
  const IndexType len = 10;

  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();
  Buffer* buff1 = ds1->createBuffer();
  Buffer* buff2 = ds1->createBuffer(INT_ID, len);
  Buffer* buff3 = ds1->createBuffer(INT_ID, len)->allocate();

  int* idata = buff3->getData();
  for(IndexType ii = 0; ii < len; ++ii)
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
  for(int ii = 0; ii < len; ++ii)
  {
    idata[ii] = ii + 200;
  }

  save_restore_buffer_association("original datastore", ds1);

  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  // Now load back in.
  for(int i = 0; i < nprotocols; ++i)
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
    EXPECT_TRUE(isequivalent);
    if(isequivalent)
    {
      save_restore_buffer_association("loaded datastore", ds2);
    }

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------

TEST(sidre_group, save_restore_other)
{
  const std::string file_path_base("sidre_save_other_");
  const int ndata = 10;
  IndexType shape1[] = {ndata, 2};
  DataStore* ds1 = new DataStore();
  Group* root1 = ds1->getRoot();

  root1->createView("empty_view");
  root1->createView("empty_described", INT_ID, ndata);
  root1->createView("empty_shape", INT_ID, 2, shape1);

  auto* view = root1->createViewAndAllocate("buffer_shape", INT_ID, 2, shape1);

  // fill the data
  {
    int* data = view->getData();
    for(int i = 0; i < ndata * 2; ++i)
    {
      data[i] = i;
    }
  }

  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  delete ds1;

  // Now load back in.
  for(int i = 0; i < nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];
    IndexType shape2[7];
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
    EXPECT_EQ(view3->getNumElements(), ndata * 2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view3->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    View* view4 = root2->getView("buffer_shape");
    EXPECT_TRUE(view4->hasBuffer());
    EXPECT_TRUE(view4->isDescribed());
    EXPECT_EQ(view4->getTypeID(), INT_ID);
    EXPECT_EQ(view4->getNumElements(), ndata * 2);
    shape2[0] = 0;
    shape2[1] = 0;
    rank = view4->getShape(7, shape2);
    EXPECT_EQ(rank, 2);
    EXPECT_TRUE(shape2[0] == ndata && shape2[1] == 2);

    delete ds2;
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_restore_complex)
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
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = i;
  }

  for(int i = 0; i < nprotocols; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    ds1->getRoot()->save(file_path, protocols[i]);
  }

  for(int i = 0; i < nprotocols; ++i)
  {
    // Only restore sidre_hdf5 protocol
    if(protocols[i] != "sidre_hdf5")
    {
      continue;
    }

    const std::string file_path = file_path_base + protocols[i];

    DataStore* ds2 = new DataStore();

    ds2->getRoot()->load(file_path, protocols[i]);

    EXPECT_TRUE(ds1->getRoot()->isEquivalentTo(ds2->getRoot()));

    flds = ds2->getRoot()->getGroup("fields");

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getData<int>(), 100.0);
    EXPECT_NEAR(flds->getGroup("a")->getView("d0")->getData<double>(),
                3000.0,
                1e-12);

    int* new_data_ptr = flds->getGroup("c")->getView("int10")->getArray();
    for(int i = 0; i < ndata; ++i)
    {
      EXPECT_TRUE(new_data_ptr[i] == i);
    }

    const char* char_ptr = flds->getView("b/s0")->getString();
    EXPECT_TRUE(std::string(char_ptr) == "foo");

    //ds2->print();

    delete ds2;
  }

  delete ds1;
}

//------------------------------------------------------------------------------
// isEquivalentTo()
//------------------------------------------------------------------------------
TEST(sidre_group, is_equivalent_to)
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

  Group* gc2 = flds2->createGroup("c");  // Note: flds2 groups added in
                                         // different order
  Group* gb2 = flds2->createGroup("b");
  Group* ga2 = flds2->createGroup("a");

  ga1->createViewScalar("i0", 1);
  gb1->createViewScalar("f0", 100.0f);
  gc1->createViewScalar("d0", 3000.00);
  gc1->createViewScalar("d1", 6000.00);
  gc1->createViewScalar("d2", 9000.00);

  ga2->createViewScalar("i0", 1);
  gb2->createViewScalar("f0", 100.0f);
  gc2->createViewScalar("d2", 9000.00);  // Note: views of gc2 added in
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
TEST(sidre_group, save_load_all_protocols)
{
  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

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
  for(int i = 0; i < ndata; ++i)
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
  for(size_t i = 0; i < protocols.size(); ++i)
  {
    SLIC_INFO("Testing protocol: " << protocols[i]);
    const std::string file_path = file_path_base + protocols[i];
    // save using current protocol
    ds.getRoot()->save(file_path, protocols[i]);

    DataStore ds_load;
    ds_load.getRoot()->load(file_path, protocols[i]);

    SLIC_INFO("Tree from protocol: " << protocols[i]);
    // show the result
    ds_load.print();

    Group* ds_load_root = ds_load.getRoot();
    // check that the sidre hierarchy is equiv
    EXPECT_TRUE(ds.getRoot()->isEquivalentTo(ds_load_root));

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView("fields/a/i0")->getData<conduit::int64>(),
              100);
    EXPECT_NEAR(ds_load_root->getView("fields/a/d0")->getData<conduit::float64>(),
                3000.00,
                1e-12);
    EXPECT_EQ(ds_load_root->getView("fields/b/s0")->getString(),
              std::string("foo"));

    conduit::int64* load_data_ptr =
      ds_load_root->getView("fields/c/int10")->getData();
    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(data_ptr[j], load_data_ptr[j]);
    }
  }

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_group, save_load_preserve_contents)
{
  // Note: This test relies on re-wiring conduit error handlers
  DataStore::setConduitSLICMessageHandlers();

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
  for(int i = 0; i < ndata; ++i)
  {
    data_ptr[i] = (conduit::int64)i;
  }

  std::vector<std::string> protocols = getAvailableSidreProtocols();
  std::string groupname;
  for(size_t i = 0; i < protocols.size(); ++i)
  {
    std::string& protocol = protocols[i];

    std::string file_path0 = file_path_tree0 + protocol;
    tree0->save(file_path0, protocol);

    Group* tree1 = tree0->createGroup("tree1");

    Group* gx = tree1->createGroup("x");
    Group* gy = tree1->createGroup("y");
    Group* gz = tree1->createGroup("z");

    gx->createViewAndAllocate("int20", DataType::int64(ndata * 2));
    conduit::int64* data_ptr20 = gx->getView("int20")->getArray();
    for(int i = 0; i < ndata * 2; ++i)
    {
      data_ptr20[i] = (conduit::int64)(i * 2);
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
    loadtree0->load(file_path0, protocol, false, groupname);
    loadtree0->load(file_path1, protocol, true, groupname);
    loadtree0->rename(groupname);

    SLIC_INFO("Tree from protocol: " << protocol);
    // show the result
    ds_load.print();

    Group* ds_load_root = ds_load.getRoot();

    // check that the values are the same
    EXPECT_EQ(ds_load_root->getView("tree1/a/i0")->getData<conduit::int64>(),
              100);
    EXPECT_NEAR(ds_load_root->getView("tree1/a/d0")->getData<conduit::float64>(),
                3000.00,
                1e-12);
    EXPECT_EQ(ds_load_root->getView("tree1/b/s0")->getString(),
              std::string("foo"));
    EXPECT_EQ(ds_load_root->getView("tree1/y/i0")->getData<conduit::int64>(),
              400);
    EXPECT_NEAR(ds_load_root->getView("tree1/z/d0")->getData<conduit::float64>(),
                17.00,
                1e-12);

    conduit::int64* load_data_ptr =
      ds_load_root->getView("tree1/c/int10")->getData();
    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(data_ptr[j], load_data_ptr[j]);
    }
    load_data_ptr = ds_load_root->getView("tree1/x/int20")->getData();
    for(int j = 0; j < ndata * 2; j++)
    {
      EXPECT_EQ(data_ptr20[j], load_data_ptr[j]);
    }

    // Destroy the group so the name can be reused by the next protocol
    tree0->destroyGroup("tree1");
  }

  // restore conduit default errors
  DataStore::setConduitDefaultMessageHandlers();
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit)
{
  conduit::Node input;

  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  DataStore ds;

  EXPECT_TRUE(ds.getRoot()->importConduitTree(input));

  EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
  EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(),
              3000.00,
              1e-12);
  EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(),
            std::string("foo"));

  conduit::int64* sidre_data_ptr =
    ds.getRoot()->getView("fields/c/int10")->getData();
  for(int j = 0; j < ndata; j++)
  {
    EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
  }
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit_external)
{
  conduit::Node input;

  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  DataStore ds;

  //Zero copy of array data
  EXPECT_TRUE(ds.getRoot()->importConduitTreeExternal(input));

  EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(), 100);
  EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(),
              3000.00,
              1e-12);
  EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(),
            std::string("foo"));

  //Scalar and string Views are never external.
  EXPECT_FALSE(ds.getRoot()->getView("fields/a/i0")->isExternal());
  EXPECT_FALSE(ds.getRoot()->getView("fields/a/d0")->isExternal());
  EXPECT_FALSE(ds.getRoot()->getView("fields/b/s0")->isExternal());

  //The View holding an array is external.
  EXPECT_TRUE(ds.getRoot()->getView("fields/c/int10")->isExternal());

  conduit::int64* sidre_data_ptr =
    ds.getRoot()->getView("fields/c/int10")->getData();
  for(int j = 0; j < ndata; j++)
  {
    EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
  }

  //Change a value in the original conduit array, then test that it's changed
  //in the Sidre external view.
  if(ndata > 3)
  {
    iarray[3] += 10;
    EXPECT_EQ(iarray[3], sidre_data_ptr[3]);
  }

  //The pointers should be the same addresses as the import treated the
  //array as an external pointer.
  EXPECT_EQ((void*)sidre_data_ptr, iarray.data_ptr());
}

//------------------------------------------------------------------------------
TEST(sidre_group, import_conduit_lists)
{
  conduit::Node input;

  //Give input some child nodes that are not lists, some that are.
  input["fields/a/i0"] = (conduit::int64)(100);
  input["fields/a/d0"] = (conduit::float64)(3000.00);
  input["fields/b/s0"] = "foo";

  int ndata = 10;
  std::vector<conduit::int64> ivec(ndata);
  input["fields/c/int10"].set(&ivec[0], ndata);
  conduit::int64_array iarray = input["fields/c/int10"].as_int64_array();
  for(int i = 0; i < ndata; ++i)
  {
    iarray[i] = (conduit::int64)i;
  }

  conduit::Node& list_item0 = input["list"].append();
  list_item0.set((conduit::int64)(12));
  conduit::Node& list_item1 = input["list"].append();
  list_item1.set((conduit::float64)(75.75));
  conduit::Node& list_item2 = input["list"].append();
  list_item2.set("test_str");
  conduit::Node& list_item3 = input["list"].append();
  list_item3["val1"] = (conduit::int64)(2);
  list_item3["val2"] = (conduit::float64)(4.0);

  DataStore ds;

  EXPECT_TRUE(ds.getRoot()->importConduitTree(input));

  {
    EXPECT_EQ(ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(),
              100);
    EXPECT_NEAR(ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(),
                3000.00,
                1e-12);
    EXPECT_EQ(ds.getRoot()->getView("fields/b/s0")->getString(),
              std::string("foo"));

    conduit::int64* sidre_data_ptr =
      ds.getRoot()->getView("fields/c/int10")->getData();

    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
    }

    Group* list = ds.getRoot()->getGroup("list");
    for(IndexType idx = list->getFirstValidGroupIndex(); indexIsValid(idx);
        idx = list->getNextValidGroupIndex(idx))
    {
      Group* child = list->getGroup(idx);
      EXPECT_EQ(child->getView("val1")->getData<conduit::int64>(), 2);
      EXPECT_NEAR(child->getView("val2")->getData<conduit::float64>(), 4.0, 1e-12);
    }

    for(IndexType idx = list->getFirstValidViewIndex(); indexIsValid(idx);
        idx = list->getNextValidViewIndex(idx))
    {
      View* view = list->getView(idx);
      const conduit::Schema& schema = view->getSchema();
      if(schema.dtype().is_int64())
      {
        EXPECT_EQ(view->getData<conduit::int64>(), 12);
      }
      if(schema.dtype().is_float64())
      {
        EXPECT_NEAR(view->getData<conduit::float64>(), 75.75, 1e-12);
      }
      if(schema.dtype().is_string())
      {
        EXPECT_EQ(view->getString(), std::string("test_str"));
      }
    }
  }

#if defined(AXOM_USE_HDF5)
  std::string lists_file = "lists.hdf5";
  std::string lists_protocol = "sidre_hdf5";
#else
  std::string lists_file = "lists.json";
  std::string lists_protocol = "sidre_json";
#endif

  ds.getRoot()->save(lists_file, lists_protocol);

  DataStore load_ds;
  load_ds.getRoot()->load(lists_file, lists_protocol);

  {
    EXPECT_EQ(
      load_ds.getRoot()->getView("fields/a/i0")->getData<conduit::int64>(),
      100);
    EXPECT_NEAR(
      load_ds.getRoot()->getView("fields/a/d0")->getData<conduit::float64>(),
      3000.00,
      1e-12);
    EXPECT_EQ(load_ds.getRoot()->getView("fields/b/s0")->getString(),
              std::string("foo"));

    conduit::int64* sidre_data_ptr =
      load_ds.getRoot()->getView("fields/c/int10")->getData();

    for(int j = 0; j < ndata; j++)
    {
      EXPECT_EQ(iarray[j], sidre_data_ptr[j]);
    }

    Group* list = load_ds.getRoot()->getGroup("list");
    for(IndexType idx = list->getFirstValidGroupIndex(); indexIsValid(idx);
        idx = list->getNextValidGroupIndex(idx))
    {
      Group* child = list->getGroup(idx);
      EXPECT_EQ(child->getView("val1")->getData<conduit::int64>(), 2);
      EXPECT_NEAR(child->getView("val2")->getData<conduit::float64>(), 4.0, 1e-12);
    }

    for(IndexType idx = list->getFirstValidViewIndex(); indexIsValid(idx);
        idx = list->getNextValidViewIndex(idx))
    {
      View* view = list->getView(idx);
      const conduit::Schema& schema = view->getSchema();
      if(schema.dtype().is_int64())
      {
        EXPECT_EQ(view->getData<conduit::int64>(), 12);
      }
      if(schema.dtype().is_float64())
      {
        EXPECT_NEAR(view->getData<conduit::float64>(), 75.75, 1e-12);
      }
      if(schema.dtype().is_string())
      {
        EXPECT_EQ(view->getString(), std::string("test_str"));
      }
    }
  }
}

#ifdef AXOM_USE_UMPIRE

class UmpireTest : public ::testing::TestWithParam<int>
{
public:
  void SetUp() override
  {
    allocID = GetParam();
    root = ds.getRoot();
  }

  void TearDown() override
  {
    axom::setDefaultAllocator(umpire::resource::Host);
  }

  static constexpr int SIZE = 100;
  DataStore ds;
  Group* root;
  umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
  int allocID;
};

TEST_P(UmpireTest, root_default_allocator)
{
  axom::setDefaultAllocator(allocID);

  DataStore dsPrime;
  ASSERT_EQ(dsPrime.getRoot()->getDefaultAllocatorID(), allocID);
}

TEST_P(UmpireTest, get_set_allocator)
{
  int defaultAllocatorID = axom::getDefaultAllocatorID();
  ASSERT_EQ(root->getDefaultAllocator().getId(), defaultAllocatorID);
  ASSERT_EQ(root->getDefaultAllocatorID(), defaultAllocatorID);

  root->setDefaultAllocator(allocID);
  defaultAllocatorID = axom::getDefaultAllocatorID();
  ASSERT_EQ(root->getDefaultAllocatorID(), allocID);

  root->setDefaultAllocator(rm.getInstance().getAllocator(defaultAllocatorID));
  ASSERT_EQ(root->getDefaultAllocator().getId(), defaultAllocatorID);
  ASSERT_EQ(root->getDefaultAllocatorID(), defaultAllocatorID);
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate)
{
  {
    View* view = root->createViewAndAllocate("v", INT_ID, SIZE, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    IndexType shape[] = {1, SIZE, 1};
    View* view = root->createViewAndAllocate("v", INT_ID, 3, shape, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    View* view = root->createViewAndAllocate("v", dtype, allocID);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

//------------------------------------------------------------------------------
TEST_P(UmpireTest, allocate_default)
{
  root->setDefaultAllocator(allocID);

  {
    View* view = root->createViewAndAllocate("v", INT_ID, SIZE);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    IndexType shape[] = {1, SIZE, 1};
    View* view = root->createViewAndAllocate("v", INT_ID, 3, shape);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }

  {
    DataType dtype = conduit::DataType::default_dtype(INT_ID);
    dtype.set_number_of_elements(SIZE);
    View* view = root->createViewAndAllocate("v", dtype);

    ASSERT_EQ(allocID, rm.getAllocator(view->getVoidPtr()).getId());
    root->destroyViewAndData("v");
  }
}

const int allocators[] = {
  axom::getUmpireResourceAllocatorID(umpire::resource::Host)
  #ifdef AXOM_USE_CUDA

    #ifdef UMPIRE_ENABLE_PINNED
    ,
  axom::getUmpireResourceAllocatorID(umpire::resource::Pinned)
    #endif

    #ifdef UMPIRE_ENABLE_DEVICE
    ,
  axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    #endif

    #ifdef UMPIRE_ENABLE_CONST
    ,
  axom::getUmpireResourceAllocatorID(umpire::resource::Constant)
    #endif

    #ifdef UMPIRE_ENABLE_UM
    ,
  axom::getUmpireResourceAllocatorID(umpire::resource::Unified)
    #endif

  #endif /* AXOM_USE_CUDA */
};

INSTANTIATE_TEST_SUITE_P(sidre_group, UmpireTest, ::testing::ValuesIn(allocators));

#endif  // AXOM_USE_UMPIRE
