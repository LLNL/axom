// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/sidre/core/sidre.hpp"

using axom::sidre::Attribute;
using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::Iterator;
using axom::sidre::View;

// Create some attribute values
const std::string name_color("color");
const std::string color_none("white");
//const std::string color_red("red");
//const std::string color_blue("blue");

const std::string name_animal("animal");
const std::string animal_none("human");
//const std::string animal_cat("cat");
//const std::string animal_dog("dog");

const std::string name_dump("dump");
const int dump_no = 0;
//const int dump_yes = 1;

//------------------------------------------------------------------------------
// Internal routine to create a datastore to be used with various iteration
// schemes
// Error checks here are minimal since it is assumed all of the routines used
// have already been tested.

DataStore* sample_datastore(void)
{
  DataStore* ds = new DataStore();

  // Create all attributes for DataStore
  Attribute* attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  Attribute* attr_color = ds->createAttributeString(name_color, color_none);
  Attribute* attr_animal = ds->createAttributeString(name_animal, animal_none);

  EXPECT_TRUE(attr_dump != nullptr);
  EXPECT_TRUE(attr_color != nullptr);
  EXPECT_TRUE(attr_animal != nullptr);

  Group* root = ds->getRoot();

  Group* grpA = root->createGroup("grpA");
  Group* grpB = root->createGroup("grpB");
  Group* grpBB = grpB->createGroup("grpBB");
  root->createGroup("grpC");  // No Views

  // Tree of empty views
  Group* grpD0 = root->createGroup("grpD0");
  Group* grpD1 = grpD0->createGroup("grpD1");
  Group* grpD2 = grpD1->createGroup("grpD2");
  Group* grpD3 = grpD1->createGroup("grpD3");
  Group* grpD4 = grpD0->createGroup("grpD4");
  Group* grpD5 = grpD4->createGroup("grpD5");

  // checks to silence compiler about unused variables
  EXPECT_TRUE(grpD2 != nullptr);
  EXPECT_TRUE(grpD3 != nullptr);
  EXPECT_TRUE(grpD5 != nullptr);

  grpA->createViewScalar("grpA_view1", 1);
  grpA->createViewScalar("grpA_view2", 2);

  grpBB->createViewScalar("grpBB_view3", 3);

  grpB->createViewScalar("grpB_view4", 4);

  root->createViewScalar("root_view5", 5);
  root->createViewScalar("root_view6", 6);
  root->createViewScalar("root_view7", 7);

  return ds;
}

//------------------------------------------------------------------------------
// Iterate thru datastore with a cursor

TEST(sidre_iterator, depth_first)
{
  enum nodeclass
  {
    GROUP,
    VIEW
  };

  struct reference
  {
    const char* name;
    nodeclass cls;
  };

  DataStore* ds = sample_datastore();

  Group* root = ds->getRoot();

  Iterator qitr(root);

  int iorder = 0;
  reference order[] = {
    {"grpA_view1", VIEW},
    {"grpA_view2", VIEW},
    {"grpA", GROUP},
    {"grpBB_view3", VIEW},
    {"grpBB", GROUP},
    {"grpB_view4", VIEW},
    {"grpB", GROUP},
    {"grpC", GROUP},
    {"grpD2", GROUP},
    {"grpD3", GROUP},
    {"grpD1", GROUP},
    {"grpD5", GROUP},
    {"grpD4", GROUP},
    {"grpD0", GROUP},
    {"root_view5", VIEW},
    {"root_view6", VIEW},
    {"root_view7", VIEW},
    {"", GROUP},
  };

  while(qitr.isValid())
  {
    // find our current path
    const std::string& name = qitr.getName();
    //std::cout << name << std::endl;

    EXPECT_EQ(order[iorder].name, name);

    if(order[iorder].cls == GROUP)
    {
      EXPECT_TRUE(qitr.isGroup());
      EXPECT_FALSE(qitr.isView());

      Group* grp_out = qitr.asGroup();
      EXPECT_EQ(order[iorder].name, grp_out->getName());

      Group const* grp_const = qitr.asGroup();
      EXPECT_EQ(order[iorder].name, grp_const->getName());

      View* view_out = qitr.asView();
      EXPECT_EQ(view_out, static_cast<void*>(nullptr));

      View const* view_const = qitr.asView();
      EXPECT_EQ(view_const, static_cast<void*>(nullptr));
    }
    else
    {
      EXPECT_FALSE(qitr.isGroup());
      EXPECT_TRUE(qitr.isView());

      Group* grp_out = qitr.asGroup();
      EXPECT_EQ(grp_out, static_cast<void*>(nullptr));

      Group const* grp_const = qitr.asGroup();
      EXPECT_EQ(grp_const, static_cast<void*>(nullptr));

      View* view_out = qitr.asView();
      EXPECT_EQ(order[iorder].name, view_out->getName());

      View const* view_const = qitr.asView();
      EXPECT_EQ(order[iorder].name, view_const->getName());
    }

    qitr.advanceToNext();
    iorder++;
  }

  // Check invalid iterator
  EXPECT_FALSE(qitr.isGroup());
  EXPECT_FALSE(qitr.isView());

  delete ds;
}
