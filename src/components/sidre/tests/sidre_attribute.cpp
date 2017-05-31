/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "sidre/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::Attribute;
using axom::sidre::Group;
using axom::sidre::View;
using axom::sidre::IndexType;

//------------------------------------------------------------------------------
TEST(sidre_attribute,create_attr)
{
  const std::string color_none("none");
  const std::string dump_none("none");

  DataStore * ds = new DataStore();

  int nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  bool has_index = ds->hasAttribute(0);
  EXPECT_FALSE( has_index );
  std::string nameattr("color");
  bool has_name = ds->hasAttribute(nameattr);
  EXPECT_FALSE( has_name );

  Attribute * color = ds->createAttribute(nameattr, color_none);
  EXPECT_TRUE( color != AXOM_NULLPTR );

  IndexType attr_index = color->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute(nameattr);
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

  Attribute * attr = ds->getAttribute(nameattr);
  EXPECT_EQ(attr, color);
  const Attribute * attrc = ds->getAttribute(nameattr);
  EXPECT_EQ(attrc, color);

  attr = ds->getAttribute(0);
  EXPECT_EQ(attr, color);
  attrc = ds->getAttribute(0);
  EXPECT_EQ(attrc, color);

  ds->destroyAttribute(color);
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);
  has_name = ds->hasAttribute(nameattr);
  EXPECT_FALSE( has_name );
  // At this point color points to deallocated memory

  // Create additional attributes
  std::string namedump1("dump1");
  Attribute * dump1 = ds->createAttribute(namedump1, dump_none);
  EXPECT_TRUE( dump1 != AXOM_NULLPTR );

  attr_index = dump1->getIndex();
  EXPECT_EQ(0, attr_index);

  std::string namedump2("dump2");
  Attribute * dump2 = ds->createAttribute(namedump2, dump_none);
  EXPECT_TRUE( dump1 != AXOM_NULLPTR );

  attr_index = dump2->getIndex();
  EXPECT_EQ(1, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(2, nattrs);

  ds->destroyAllAttributes();
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

}

//------------------------------------------------------------------------------

TEST(sidre_attribute,view_attr)
{
  // Create some attribute values
  const std::string color_none("white");
  const std::string color_red("red");
  const std::string color_blue("blue");

  const std::string animal_none("human");
  const std::string animal_cat("cat");
  const std::string animal_dog("dog");

  bool ok;

  DataStore * ds = new DataStore();

  // Create all attributes for DataStore
  std::string nameattr1("color");
  Attribute * attr_color = ds->createAttribute(nameattr1, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );

  std::string nameattr2("animal");
  Attribute * attr_animal = ds->createAttribute(nameattr2, animal_none);
  EXPECT_TRUE( attr_animal != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView("a");
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  // Check values of unset attributes
  const std::string & out1x = view1a->getAttribute(attr_color);
  EXPECT_EQ(color_none, out1x);

  const std::string & out1y = view1a->getAttribute(attr_animal);
  EXPECT_EQ(animal_none, out1y);

  ok = view1a->setAttribute(attr_color, color_red);
  EXPECT_TRUE( ok );

  const std::string & out = view1a->getAttribute(attr_color);
  EXPECT_EQ(color_red, out);

  // reset attribute value
  ok = view1a->setAttribute(attr_color, color_blue);
  EXPECT_TRUE( ok );

  const std::string & out1b = view1a->getAttribute(attr_color);
  EXPECT_EQ(color_blue, out1b);

  // Now set second attribute
  ok = view1a->setAttribute(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  const std::string & out1c = view1a->getAttribute(attr_animal);
  EXPECT_EQ(animal_dog, out1c);

  //----------------------------------------
  // Set the second attribute in a Group
  Group * grp2 = root->createGroup("grp2");

  View  * view2a = grp2->createView("a");
  EXPECT_TRUE( view2a != AXOM_NULLPTR );

  ok = view2a->setAttribute(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  const std::string & out2a = view2a->getAttribute(attr_animal);
  EXPECT_EQ(animal_dog, out2a);

  // Now set first attribute
  ok = view2a->setAttribute(attr_color, color_red);
  EXPECT_TRUE( ok );

  const std::string & out2b = view2a->getAttribute(attr_color);
  EXPECT_EQ(color_red, out2b);


  //----------------------------------------
  // Set attribute on second View in a Group
  Group * grp3 = root->createGroup("grp3");
  View  * view3a = grp3->createView("a");
  EXPECT_TRUE( view3a != AXOM_NULLPTR );
  View  * view3b = grp3->createView("b");
  EXPECT_TRUE( view3b != AXOM_NULLPTR );

  ok = view3b->setAttribute(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  const std::string & out3a = view3b->getAttribute(attr_animal);
  EXPECT_EQ(animal_dog, out3a);

  //----------------------------------------
  // Moving a view should preserve attributes
  Group * grp4 = root->createGroup("grp4");

  grp4->moveView(view3b);

#if 0
  const std::string & out4a = view3b->getAttribute(attr_animal);
  EXPECT_EQ(animal_dog, out4a);

  grp3->destroyView("a");
  grp4->destroyView("b");
#endif

#if 0
  IndexType icolor = color->getIndex();

  Group * root = ds.getRoot();
  View * view = root->createView("var1");

  view->setAttribute(color, "red");
  view->setAttribute(icolor, "red");
  view->setAttribute("color", "red");

  const char * attr1 = view->getAttribute(color);
  const char * attr2 = view->getAttribute(icolor);
  const char * attr3 = view->getAttribute("color");


  view = root->createView("var2");
  // Get attributes without setting returns default value
  const char * attr1 = view->getAttribute(color);
#endif
}

//------------------------------------------------------------------------------
