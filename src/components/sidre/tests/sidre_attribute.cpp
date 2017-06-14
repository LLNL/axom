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
using axom::sidre::Node;

// Create some attribute values
const std::string name_color("color");
const std::string color_none("white");
const std::string color_red("red");
const std::string color_blue("blue");

const std::string name_animal("animal");
const std::string animal_none("human");
const std::string animal_cat("cat");
const std::string animal_dog("dog");

const std::string namea("a");
const std::string nameb("b");

const std::string name_dump("dump");
const int dump_no = 0;
const int dump_yes = 1;

const std::string name_size("size");
const double size_small = 1.2;
const double size_medium = 2.3;
const double size_large = 3.4;

//------------------------------------------------------------------------------
// Create attribute in a Datastore
//
TEST(sidre_attribute,create_attr)
{
  const std::string dump_none("none");
  bool ok;

  DataStore * ds = new DataStore();

  int nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  bool has_index = ds->hasAttribute(0);
  EXPECT_FALSE( has_index );
  bool has_name = ds->hasAttribute(name_color);
  EXPECT_FALSE( has_name );

  // Create string attribute
  Attribute * color = ds->createAttribute(name_color, color_none);
  EXPECT_TRUE( color != AXOM_NULLPTR );

  IndexType attr_index = color->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute(name_color);
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

  // Try to change default to a different type.
  // Check template of setDefault.
  ok = color->setDefault(1);
  EXPECT_FALSE(ok);
  ok = color->setDefault(3.14);
  EXPECT_FALSE(ok);

  // Change to legal values.
  ok = color->setDefault("unknown");
  EXPECT_TRUE(ok);
#if 0
  // XXX - this causes the default type to change to int.
  ok = color->setDefault('u');
  EXPECT_TRUE(ok);
#endif
  ok = color->setDefault(std::string("string"));  // non-const string
  EXPECT_TRUE(ok);

  Attribute * attr = ds->getAttribute(name_color);
  EXPECT_EQ(attr, color);
  const Attribute * attrc = ds->getAttribute(name_color);
  EXPECT_EQ(attrc, color);

  attr = ds->getAttribute(0);
  EXPECT_EQ(attr, color);
  attrc = ds->getAttribute(0);
  EXPECT_EQ(attrc, color);

  ds->destroyAttribute(color);
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);
  has_name = ds->hasAttribute(name_color);
  EXPECT_FALSE( has_name );
  // At this point color points to deallocated memory

  // Create additional attributes
  Attribute * dump = ds->createAttribute(name_dump, dump_no);
  EXPECT_TRUE( dump != AXOM_NULLPTR );

  attr_index = dump->getIndex();
  EXPECT_EQ(0, attr_index);

  Attribute * size = ds->createAttribute(name_size, size_small);
  EXPECT_TRUE( dump != AXOM_NULLPTR );

  attr_index = size->getIndex();
  EXPECT_EQ(1, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(2, nattrs);

  ok = dump->setDefault(1);
  EXPECT_TRUE(ok);
  // try to change default to a different type
  ok = dump->setDefault(name_dump);
  EXPECT_FALSE(ok);
  ok = dump->setDefault("yes");
  EXPECT_FALSE(ok);
  ok = dump->setDefault(3.1415);
  EXPECT_FALSE(ok);

  ds->destroyAllAttributes();
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  delete ds;
}

//------------------------------------------------------------------------------
// Set attributes on a view

TEST(sidre_attribute,view_attr)
{
  bool ok;

  DataStore * ds = new DataStore();

  // Create all attributes for DataStore
  Attribute * attr_color = ds->createAttribute(name_color, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );

  Attribute * attr_animal = ds->createAttribute(name_animal, animal_none);
  EXPECT_TRUE( attr_animal != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  EXPECT_FALSE(view1a->hasAttributeValue(AXOM_NULLPTR));
  EXPECT_FALSE(view1a->hasAttributeValue(attr_color));

  // Check values of unset attributes
  const std::string out1x = view1a->getAttributeValueString(attr_color);
  EXPECT_EQ(color_none, out1x);

  const std::string out1y = view1a->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_none, out1y);

  ok = view1a->setAttributeValue(attr_color, color_red);
  EXPECT_TRUE( ok );

  EXPECT_TRUE(view1a->hasAttributeValue(attr_color));

  const std::string out = view1a->getAttributeValueString(attr_color);
  EXPECT_EQ(color_red, out);

  // reset attribute value
  ok = view1a->setAttributeValue(attr_color, color_blue);
  EXPECT_TRUE( ok );

  const std::string out1b = view1a->getAttributeValueString(attr_color);
  EXPECT_EQ(color_blue, out1b);

  // Check second, unset attribute. Should be default value
  EXPECT_FALSE(view1a->hasAttributeValue(attr_animal));
  const std::string out1d = view1a->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_none, out1d);

  // Now set second attribute
  ok = view1a->setAttributeValue(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  const std::string out1c = view1a->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_dog, out1c);

  //----------------------------------------
  // Set the second attribute in a Group
  Group * grp2 = root->createGroup("grp2");

  View  * view2a = grp2->createView(namea);
  EXPECT_TRUE( view2a != AXOM_NULLPTR );

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_FALSE(view2a->hasAttributeValue(attr_animal));

  ok = view2a->setAttributeValue(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2a = view2a->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_dog, out2a);

  // Get the first, unset, attribute
  const std::string out2b = view2a->getAttributeValueString(attr_color);
  EXPECT_EQ(color_none, out2b);

  // Now set first attribute
  ok = view2a->setAttributeValue(attr_color, color_red);
  EXPECT_TRUE( ok );

  EXPECT_TRUE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2c = view2a->getAttributeValueString(attr_color);
  EXPECT_EQ(color_red, out2c);


  //----------------------------------------
  // Set attribute on second View in a Group
  Group * grp3 = root->createGroup("grp3");
  View  * view3a = grp3->createView(namea);
  EXPECT_TRUE( view3a != AXOM_NULLPTR );
  View  * view3b = grp3->createView(nameb);
  EXPECT_TRUE( view3b != AXOM_NULLPTR );

  ok = view3b->setAttributeValue(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  EXPECT_FALSE(view3b->hasAttributeValue(attr_color));
  EXPECT_TRUE(view3b->hasAttributeValue(attr_animal));

  const std::string & out3a = view3b->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_dog, out3a);

  //----------------------------------------
  // Moving a view should preserve attributes
  Group * grp4 = root->createGroup("grp4");

  grp4->moveView(view3b);

  const std::string & out4a = view3b->getAttributeValueString(attr_animal);
  EXPECT_EQ(animal_dog, out4a);

  // Create an attribute which will be destroyed
  view3a->setAttributeValue(attr_animal, animal_dog);

  grp3->destroyView(namea);
  grp4->destroyView(nameb);

  delete ds;
}

//------------------------------------------------------------------------------
// Use different type of attributes

TEST(sidre_attribute,view_int_and_double)
{
  bool ok;

  DataStore * ds = new DataStore();

  // Create all attributes for DataStore
  Attribute * attr_dump = ds->createAttribute(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );

  Attribute * attr_size = ds->createAttribute(name_size, size_small);
  EXPECT_TRUE( attr_size != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  //----------------------------------------
  // Create a View
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  // Get default values
  int dump = view1a->getAttributeValue(attr_dump);
  EXPECT_EQ( dump_no, dump );

  double size = view1a->getAttributeValue(attr_size);
  EXPECT_EQ( size_small, size );

  // Set values
  ok = view1a->setAttributeValue(attr_dump, dump_yes);
  EXPECT_TRUE( ok );
  dump = -1; // clear value
  dump = view1a->getAttributeValue(attr_dump);
  EXPECT_EQ( dump_yes, dump );

  ok = view1a->setAttributeValue(attr_size, size_medium);
  EXPECT_TRUE( ok );
  size = 0.0;  // clear value
  size = view1a->getAttributeValue(attr_size);
  EXPECT_EQ( size_medium, size );

  // Set values with incorrect types
  ok = view1a->setAttributeValue(attr_dump, size_small);
  EXPECT_FALSE( ok );
  ok = view1a->setAttributeValue(attr_dump, namea);
  EXPECT_FALSE( ok );
#if 0
  ok = view1a->setAttributeValue(attr_dump, 'a');
  EXPECT_FALSE( ok );
#endif
  ok = view1a->setAttributeValue(attr_dump, "namea");
  EXPECT_FALSE( ok );

  // Error checks
#if 0
  // AttrValue::getValueNodeRef will fail an assert
  int err = -1;
  err = view1a->getAttributeValue(AXOM_NULLPTR);
  EXPECT_EQ(-1, err);
#endif

  delete ds;
}

//------------------------------------------------------------------------------
// Access attributes by name or index

#if 0
TEST(sidre_attribute,overloads)
{
  IndexType icolor = color->getIndex();

  Group * root = ds.getRoot();
  View * view = root->createView("var1");

  view->setAttributeValue(color, color_red);
  view->setAttributeValue(icolor, color_red);
  view->setAttributeValue("color", color_red);

  const char * attr1 = view->getAttributeValueString(color);
  const char * attr2 = view->getAttributeValueString(icolor);
  const char * attr3 = view->getAttributeValueString("color");


  view = root->createView("var2");
  // Get attributes without setting returns default value
  const char * attr1 = view->getAttributeValueString(color);

  delete ds;
}
#endif

//------------------------------------------------------------------------------

TEST(sidre_attribute,as_node)
{
  bool ok;

  DataStore * ds = new DataStore();

  // Create attributes for DataStore
  Attribute * attr_color = ds->createAttribute(name_color, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );

  Attribute * attr_dump = ds->createAttribute(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  ok = view1a->setAttributeValue(attr_color, color_red);
  EXPECT_TRUE( ok );

  const Node & node = view1a->getAttributeValueNodeRef(attr_color);
  EXPECT_EQ(color_red, node.as_string());

  const Node & node2 = view1a->getAttributeValueNodeRef(attr_dump);
  EXPECT_EQ(dump_no, node2.as_int32());

  delete ds;
}

//------------------------------------------------------------------------------
