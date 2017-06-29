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
using axom::sidre::DOUBLE_ID;
using axom::sidre::INT_ID;
using axom::sidre::CHAR8_STR_ID;

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
  Attribute * color = ds->createAttributeString(name_color, color_none);
  EXPECT_TRUE( color != AXOM_NULLPTR );
  EXPECT_EQ( CHAR8_STR_ID, color->getTypeID());

  IndexType attr_index = color->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute(name_color);
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

  // Try to change default to a different type.
  // Check template of setDefaultScalar.
  ok = color->setDefaultScalar(1);
  EXPECT_FALSE(ok);
  ok = color->setDefaultScalar(3.14);
  EXPECT_FALSE(ok);

  // Change to legal values.
  ok = color->setDefaultString("unknown");
  EXPECT_TRUE(ok);
#if 0
  // XXX - unsupported overload
  ok = color->setDefaultString('u');
  EXPECT_TRUE(ok);
#endif
  ok = color->setDefaultString(std::string("string"));  // non-const string
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
  Attribute * dump = ds->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( dump != AXOM_NULLPTR );

  attr_index = dump->getIndex();
  EXPECT_EQ(0, attr_index);

  Attribute * size = ds->createAttributeScalar(name_size, size_small);
  EXPECT_TRUE( dump != AXOM_NULLPTR );

  attr_index = size->getIndex();
  EXPECT_EQ(1, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(2, nattrs);

  ok = dump->setDefaultScalar(1);
  EXPECT_TRUE(ok);
  // try to change default to a different type
  ok = dump->setDefaultString(name_dump);
  EXPECT_FALSE(ok);
  ok = dump->setDefaultString("yes");
  EXPECT_FALSE(ok);
  ok = dump->setDefaultScalar(3.1415);
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
  Attribute * attr_color = ds->createAttributeString(name_color, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );

  Attribute * attr_animal = ds->createAttributeString(name_animal, animal_none);
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
  const std::string out1x = view1a->getAttributeString(attr_color);
  EXPECT_EQ(color_none, out1x);

  const std::string out1y = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(animal_none, out1y);

  ok = view1a->setAttributeString(attr_color, color_red);
  EXPECT_TRUE( ok );

  EXPECT_TRUE(view1a->hasAttributeValue(attr_color));

  const std::string out = view1a->getAttributeString(attr_color);
  EXPECT_EQ(color_red, out);

  // reset attribute value
  ok = view1a->setAttributeString(attr_color, color_blue);
  EXPECT_TRUE( ok );

  const std::string out1b = view1a->getAttributeString(attr_color);
  EXPECT_EQ(color_blue, out1b);

  // Check second, unset attribute. Should be default value
  EXPECT_FALSE(view1a->hasAttributeValue(attr_animal));
  const std::string out1d = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(animal_none, out1d);

  // Now set second attribute
  ok = view1a->setAttributeString(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  const std::string out1c = view1a->getAttributeString(attr_animal);
  EXPECT_EQ(animal_dog, out1c);

  //----------------------------------------
  // Set the second attribute in a Group
  Group * grp2 = root->createGroup("grp2");

  View  * view2a = grp2->createView(namea);
  EXPECT_TRUE( view2a != AXOM_NULLPTR );

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_FALSE(view2a->hasAttributeValue(attr_animal));

  ok = view2a->setAttributeString(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  EXPECT_FALSE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2a = view2a->getAttributeString(attr_animal);
  EXPECT_EQ(animal_dog, out2a);

  // Get the first, unset, attribute
  const std::string out2b = view2a->getAttributeString(attr_color);
  EXPECT_EQ(color_none, out2b);

  // Now set first attribute
  ok = view2a->setAttributeString(attr_color, color_red);
  EXPECT_TRUE( ok );

  EXPECT_TRUE(view2a->hasAttributeValue(attr_color));
  EXPECT_TRUE(view2a->hasAttributeValue(attr_animal));

  const std::string out2c = view2a->getAttributeString(attr_color);
  EXPECT_EQ(color_red, out2c);

  // Try to get a scalar from string
  int novalue = view2a->getAttributeScalar(attr_color);
  EXPECT_EQ(0, novalue);

  //----------------------------------------
  // Set attribute on second View in a Group
  Group * grp3 = root->createGroup("grp3");
  View  * view3a = grp3->createView(namea);
  EXPECT_TRUE( view3a != AXOM_NULLPTR );
  View  * view3b = grp3->createView(nameb);
  EXPECT_TRUE( view3b != AXOM_NULLPTR );

  ok = view3b->setAttributeString(attr_animal, animal_dog);
  EXPECT_TRUE( ok );

  EXPECT_FALSE(view3b->hasAttributeValue(attr_color));
  EXPECT_TRUE(view3b->hasAttributeValue(attr_animal));

  const std::string & out3a = view3b->getAttributeString(attr_animal);
  EXPECT_EQ(animal_dog, out3a);

  //----------------------------------------
  // Moving a view should preserve attributes
  Group * grp4 = root->createGroup("grp4");

  grp4->moveView(view3b);

  const std::string & out4a = view3b->getAttributeString(attr_animal);
  EXPECT_EQ(animal_dog, out4a);

  // Create an attribute which will be destroyed
  view3a->setAttributeString(attr_animal, animal_dog);

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
  Attribute * attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );
  EXPECT_EQ( INT_ID, attr_dump->getTypeID());

  Attribute * attr_size = ds->createAttributeScalar(name_size, size_small);
  EXPECT_TRUE( attr_size != AXOM_NULLPTR );
  EXPECT_EQ( DOUBLE_ID, attr_size->getTypeID());

  Group * root = ds->getRoot();

  //----------------------------------------
  // Create a View
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  // Get default values
  int dump = view1a->getAttributeScalar(attr_dump);
  EXPECT_EQ( dump_no, dump );

  double size = view1a->getAttributeScalar(attr_size);
  EXPECT_EQ( size_small, size );

  // Set values
  ok = view1a->setAttributeScalar(attr_dump, dump_yes);
  EXPECT_TRUE( ok );
  dump = -1; // clear value
  dump = view1a->getAttributeScalar(attr_dump);
  EXPECT_EQ( dump_yes, dump );

  ok = view1a->setAttributeScalar(attr_size, size_medium);
  EXPECT_TRUE( ok );
  size = 0.0;  // clear value
  size = view1a->getAttributeScalar(attr_size);
  EXPECT_EQ( size_medium, size );

  // Set values with incorrect types
  ok = view1a->setAttributeScalar(attr_dump, size_small);
  EXPECT_FALSE( ok );
  ok = view1a->setAttributeString(attr_dump, namea);
  EXPECT_FALSE( ok );
#if 0
  ok = view1a->setAttributeString(attr_dump, 'a');
  EXPECT_FALSE( ok );
#endif
  ok = view1a->setAttributeString(attr_dump, "namea");
  EXPECT_FALSE( ok );

  // Try to get a string from a scalar
  const char * nostr = view1a->getAttributeString(attr_dump);
  EXPECT_EQ(AXOM_NULLPTR, nostr);

  int i = -1;
  i = view1a->getAttributeScalar(AXOM_NULLPTR);
  EXPECT_EQ(0, i);

  delete ds;
}

//------------------------------------------------------------------------------
// Reset attribute to default

TEST(sidre_attribute,set_default)
{
  bool ok;

  DataStore * ds = new DataStore();

  // Create all attributes for DataStore
  Attribute * attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );
  EXPECT_EQ( INT_ID, attr_dump->getTypeID());

  Attribute * attr_size = ds->createAttributeScalar(name_size, size_small);
  EXPECT_TRUE( attr_size != AXOM_NULLPTR );
  EXPECT_EQ( DOUBLE_ID, attr_size->getTypeID());

  Group * root = ds->getRoot();

  //----------------------------------------
  // Create a View
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  // reset unset attribute 1
  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  ok = view1a->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  // Set value
  ok = view1a->setAttributeScalar(attr_dump, dump_yes);
  EXPECT_TRUE( ok );
  EXPECT_TRUE(view1a->hasAttributeValue(attr_dump));

  // reset set attribute 1
  ok = view1a->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);
  EXPECT_FALSE(view1a->hasAttributeValue(attr_dump));

  // reset unset attribute 2
  EXPECT_FALSE(view1a->hasAttributeValue(attr_size));

  ok = view1a->setAttributeToDefault(attr_size);
  EXPECT_TRUE(ok);

  EXPECT_FALSE(view1a->hasAttributeValue(attr_size));

  // Check errors
  ok = view1a->setAttributeToDefault(AXOM_NULLPTR);
  EXPECT_FALSE( ok );

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

  view->setAttributeString(color, color_red);
  view->setAttributeString(icolor, color_red);
  view->setAttributeString("color", color_red);

  const char * attr1 = view->getAttributeString(color);
  const char * attr2 = view->getAttributeString(icolor);
  const char * attr3 = view->getAttributeString("color");


  view = root->createView("var2");
  // Get attributes without setting returns default value
  const char * attr1 = view->getAttributeString(color);

  delete ds;
}
#endif

//------------------------------------------------------------------------------

TEST(sidre_attribute,as_node)
{
  bool ok;

  DataStore * ds = new DataStore();

  // Create attributes for DataStore
  Attribute * attr_color = ds->createAttributeString(name_color, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );

  Attribute * attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  //----------------------------------------
  // Set the first attribute in a Group
  Group * grp1 = root->createGroup("grp1");
  View  * view1a = grp1->createView(namea);
  EXPECT_TRUE( view1a != AXOM_NULLPTR );

  ok = view1a->setAttributeString(attr_color, color_red);
  EXPECT_TRUE( ok );

  const Node & node1 = view1a->getAttributeNodeRef(attr_color);
  EXPECT_EQ(color_red, node1.as_string());

  const Node & node2 = view1a->getAttributeNodeRef(attr_dump);
  EXPECT_EQ(dump_no, node2.as_int32());

  const Node & node3 = view1a->getAttributeNodeRef(AXOM_NULLPTR);
  EXPECT_TRUE(node3.schema().dtype().is_empty());

  delete ds;
}

//------------------------------------------------------------------------------
