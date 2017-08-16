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
using axom::sidre::QueryIterator;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
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

// intel has a problem overloading 'Attribute *' and 'IndexType'.
Attribute * attr_null = AXOM_NULLPTR;

// Test protocols
int nprotocols = 3;
std::string const protocols[] = { "sidre_json", "sidre_hdf5", "json" };

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
  EXPECT_TRUE( size != AXOM_NULLPTR );

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

  EXPECT_FALSE(view1a->hasAttributeValue(attr_null));
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
  i = view1a->getAttributeScalar(attr_null);
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
  ok = view1a->setAttributeToDefault(attr_null);
  EXPECT_FALSE( ok );

  delete ds;
}

//------------------------------------------------------------------------------
// get attribute as Conduit::Node

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

  const Node & node3 = view1a->getAttributeNodeRef(attr_null);
  EXPECT_TRUE(node3.schema().dtype().is_empty());

  delete ds;
}

//------------------------------------------------------------------------------
// Access attributes by name or index

TEST(sidre_attribute,overloads)
{
   bool ok;
   DataStore * ds = new DataStore();

  // Create string and scalar attributes
  Attribute * attr_color = ds->createAttributeString(name_color, color_none);
  EXPECT_TRUE( attr_color != AXOM_NULLPTR );
  IndexType icolor = attr_color->getIndex();
  EXPECT_EQ(0, icolor);

  Attribute * attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( attr_dump != AXOM_NULLPTR );
  IndexType idump = attr_dump->getIndex();
  EXPECT_EQ(1, idump);

  EXPECT_EQ(attr_color, ds->getAttribute(name_color));
  EXPECT_EQ(attr_color, ds->getAttribute(icolor));

  //----------------------------------------
  Group * root = ds->getRoot();
  View * view = root->createView("view1");

  // string
  ok = view->setAttributeString(attr_color, color_red);
  EXPECT_TRUE(ok);
  ok = view->setAttributeString(icolor, color_red);
  EXPECT_TRUE(ok);
  ok = view->setAttributeString(name_color, color_red);
  EXPECT_TRUE(ok);

  const char * attr1a = view->getAttributeString(attr_color);
  EXPECT_EQ(color_red, attr1a);
  const char * attr2a = view->getAttributeString(icolor);
  EXPECT_EQ(color_red, attr2a);
  const char * attr3a = view->getAttributeString(name_color);
  EXPECT_EQ(color_red, attr3a);

  // scalar
  ok = view->setAttributeScalar(attr_dump, dump_yes);
  EXPECT_TRUE(ok);
  ok = view->setAttributeScalar(idump, dump_yes);
  EXPECT_TRUE(ok);
  ok = view->setAttributeScalar(name_dump, dump_yes);
  EXPECT_TRUE(ok);

  int attr1b = view->getAttributeScalar(attr_dump);
  EXPECT_EQ(dump_yes, attr1b);
  int attr2b = view->getAttributeScalar(idump);
  EXPECT_EQ(dump_yes, attr2b);
  int attr3b = view->getAttributeScalar(name_dump);
  EXPECT_EQ(dump_yes, attr3b);

  EXPECT_EQ(dump_yes, view->getAttributeScalar<int>(attr_dump));
  EXPECT_EQ(dump_yes, view->getAttributeScalar<int>(idump));
  EXPECT_EQ(dump_yes, view->getAttributeScalar<int>(name_dump));

  const Node & node1 = view->getAttributeNodeRef(attr_dump);
  EXPECT_EQ(dump_yes, node1.as_int());
  const Node & node2 = view->getAttributeNodeRef(idump);
  EXPECT_EQ(dump_yes, node2.as_int());
  const Node & node3 = view->getAttributeNodeRef(name_dump);
  EXPECT_EQ(dump_yes, node3.as_int());

  EXPECT_TRUE(view->hasAttributeValue(attr_dump));
  EXPECT_TRUE(view->hasAttributeValue(idump));
  EXPECT_TRUE(view->hasAttributeValue(name_dump));

  ok = view->setAttributeToDefault(attr_dump);
  EXPECT_TRUE(ok);
  ok = view->setAttributeToDefault(idump);
  EXPECT_TRUE(ok);
  ok = view->setAttributeToDefault(name_dump);
  EXPECT_TRUE(ok);

  // Attribute no longer set
  EXPECT_FALSE(view->hasAttributeValue(attr_dump));
  EXPECT_FALSE(view->hasAttributeValue(idump));
  EXPECT_FALSE(view->hasAttributeValue(name_dump));

  // Check some errors
  EXPECT_EQ(0, view->getAttributeScalar<int>(attr_null));
  EXPECT_EQ(0, view->getAttributeScalar<int>(InvalidIndex));
  EXPECT_EQ(0, view->getAttributeScalar<int>("noname"));

  delete ds;
}

//------------------------------------------------------------------------------
// Test looping over Attributes and Attribute Values.

TEST(sidre_attribute, loop_attributes)
{
  DataStore * ds = new DataStore();

  // Create attributes for DataStore
  Attribute * color = ds->createAttributeString(name_color, color_none);
  IndexType icolor = color->getIndex();
  EXPECT_EQ(0, icolor);

  Attribute * dump = ds->createAttributeScalar(name_dump, dump_no);
  IndexType idump = dump->getIndex();
  EXPECT_EQ(1, idump);

  Attribute * size = ds->createAttributeScalar(name_size, size_small);
  IndexType isize = size->getIndex();
  EXPECT_EQ(2, isize);

  {
    IndexType idx1 = ds->getFirstValidAttributeIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = ds->getNextValidAttributeIndex(idx1);
    EXPECT_EQ(1, idx2);
    IndexType idx3 = ds->getNextValidAttributeIndex(idx2);
    EXPECT_EQ(2, idx3);
    IndexType idx4 = ds->getNextValidAttributeIndex(idx3);
    EXPECT_EQ(InvalidIndex, idx4);
    IndexType idx5 = ds->getNextValidAttributeIndex(idx4);
    EXPECT_EQ(InvalidIndex, idx5);
  }

  //----------------------------------------
  Group * root = ds->getRoot();

  // set all attributes
  View * view1 = root->createView("view1");
  view1->setAttributeString(color, color_red);
  view1->setAttributeScalar(dump, dump_yes);
  view1->setAttributeScalar(size, size_large);

  {
    IndexType idx1 = view1->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view1->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(1, idx2);
    IndexType idx3 = view1->getNextValidAttrValueIndex(idx2);
    EXPECT_EQ(2, idx3);
    IndexType idx4 = view1->getNextValidAttrValueIndex(idx3);
    EXPECT_EQ(InvalidIndex, idx4);
  }

  // set first attribute
  View * view2 = root->createView("view2");
  view2->setAttributeString(color, color_red);

  {
    IndexType idx1 = view2->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view2->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // set last attribute
  View * view3 = root->createView("view3");
  view3->setAttributeScalar(size, size_large);

  {
    IndexType idx1 = view3->getFirstValidAttrValueIndex();
    EXPECT_EQ(2, idx1);
    IndexType idx2 = view3->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // set first and last attributes
  View * view4 = root->createView("view4");
  view4->setAttributeString(color, color_red);
  view4->setAttributeScalar(size, size_large);

  {
    IndexType idx1 = view4->getFirstValidAttrValueIndex();
    EXPECT_EQ(0, idx1);
    IndexType idx2 = view4->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(2, idx2);
    IndexType idx3 = view4->getNextValidAttrValueIndex(idx2);
    EXPECT_EQ(InvalidIndex, idx3);
  }

  // no attributes
  View * view5 = root->createView("view5");

  {
    IndexType idx1 = view5->getFirstValidAttrValueIndex();
    EXPECT_EQ(InvalidIndex, idx1);
    IndexType idx2 = view5->getNextValidAttrValueIndex(idx1);
    EXPECT_EQ(InvalidIndex, idx2);
  }

  // XXX - now delete attribute and check again
}

//------------------------------------------------------------------------------
// save and load attributes from a file

TEST(sidre_attribute,save_attributes)
{
  //  bool ok;
  int idata[5], *bdata;

  const std::string file_path_base("sidre_attribute_datastore_");
  DataStore * ds1 = new DataStore();

  // Create attributes for DataStore
  Attribute * color = ds1->createAttributeString(name_color, color_none);
  EXPECT_TRUE( color != AXOM_NULLPTR );

  Attribute * dump = ds1->createAttributeScalar(name_dump, dump_no);
  EXPECT_TRUE( dump != AXOM_NULLPTR );

  Attribute * size = ds1->createAttributeScalar(name_size, size_small);
  EXPECT_TRUE( size != AXOM_NULLPTR );

  EXPECT_EQ(3, ds1->getNumAttributes());

  Group * root1 = ds1->getRoot();

  // empty
  View * view1a = root1->createView("empty");
  view1a->setAttributeString(color, "color-empty");
  view1a->setAttributeScalar(dump, dump_yes);
  view1a->setAttributeScalar(size, size_small);

  // buffer
  View * view1b = root1->createViewAndAllocate("buffer", INT_ID, 5);
  bdata = view1b->getData();
  view1b->setAttributeString(color, "color-buffer");
  view1b->setAttributeScalar(size, size_medium);

  // external
  View * view1c = root1->createView("external", INT_ID, 5, idata);
  view1c->setAttributeScalar(size, size_large);

  // scalar
  View * view1d = root1->createViewScalar("scalar", 1);
  view1d->setAttributeString(color, "color-scalar");

  // string
  View * view1e = root1->createViewString("string", "value");
  view1e->setAttributeString(color, "color-string");

  // empty without attributes
  root1->createView("empty-no-attributes");

  for (int i=0; i < 5; i++) 
  {
    idata[i] = i;
    bdata[i] = i;
  }

  //----------------------------------------

  for (int i = 0 ; i < nprotocols ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];
    root1->save(file_path, protocols[i]);
  }

  delete ds1;

  //----------------------------------------
  // Only restore conduit_hdf5
  for (int i = 1 ; i < 2 ; ++i)
  {
    const std::string file_path = file_path_base + protocols[i];

    DataStore * ds2 = new DataStore();
    Group * root2 = ds2->getRoot();

    root2->load(file_path, protocols[i]);
    EXPECT_EQ(3, ds2->getNumAttributes());

    // Check available attributes

    Attribute *attr_color = ds2->getAttribute(name_color);
    EXPECT_EQ(color_none, attr_color->getDefaultNodeRef().as_string());

    Attribute *attr_dump = ds2->getAttribute(name_dump);
    EXPECT_EQ(dump_no, attr_dump->getDefaultNodeRef().as_int());

    Attribute *attr_size = ds2->getAttribute(name_size);
    EXPECT_EQ(size_small, attr_size->getDefaultNodeRef().as_double());

    // Check attributes assigned to Views

    View * view2a = root2->getView("empty");
    EXPECT_TRUE(view2a->hasAttributeValue(name_color));
    EXPECT_TRUE(view2a->hasAttributeValue(name_dump));
    EXPECT_TRUE(view2a->hasAttributeValue(name_size));
    EXPECT_TRUE(strcmp("color-empty",
		       view2a->getAttributeString(attr_color)) == 0);
    EXPECT_EQ(dump_yes, view2a->getAttributeScalar<int>(attr_dump));
    EXPECT_EQ(size_small, view2a->getAttributeScalar<double>(attr_size));
    
    View * view2b = root2->getView("buffer");
    EXPECT_TRUE(view2b->hasAttributeValue(name_color));
    EXPECT_FALSE(view2b->hasAttributeValue(name_dump));
    EXPECT_TRUE(view2b->hasAttributeValue(name_size));
    EXPECT_TRUE(strcmp("color-buffer",
		       view2b->getAttributeString(attr_color)) == 0);
    EXPECT_EQ(size_medium, view2b->getAttributeScalar<double>(attr_size));

    View * view2c = root2->getView("external");
    EXPECT_FALSE(view2c->hasAttributeValue(name_color));
    EXPECT_FALSE(view2c->hasAttributeValue(name_dump));
    EXPECT_TRUE(view2c->hasAttributeValue(name_size));
    EXPECT_EQ(size_large, view2c->getAttributeScalar<double>(attr_size));

    View * view2d = root2->getView("scalar");
    EXPECT_TRUE(view2d->hasAttributeValue(name_color));
    EXPECT_FALSE(view2d->hasAttributeValue(name_dump));
    EXPECT_FALSE(view2d->hasAttributeValue(name_size));
    EXPECT_TRUE(strcmp("color-scalar",
		       view2d->getAttributeString(attr_color)) == 0);

    View * view2e = root2->getView("string");
    EXPECT_TRUE(view2e->hasAttributeValue(name_color));
    EXPECT_FALSE(view2e->hasAttributeValue(name_dump));
    EXPECT_FALSE(view2e->hasAttributeValue(name_size));
    EXPECT_TRUE(strcmp("color-string",
		       view2e->getAttributeString(attr_color)) == 0);

    View * view2f = root2->getView("empty-no-attributes");
    EXPECT_FALSE(view2f->hasAttributeValue(name_color));
    EXPECT_FALSE(view2f->hasAttributeValue(name_dump));
    EXPECT_FALSE(view2f->hasAttributeValue(name_size));

    delete ds2;
  }

}

//------------------------------------------------------------------------------
// Internal routine to create a datastore to be used with various iteration schemes
// Error checks here are minimal since it is assumed all of the routines used
// have already been tested.

DataStore *sample_datastore(void)
{
  DataStore * ds = new DataStore();

  // Create all attributes for DataStore
  Attribute * attr_dump = ds->createAttributeScalar(name_dump, dump_no);
  Attribute * attr_color = ds->createAttributeString(name_color, color_none);
  Attribute * attr_animal = ds->createAttributeString(name_animal, animal_none);

  EXPECT_TRUE( attr_dump   != AXOM_NULLPTR );
  EXPECT_TRUE( attr_color  != AXOM_NULLPTR );
  EXPECT_TRUE( attr_animal != AXOM_NULLPTR );

  Group * root = ds->getRoot();

  Group * grpA = root->createGroup("grpA");
  Group * grpB = root->createGroup("grpB");
  Group * grpBB = grpB->createGroup("grpBB");
                 root->createGroup("grpC");  // No Views

  // Tree of empty views
  Group * grpD0 = root->createGroup("grpD0");
  Group * grpD1 = grpD0->createGroup("grpD1");
  Group * grpD2 = grpD1->createGroup("grpD2");
  Group * grpD3 = grpD1->createGroup("grpD3");
  Group * grpD4 = grpD0->createGroup("grpD4");
  Group * grpD5 = grpD4->createGroup("grpD5");

  // checks to silence compiler about unused variables
  EXPECT_TRUE( grpD2 != AXOM_NULLPTR );
  EXPECT_TRUE( grpD3 != AXOM_NULLPTR );
  EXPECT_TRUE( grpD5 != AXOM_NULLPTR );

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

TEST(sidre_attribute,depth_first)
{
  enum nodeclass { GROUP, VIEW };

  struct reference {
    const char *name;
    nodeclass cls;
  };

  DataStore * ds = sample_datastore();

  Group * root = ds->getRoot();

  QueryIterator qitr(root);

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
#if 0
    // check if I have a view and access it :
    qitr.isView();
    View *v = qitr.asView();
 
    // check if I have a group and access it :
    qitr.isGroup();
    Group *g = qitr.asGroup();
 #endif

    // find our current path
    const std::string & name = qitr.getName();
    //std::cout << name << std::endl;

    EXPECT_EQ(order[iorder].name, name);

    if (order[iorder].cls == GROUP)
    {
      EXPECT_TRUE(qitr.isGroup());
      EXPECT_FALSE(qitr.isView());
    }
    else
    {
      EXPECT_FALSE(qitr.isGroup());
      EXPECT_TRUE(qitr.isView());
    }

    qitr.getNext();
    iorder++;
  }

  // Check invalid iterator
  EXPECT_FALSE(qitr.isGroup());
  EXPECT_FALSE(qitr.isView());

  delete ds;
}

//------------------------------------------------------------------------------
