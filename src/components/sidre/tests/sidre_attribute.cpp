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
  DataStore * ds = new DataStore();

  int nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  bool has_index = ds->hasAttribute(0);
  EXPECT_FALSE( has_index );
  std::string nameattr("units");
  bool has_name = ds->hasAttribute(nameattr);
  EXPECT_FALSE( has_name );

  Attribute * units = ds->createAttribute(nameattr); // , STRING, " ");
  EXPECT_TRUE( units != AXOM_NULLPTR );

  IndexType attr_index = units->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute(nameattr);
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

  Attribute * attr = ds->getAttribute(nameattr);
  EXPECT_EQ(attr, units);
  const Attribute * attrc = ds->getAttribute(nameattr);
  EXPECT_EQ(attrc, units);

  attr = ds->getAttribute(0);
  EXPECT_EQ(attr, units);
  attrc = ds->getAttribute(0);
  EXPECT_EQ(attrc, units);

  ds->destroyAttribute(units);
  nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);
  has_name = ds->hasAttribute(nameattr);
  EXPECT_FALSE( has_name );
  // At this point units points to deallocated memory

  // Create additional attributes
  std::string namedump1("dump1");
  Attribute * dump1 = ds->createAttribute(namedump1); // , STRING, " ");
  EXPECT_TRUE( dump1 != AXOM_NULLPTR );

  attr_index = dump1->getIndex();
  EXPECT_EQ(0, attr_index);

  std::string namedump2("dump2");
  Attribute * dump2 = ds->createAttribute(namedump2); // , STRING, " ");
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
  DataStore * ds = new DataStore();

  std::string nameattr("units");
  Attribute * units = ds->createAttribute(nameattr); // , STRING, " ");
  EXPECT_TRUE( units != AXOM_NULLPTR );

  Group * root = ds->getRoot();
  View  * view1 = root->createView("v1");
  EXPECT_TRUE( view1 != AXOM_NULLPTR );

  bool ok = view1->setAttribute(units, "miles");
  EXPECT_TRUE( ok );

  view1->getAttribute(units);
  //std::string & sattr = view1->getAttribute(units);

#if 0
  IndexType iunits = units->getIndex();

  Group * root = ds.getRoot();
  View * view = root->createView("var1");

  view->setAttribute(units, "miles");
  view->setAttribute(iunits, "miles");
  view->setAttribute("units", "miles");

  const char * attr1 = view->getAttribute(units);
  const char * attr2 = view->getAttribute(iunits);
  const char * attr3 = view->getAttribute("units");


  view = root->createView("var2");
  // Get attributes without setting returns default value
  const char * attr1 = view->getAttribute(units);
#endif
}

//------------------------------------------------------------------------------
