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
  bool has_name = ds->hasAttribute("units");
  EXPECT_FALSE( has_name );

  Attribute * units = ds->createAttribute("units"); // , STRING, " ");
  EXPECT_TRUE( units != AXOM_NULLPTR );

  IndexType attr_index = units->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute("units");
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

  Attribute * attr = ds->getAttribute("units");
  EXPECT_EQ(attr, units);
  const Attribute * attrc = ds->getAttribute("units");
  EXPECT_EQ(attrc, units);

  attr = ds->getAttribute(0);
  EXPECT_EQ(attr, units);
  attrc = ds->getAttribute(0);
  EXPECT_EQ(attrc, units);

  ds->destroyAttribute(units);

}

//------------------------------------------------------------------------------

#if 0
TEST(sidre_attribute,view_attr)
{
  DataStore * ds = new DataStore();

  int nattrs = ds->getNumAttributes();
  EXPECT_EQ(0, nattrs);

  bool has_index = ds->hasAttribute(0);
  EXPECT_FALSE( has_index );
  bool has_name = ds->hasAttribute("units");
  EXPECT_FALSE( has_name );

  Attribute * units = ds->createAttribute("units"); // , STRING, " ");
  EXPECT_TRUE( units != AXOM_NULLPTR );

  IndexType attr_index = units->getIndex();
  EXPECT_EQ(0, attr_index);

  nattrs = ds->getNumAttributes();
  EXPECT_EQ(1, nattrs);

  has_name = ds->hasAttribute("units");
  EXPECT_TRUE( has_name );
  has_index = ds->hasAttribute(0);
  EXPECT_TRUE( has_index );

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
#endif
//------------------------------------------------------------------------------
