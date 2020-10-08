// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>
#include <unordered_map>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

std::shared_ptr<Inlet> createBasicInlet(DataStore* ds,
                                        const std::string& luaString,
                                        bool enableDocs = true)
{
  auto lr = std::make_shared<LuaReader>();
  lr->parseString(luaString);

  return std::make_shared<Inlet>(lr, ds->getRoot(), enableDocs);
}

TEST(inlet_Inlet_basic, getTopLevelBools)
{
  std::string testString = "foo = true; bar = false";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo", "foo's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("bar", "bar's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addBool("non/existant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get_to("bar", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get_to("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getNestedBools)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addBool("foo/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get_to("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get_to("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getDoublyNestedBools)
{
  std::string testString = "foo = { quux = { bar = true; baz = false } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/quux/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/quux/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addBool("foo/quux/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/quux/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get_to("foo/quux/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get_to("foo/quux/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getDeeplyNestedBools)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/quux/corge/quuz/grault/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get_to("foo/quux/corge/quuz/grault/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get_to("foo/quux/corge/quuz/grault/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getNestedBoolsThroughTable)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addBool("foo/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Grab the subtable
  auto table = inlet->getTable("foo");

  // Check for existing fields
  found = table->get_to("bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = table->get_to("baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = table->get_to("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getDeeplyNestedBoolsThroughTable)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  auto table = inlet->getTable("foo/quux/corge");

  // Check for existing fields
  found = table->get_to("quuz/grault/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = table->get_to("quuz/grault/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = table->get_to("quuz/grault/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getDeeplyNestedBoolsThroughField)
{
  std::string testString =
    "foo = { quux = { corge = { quuz = { grault = { bar = true; baz = false } "
    "} } } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  auto table = inlet->getTable("foo/quux/corge");

  // Check for existing fields
  auto bar_field = table->getField("quuz/grault/bar");
  found = bar_field->get_to(value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  auto baz_field = table->getField("quuz/grault/baz");
  found = baz_field->get_to(value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  auto nonexistant_field = table->getField("quuz/grault/nonexistant");
  found = nonexistant_field->get_to(value);
  EXPECT_FALSE(found);
  EXPECT_TRUE(value);
}

TEST(inlet_Inlet_basic, getTopLevelDoubles)
{
  std::string testString = "foo = 5.05; bar = 15.1";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addDouble("foo", "foo's description");
  EXPECT_TRUE(currField);

  currField = inlet->addDouble("bar", "bar's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addDouble("nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  double value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5.05);

  found = inlet->get_to("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15.1);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get_to("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getNestedDoubles)
{
  std::string testString = "foo = { bar = 200.5; baz = 100.987654321 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addDouble("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addDouble("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addDouble("foo/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  double value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200.5);

  found = inlet->get_to("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100.987654321);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get_to("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getTopLevelInts)
{
  std::string testString = "foo = 5; bar = 15";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addInt("foo", "foo's description");
  EXPECT_TRUE(currField);

  currField = inlet->addInt("bar", "bar's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addInt("nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  int value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5);

  found = inlet->get_to("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get_to("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getNestedInts)
{
  std::string testString = "foo = { bar = 200; baz = 100 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addInt("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addInt("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addInt("foo/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  int value = -1;
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200);

  found = inlet->get_to("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get_to("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, 1000);
}

TEST(inlet_Inlet_basic, getTopLevelStrings)
{
  std::string testString = "foo = 'test string'; bar = '15'";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addString("foo", "foo's description");
  EXPECT_TRUE(currField);

  currField = inlet->addString("bar", "bar's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addString("nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  std::string value = "";
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "test string");

  found = inlet->get_to("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "15");

  // Check one that doesn't exist and doesn't have a default value
  value = "don't change";
  found = inlet->get_to("nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, "don't change");
}

TEST(inlet_Inlet_basic, getNestedStrings)
{
  std::string testString = "foo = { bar = 'yet another string'; baz = '' }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addString("foo/bar", "bar's description");
  EXPECT_TRUE(currField);

  currField = inlet->addString("foo/baz", "baz's description");
  EXPECT_TRUE(currField);

  // Check one that doesn't exist and doesn't have a default value
  currField = inlet->addString("foo/nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  std::string value = "";
  bool found = false;

  // Check for existing fields
  found = inlet->get_to("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "yet another string");

  found = inlet->get_to("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "");

  // Check one that doesn't exist and doesn't have a default value
  value = "1000";
  found = inlet->get_to("foo/nonexistant", value);
  EXPECT_FALSE(found);
  EXPECT_EQ(value, "1000");
}

TEST(inlet_Inlet_basic, getNestedValuesAddedUsingTable)
{
  std::string strVal = "";
  int intVal = 0;
  double doubleVal = 0;
  bool boolVal = false;
  bool found = false;

  std::string testString =
    "foo = { bar = 'yet another string'; so = 3.5; re = 9; mi = true }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  auto table = inlet->addTable("foo", "A table called foo");
  table->required(true);

  currField = table->addString("bar", "bar's description");
  EXPECT_TRUE(currField);
  currField->required(true);

  currField = table->addDouble("so", "so's description");
  EXPECT_TRUE(currField);

  currField = table->addInt("re", "re's description");
  EXPECT_TRUE(currField);

  currField = table->addBool("mi", "mi's description");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  found = inlet->get_to("foo/bar", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "yet another string");

  found = inlet->get_to("foo/mi", boolVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(boolVal, true);

  found = inlet->get_to("foo/so", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 3.5);

  found = inlet->get_to("foo/re", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 9);
}

TEST(inlet_Inlet_views, NestedTableViewCheck1)
{
  std::string testString =
    "field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  std::shared_ptr<axom::inlet::Field> currField;
  currField = inlet->addBool("field1", "this is field #1, a boolean value");
  currField->required(true);
  currField = inlet->addInt("field2", "this is field #2, an integer");
  currField->required(false);
  auto t = inlet->addTable("NewTable", "It's blue");
  t->required(false);
  currField = t->addString("str", "str's description");
  currField->required(true);
  currField = t->addInt("integer", "a whole number");
  currField->required(false);

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  bool found;
  bool boolVal;
  int intVal;
  std::string strVal;

  found = inlet->get_to("field1", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("field2", intVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("NewTable/str", strVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("NewTable/integer", intVal);
  EXPECT_TRUE(found);

  EXPECT_TRUE(sidreGroup->hasView("field1/required"));
  EXPECT_TRUE(sidreGroup->hasView("field2/required"));
  EXPECT_TRUE(sidreGroup->hasView("NewTable/str/required"));
  EXPECT_TRUE(sidreGroup->hasView("NewTable/integer/required"));

  EXPECT_TRUE(sidreGroup->hasView("field1/description"));
  EXPECT_TRUE(sidreGroup->hasView("field2/description"));
  EXPECT_TRUE(sidreGroup->hasView("NewTable/str/description"));
  EXPECT_TRUE(sidreGroup->hasView("NewTable/integer/description"));
}

TEST(inlet_Inlet_views, NestedTableViewCheck2)
{
  std::string testString =
    "foo = false; bar = true; Table1 = { float1 = 3.14; Table11 = { Table111 = "
    "{ x = 4 } } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  std::shared_ptr<axom::inlet::Field> currField;
  currField = inlet->addBool("foo", "foo's description");
  currField->required(true);
  currField = inlet->addBool("bar", "bar's description");
  currField->required(false);

  auto t = inlet->addTable("Table1", "The first table");
  t->required(false);
  currField = t->addDouble("float1", "floating point number within table 1");
  currField->required(true);
  t = t->addTable("Table11", "Table within Table 1");
  t = t->addTable("Table111", "Table within Table 11");
  t->addInt("x", "A variable");

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  bool found;
  bool boolVal;
  double doubleVal;

  found = inlet->get_to("foo", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("bar", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("Table1/float1", doubleVal);
  EXPECT_TRUE(found);

  EXPECT_TRUE(sidreGroup->hasView("foo/required"));
  EXPECT_TRUE(sidreGroup->hasView("bar/required"));
  EXPECT_TRUE(sidreGroup->hasView("Table1/float1/required"));

  EXPECT_TRUE(sidreGroup->hasView("foo/description"));
  EXPECT_TRUE(sidreGroup->hasView("bar/description"));
  EXPECT_TRUE(sidreGroup->hasView("Table1/float1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Table1/Table11/Table111/x/description"));
}

TEST(inlet_Inlet_views, NestedTableViewCheck3)
{
  std::string testString =
    "Table1 = { float1 = 5.6 }; Table2 = { int1 = 95 }; Table3 = { bool1 = "
    "true }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto t = inlet->addTable("Table1", "The first table");
  t->addDouble("float1", " A floating point number in Table 1");
  t = inlet->addTable("Table2", "The second table");
  t->addInt("int1", "An integer in Table 2");
  t = inlet->addTable("Table3", "The third table");
  t->addBool("bool1", "A boolean value in Table 3");

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  bool found;
  bool boolVal;
  int intVal;
  double doubleVal;

  found = inlet->get_to("Table1/float1", doubleVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("Table2/int1", intVal);
  EXPECT_TRUE(found);
  found = inlet->get_to("Table3/bool1", boolVal);
  EXPECT_TRUE(found);

  EXPECT_TRUE(sidreGroup->hasView("Table1/float1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Table2/int1/description"));
  EXPECT_TRUE(sidreGroup->hasView("Table3/bool1/description"));
}

TEST(inlet_Inlet, mixLevelTables)
{
  // regression test for the lua stack being at the wrong
  // level with nested tables with non-nested values

  axom::sidre::DataStore dataStore;
  std::string input =
    "thermal_solver={"
    "   u0 = { type = 'function', func = 'BoundaryTemperature'},"
    "   kappa = { type = 'constant', constant = 0.5},"
    "   solver = {"
    "     rel_tol = 1.e-6,"
    "     abs_tol = 1.e-12,"
    "     print_level = 0,"
    "     max_iter = 100,"
    "     dt = 1.0,"
    "     steps = 1 "
    "   }"
    "}";

  auto inlet = createBasicInlet(&dataStore, input);

  //
  // Define input file schema
  //

  // fields that are expected to not be present in above string
  inlet->addString("thermal_solver/mesh/filename", "filename");
  inlet->addInt("thermal_solver/mesh/serial", "serial");
  inlet->addInt("thermal_solver/mesh/parallel", "parallel");
  inlet->addInt("thermal_solver/order", "order");
  inlet->addString("thermal_solver/timestepper", "timestepper");

  // Rest of the fields are present
  inlet->addString("thermal_solver/u0/type", "type");
  inlet->addString("thermal_solver/u0/func", "func");

  inlet->addString("thermal_solver/kappa/type", "type");
  inlet->addDouble("thermal_solver/kappa/constant", "constant");

  inlet->addDouble("thermal_solver/solver/rel_tol", "rel_tol");
  inlet->addDouble("thermal_solver/solver/abs_tol", "abs_tol");
  inlet->addInt("thermal_solver/solver/print_level", "print_level");
  inlet->addInt("thermal_solver/solver/max_iter", "max_iter");
  inlet->addDouble("thermal_solver/solver/dt", "dt");
  inlet->addInt("thermal_solver/solver/steps", "steps");

  //
  //  Verify values found in input file
  //
  bool found = false;
  std::string strVal;
  int intVal;
  double doubleVal;

  // u0
  found = inlet->get_to("thermal_solver/u0/type", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "function");

  found = inlet->get_to("thermal_solver/u0/func", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "BoundaryTemperature");

  // kappa
  found = inlet->get_to("thermal_solver/kappa/type", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "constant");

  found = inlet->get_to("thermal_solver/kappa/constant", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 0.5);

  // solver
  found = inlet->get_to("thermal_solver/solver/rel_tol", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.e-6);

  found = inlet->get_to("thermal_solver/solver/abs_tol", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.e-12);

  found = inlet->get_to("thermal_solver/solver/print_level", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 0);

  found = inlet->get_to("thermal_solver/solver/max_iter", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 100);

  found = inlet->get_to("thermal_solver/solver/dt", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.0);

  found = inlet->get_to("thermal_solver/solver/steps", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 1);
}

TEST(inlet_Field, defaultValuesDocsEnabled)
{
  std::string testString =
    "Table1 = { float1 = 5.6; Table2 = { int1 = 95; Table4 = { str1= 'hi' } } "
    "}; Table3 = { bool1 = true }";
  DataStore ds;
  // Creating basic inlet with documentation enabled (indicated by the last param)
  auto inlet = createBasicInlet(&ds, testString, true);

  // new fields
  inlet->addDouble("field1")->defaultValue(
    2);  // int argument will get casted to double
  inlet->addInt("Table1/Table2/field2")->defaultValue(5);
  inlet->addBool("Table1/field3")->defaultValue(true);
  inlet->addString("Table3/field4")->defaultValue("default for new string");

  // existing fields
  inlet->addInt("Table1/Table2/int1")->defaultValue(100);
  inlet->addBool("Table3/bool1")->defaultValue(false);
  inlet->addDouble("Table1/float1")->defaultValue(3.14);
  inlet->addString("Table1/Table2/Table4/str1")
    ->defaultValue("default for old string");

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  EXPECT_TRUE(sidreGroup->hasView("field1/defaultValue"));
  double doubleVal = sidreGroup->getView("field1/defaultValue")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  EXPECT_TRUE(sidreGroup->hasView("Table1/Table2/field2/defaultValue"));
  int intVal =
    sidreGroup->getView("Table1/Table2/field2/defaultValue")->getScalar();
  EXPECT_EQ(intVal, 5);

  EXPECT_TRUE(sidreGroup->hasView("Table1/field3/defaultValue"));
  int8_t boolVal = sidreGroup->getView("Table1/field3/defaultValue")->getScalar();
  EXPECT_EQ(boolVal, 1);

  EXPECT_TRUE(sidreGroup->hasView("Table3/field4/defaultValue"));
  std::string strVal =
    sidreGroup->getView("Table3/field4/defaultValue")->getString();
  EXPECT_EQ(strVal, "default for new string");

  EXPECT_TRUE(sidreGroup->hasView("Table1/Table2/int1/defaultValue"));
  intVal = sidreGroup->getView("Table1/Table2/int1/defaultValue")->getScalar();
  EXPECT_EQ(intVal, 100);

  EXPECT_TRUE(sidreGroup->hasView("Table3/bool1/defaultValue"));
  boolVal = sidreGroup->getView("Table3/bool1/defaultValue")->getScalar();
  EXPECT_EQ(boolVal, 0);

  EXPECT_TRUE(sidreGroup->hasView("Table1/float1/defaultValue"));
  doubleVal = sidreGroup->getView("Table1/float1/defaultValue")->getScalar();
  EXPECT_EQ(doubleVal, 3.14);

  EXPECT_TRUE(sidreGroup->hasView("Table1/Table2/Table4/str1/defaultValue"));
  strVal =
    sidreGroup->getView("Table1/Table2/Table4/str1/defaultValue")->getString();
  EXPECT_EQ(strVal, "default for old string");

  doubleVal = sidreGroup->getView("field1/value")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  intVal = sidreGroup->getView("Table1/Table2/field2/value")->getScalar();
  EXPECT_EQ(intVal, 5);

  boolVal = sidreGroup->getView("Table1/field3/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  strVal = sidreGroup->getView("Table3/field4/value")->getString();
  EXPECT_EQ(strVal, "default for new string");

  intVal = sidreGroup->getView("Table1/Table2/int1/value")->getScalar();
  EXPECT_EQ(intVal, 95);

  boolVal = sidreGroup->getView("Table3/bool1/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  doubleVal = sidreGroup->getView("Table1/float1/value")->getScalar();
  EXPECT_EQ(doubleVal, 5.6);

  strVal = sidreGroup->getView("Table1/Table2/Table4/str1/value")->getString();
  EXPECT_EQ(strVal, "hi");
}

TEST(inlet_Field, defaultValuesDocsDisabled)
{
  std::string testString =
    "Table1 = { float1 = 5.6; Table2 = { int1 = 95; Table4 = { str1= 'hi' } } "
    "}; Table3 = { bool1 = true }";
  DataStore ds;
  // Creating basic inlet with documentation disabled (indicated by the last param)
  auto inlet = createBasicInlet(&ds, testString, false);

  // new fields
  inlet->addDouble("field1")->defaultValue(2.0);
  inlet->addInt("Table1/Table2/field2")->defaultValue(5);
  inlet->addBool("Table1/field3")->defaultValue(true);
  inlet->addString("Table3/field4")->defaultValue("default for new string");

  // existing fields
  inlet->addInt("Table1/Table2/int1")->defaultValue(100);
  inlet->addBool("Table3/bool1")->defaultValue(false);
  inlet->addDouble("Table1/float1")->defaultValue(3.14);
  inlet->addString("Table1/Table2/Table4/str1")
    ->defaultValue("default for old string");

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  EXPECT_FALSE(sidreGroup->hasView("field1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table1/Table2/field2/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table1/field3/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table3/field4/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table1/Table2/int1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table3/bool1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table1/float1/defaultValue"));
  EXPECT_FALSE(sidreGroup->hasView("Table1/Table2/Table4/str1/defaultValue"));

  double doubleVal = sidreGroup->getView("field1/value")->getScalar();
  EXPECT_EQ(doubleVal, 2.0);

  int intVal = sidreGroup->getView("Table1/Table2/field2/value")->getScalar();
  EXPECT_EQ(intVal, 5);

  int8_t boolVal = sidreGroup->getView("Table1/field3/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  std::string strVal = sidreGroup->getView("Table3/field4/value")->getString();
  EXPECT_EQ(strVal, "default for new string");

  intVal = sidreGroup->getView("Table1/Table2/int1/value")->getScalar();
  EXPECT_EQ(intVal, 95);

  boolVal = sidreGroup->getView("Table3/bool1/value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  doubleVal = sidreGroup->getView("Table1/float1/value")->getScalar();
  EXPECT_EQ(doubleVal, 5.6);

  strVal = sidreGroup->getView("Table1/Table2/Table4/str1/value")->getString();
  EXPECT_EQ(strVal, "hi");
}

TEST(inlet_Field, ranges)
{
  std::string testString =
    "Table1 = { float1 = 5.6; Table2 = { int1 = 95; Table4 = { str1= 'hi' } } "
    "}; Table3 = { bool1 = true }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  inlet->addInt("Table1/set")->validValues({2, 4, 6});
  EXPECT_TRUE(sidreGroup->hasView("Table1/set/validValues"));
  int* bufferArr1 = sidreGroup->getView("Table1/set/validValues")->getArray();
  EXPECT_TRUE(bufferArr1);
  EXPECT_EQ(
    sidreGroup->getView("Table1/set/validValues")->getBuffer()->getNumElements(),
    3);
  EXPECT_EQ(bufferArr1[0], 2);
  EXPECT_EQ(bufferArr1[1], 4);
  EXPECT_EQ(bufferArr1[2], 6);

  inlet->addDouble("Table1/Table2/Table4/double1")->range(2.0, 5.0);
  EXPECT_TRUE(sidreGroup->hasView("Table1/Table2/Table4/double1/range"));
  double* bufferArr2 =
    sidreGroup->getView("Table1/Table2/Table4/double1/range")->getArray();
  EXPECT_TRUE(bufferArr2);
  EXPECT_EQ(sidreGroup->getView("Table1/Table2/Table4/double1/range")
              ->getBuffer()
              ->getNumElements(),
            2);
  EXPECT_EQ(bufferArr2[0], 2.0);
  EXPECT_EQ(bufferArr2[1], 5.0);

  inlet->addInt("Table1/Table2/int1")->range(1, 50);
  EXPECT_TRUE(sidreGroup->hasView("Table1/Table2/int1/range"));
  int* bufferArr3 = sidreGroup->getView("Table1/Table2/int1/range")->getArray();
  EXPECT_TRUE(bufferArr3);
  EXPECT_EQ(
    sidreGroup->getView("Table1/Table2/int1/range")->getBuffer()->getNumElements(),
    2);
  EXPECT_EQ(bufferArr3[0], 1);
  EXPECT_EQ(bufferArr3[1], 50);
}

TEST(inlet_Inlet_verify, verifyRequired)
{
  std::string testString =
    "field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  inlet->addString("NewTable/str")->required(true);
  inlet->addInt("NewTable/int")->required(false);
  inlet->addBool("field1")->required(true);
  inlet->addInt("field2")->required(true);
  EXPECT_TRUE(inlet->verify());

  inlet->addBool("NewTable/field3")->required(true);
  inlet->addTable("NewTable/LastTable")->addDouble("field4")->required(true);
  EXPECT_FALSE(inlet->verify());
}

TEST(inlet_Inlet_verify, verifyDoubleRange)
{
  // For checking values
  std::string testString =
    "field1 = true; field2 = 56.32; NewTable = { str = 'hello'; field4 = 22.19 "
    "}";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field = inlet->addDouble("field2");
  field->range(1.0, 57.2);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addDouble("field3");
  field->range(-10.5, 13.23);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addDouble("NewTable/field4");
  field->range(1.0, 22.0);
  EXPECT_FALSE(inlet->verify());

  // For checking default values
  std::string testString1 =
    "field1 = true; field2 = 56.32; NewTable = { str = 'hello'; field4 = 22.19 "
    "}";
  DataStore ds1;
  auto inlet1 = createBasicInlet(&ds1, testString1);

  auto field1 = inlet1->addDouble("field2");
  field1->defaultValue(5)->range(0, 60);  // ints will get casted to double
  EXPECT_TRUE(inlet1->verify());

  field1 = inlet1->addDouble("field3");
  field1->defaultValue(-12)->range(-10.5, 13.23);
  EXPECT_FALSE(inlet1->verify());
}

TEST(inlet_Inlet_verify, verifyIntRange)
{
  // For checking values
  std::string testString =
    "field1 = true; field2 = 56; NewTable = { field4 = 22; field5 = 48 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field = inlet->addInt("field2");
  field->range(0, 56);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addInt("field3");
  field->range(-12, 13);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addInt("NewTable/field4");
  field->range(1, 23);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addInt("NewTable/field5");
  field->range(1, 7);
  EXPECT_FALSE(inlet->verify());

  // For checking default values
  std::string testString1 =
    "field1 = true; field2 = 56; NewTable = { field4 = 22; field5 = 48 }";
  DataStore ds1;
  auto inlet1 = createBasicInlet(&ds1, testString1);
  field = inlet1->addInt("field2");
  field->range(0, 56)->defaultValue(32);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("field3");
  field->range(-12, 13)->defaultValue(-3);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field4");
  field->range(1, 23)->defaultValue(24);
  EXPECT_FALSE(inlet1->verify());
}

TEST(inlet_Inlet_verify, verifyValidIntValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 56; NewTable = { field4 = 22; field5 = 48 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field = inlet->addInt("field2");
  field->validValues({1, 2, 3, 56, 57, 58});
  EXPECT_TRUE(inlet->verify());

  std::vector<int> nums = {-1, -2, -6, -18, 21};

  field = inlet->addInt("field3");
  field->validValues(nums);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addInt("NewTable/field4");
  field->validValues({44, 40, 48, 22});
  EXPECT_TRUE(inlet->verify());

  field = inlet->addInt("NewTable/field5");
  field->validValues(nums);
  EXPECT_FALSE(inlet->verify());

  // check default values
  std::string testString1 =
    "field1 = true; field2 = 56; NewTable = { field4 = 22; field5 = 48 }";
  DataStore ds1;
  auto inlet1 = createBasicInlet(&ds1, testString1);

  field = inlet1->addInt("field2");
  field->validValues({1, 2, 3, 56, 57, 58})->defaultValue(2);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("field3");
  field->validValues(nums)->defaultValue(21);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field4");
  field->validValues({44, 40, 48, 22})->defaultValue(48);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field5");
  field->validValues(nums)->defaultValue(90);
  EXPECT_FALSE(inlet1->verify());
}

TEST(inlet_Inlet_verify, verifyValidDoubleValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 56.0; NewTable = { field4 = 22.0; field5 = 48.23 "
    "}";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field = inlet->addDouble("field2");
  field->validValues({1, 2, 3, 56, 57, 58});
  EXPECT_TRUE(inlet->verify());

  std::vector<int> nums = {-1, -2, -6, -18, 21};

  field = inlet->addDouble("field3");
  field->validValues(nums);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addDouble("NewTable/field4");
  field->validValues({44.1, 40., 48., 22.});
  EXPECT_TRUE(inlet->verify());

  field = inlet->addDouble("NewTable/field5");
  field->validValues(nums);
  EXPECT_FALSE(inlet->verify());

  // check default values
  std::string testString1 =
    "field1 = true; field2 = 56.0; NewTable = { field4 = 22.0; field5 = 48.23 "
    "}";
  DataStore ds1;
  auto inlet1 = createBasicInlet(&ds1, testString1);

  field = inlet1->addDouble("field2");
  field->validValues({1, 2, 3, 56, 57, 58})->defaultValue(2.);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("field3");
  field->validValues(nums)->defaultValue(21);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("NewTable/field4");
  field->validValues({44.05, 40.13, 48.28, 22.})->defaultValue(48.28);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("NewTable/field5");
  field->validValues(nums)->defaultValue(90.1);
  EXPECT_FALSE(inlet1->verify());
}

TEST(inlet_Inlet_verify, verifyValidStringValues)
{
  // check values
  std::string testString =
    "field1 = true; field2 = 'abc'; NewTable = { field3 = 'xyz'; field4 = "
    "'yes' }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field = inlet->addString("field2");
  field->validValues({"abc", "defg", "hijk", "lm"});
  EXPECT_TRUE(inlet->verify());

  std::vector<std::string> strs = {"nop", "qrstuv", "xyz", "wx"};

  field = inlet->addString("NewTable/field3");
  field->validValues(strs);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addString("Table1/field5");
  field->validValues(strs);
  EXPECT_TRUE(inlet->verify());

  field = inlet->addString("NewTable/field4");
  field->validValues(strs);
  EXPECT_FALSE(inlet->verify());

  // check default values
  std::string testString1 = "field1 = true; NewTable = { field5 = 'nop' }";
  DataStore ds1;
  auto inlet1 = createBasicInlet(&ds1, testString1);

  field = inlet1->addString("field2");
  field->validValues({"abc", "defg", "hijk", "lm"})->defaultValue("defg");
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("field3");
  field->validValues(strs)->defaultValue("wx");
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("NewTable/field4");
  field->validValues(strs);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("NewTable/field5");
  field->validValues(strs)->defaultValue("zyx");
  EXPECT_FALSE(inlet1->verify());
}

TEST(inlet_verify, verifyFieldLambda)
{
  std::string testString =
    "field1 = true; field2 = 'abc'; NewTable = { field3 = 'xyz'; field4 = "
    "'yes' }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field1 = inlet->addBool("field1");
  auto field2 = inlet->addString("field2");
  auto field3 = inlet->addString("NewTable/field3");

  field2->registerVerifier([&]() -> bool {
    std::string str = "";
    inlet->get_to("field2", str);
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet->verify());

  field3->registerVerifier([&]() -> bool {
    std::string str = "";
    inlet->get_to("NewTable/field3", str);
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_FALSE(inlet->verify());
}

TEST(inlet_verify, verifyTableLambda1)
{
  std::string testString =
    "field1 = true; field2 = 'abc'; NewTable = { field3 = 'xyz'; field4 = "
    "'yes' }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  auto field1 = inlet->addBool("field1");
  auto field2 = inlet->addString("field2");
  auto table1 = inlet->addTable("NewTable");
  auto field3 = table1->addString("field3");

  field2->registerVerifier([&]() -> bool {
    std::string str = "";
    inlet->get_to("field2", str);
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet->verify());

  field3->registerVerifier([&]() -> bool {
    std::string str = "";
    inlet->get_to("NewTable/field3", str);
    return (str.size() >= 1 && str[0] == 'x');
  });
  EXPECT_TRUE(inlet->verify());

  EXPECT_TRUE(table1->hasField("field3"));

  table1->registerVerifier([&]() -> bool { return table1->hasField("field3"); });
  EXPECT_TRUE(inlet->verify());

  table1->registerVerifier([&]() -> bool { return table1->hasField("field22"); });
  EXPECT_FALSE(inlet->verify());
}

TEST(inlet_verify, verifyTableLambda2)
{
  std::string testString =
    "thermal_solver={}\n"
    "thermal_solver.order = 2\n"
    "thermal_solver.timestepper = 'quasistatic'\n"
    "solid_solver={}\n"
    "solid_solver.order = 3\n"
    "solid_solver.timestepper = 'BackwardEuler'\n"
    "material={}\n"
    "material.attribute = 1\n"
    "material.thermalview = 'isotropic'\n"
    "material.solidview = 'hyperelastic'";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  auto thermalOrder = inlet->addInt("thermal_solver/order");
  auto thermalTimeStep = inlet->addString("thermal_solver/timestepper");
  auto solidOrder = inlet->addInt("solid_solver/order");
  auto soldTimestep = inlet->addString("solid_solver/timestepper");
  auto attrib = inlet->addInt("material/attribute");
  auto globalTable = inlet->getGlobalTable();
  auto material = globalTable->getTable("material");

  globalTable->registerVerifier([&]() -> bool {
    bool verifySuccess = true;
    if(globalTable->hasTable("thermal_solver") &&
       !material->hasField("thermalview"))
    {
      verifySuccess = false;
    }
    if(globalTable->hasTable("solid_solver") && !material->hasField("solidview"))
    {
      verifySuccess = false;
    }
    return verifySuccess;
  });

  EXPECT_FALSE(inlet->verify());

  auto thermalView = inlet->addString("material/thermalview");
  auto solidView = inlet->addString("material/solidview");

  EXPECT_TRUE(material->hasField("solidview"));
  EXPECT_TRUE(material->hasField("thermalview"));

  EXPECT_TRUE(inlet->verify());

  auto thing = inlet->addTable("test");
  auto thing2 = thing->addTable("test2");
  auto thing5 = thing2->addTable("test3/test4/test5");

  EXPECT_EQ(globalTable->getTable("test")->name(), "test");
  EXPECT_EQ(thing->getTable("test2")->name(), "test/test2");
  EXPECT_EQ(thing2->getTable("test3/test4/test5")->name(),
            "test/test2/test3/test4/test5");
}

TEST(inlet_verify, requiredTable)
{
  std::string testString =
    "thermal_solver={}\n"
    "thermal_solver.order = 2\n"
    "thermal_solver.timestepper = 'quasistatic'\n"
    "solid_solver={}\n"
    "solid_solver.order = 3\n"
    "solid_solver.timestepper = 'BackwardEuler'\n"
    "material={}\n"
    "material.attribute = 1\n"
    "material.thermalview = 'isotropic'\n"
    "material.solidview = 'hyperelastic'";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  auto material = inlet->addTable("material");
  material->required(true);

  EXPECT_FALSE(inlet->verify());

  auto thermalView = inlet->addString("material/thermalview");
  auto solidView = inlet->addString("material/solidview");

  EXPECT_TRUE(inlet->hasField("material/thermalview"));
  EXPECT_TRUE(inlet->hasField("material/solidview"));

  EXPECT_TRUE(inlet->verify());
}

TEST(inlet_verify, verifyTableLambda3)
{
  std::string testString = "dimensions = 2; vector = { x = 1; y = 2; z = 3; }";
  DataStore ds;
  auto myInlet = createBasicInlet(&ds, testString);
  myInlet->addInt("dimensions")->required(true);
  auto v = myInlet->addTable("vector")->required(true);
  v->addInt("x");

  v->registerVerifier([&]() -> bool {
    int dim;
    myInlet->get_to("dimensions", dim);
    int value;  // field value doesnt matter just that it is present in input file
    bool x_present = v->hasField("x") && myInlet->get_to("vector/x", value);
    bool y_present = v->hasField("y") && myInlet->get_to("vector/y", value);
    bool z_present = v->hasField("z") && myInlet->get_to("vector/z", value);
    if(dim == 1 && x_present)
    {
      return true;
    }
    else if(dim == 2 && x_present && y_present)
    {
      return true;
    }
    else if(dim == 3 && x_present && y_present && z_present)
    {
      return true;
    }
    return false;
  });

  EXPECT_FALSE(myInlet->verify());

  v->addInt("y");
  v->addInt("z");

  EXPECT_TRUE(myInlet->verify());
}

// Checks all of the Table::getArray functions
TEST(inletArrays, getArray)
{
  DataStore ds;
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4}, "
    "              arr2 = {[4] = true, [8] = false}, "
    "              arr3 = {[33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  auto inlet = createBasicInlet(&ds, testString);

  std::unordered_map<int, bool> boolMap;
  std::unordered_map<int, int> intMap;
  std::unordered_map<int, double> doubleMap;
  std::unordered_map<int, std::string> strMap;
  auto arr1 = inlet->getGlobalTable()->addIntArray("luaArrays/arr1");
  auto arr2 = inlet->getGlobalTable()->addBoolArray("luaArrays/arr2");
  auto arr3 = inlet->getGlobalTable()->addStringArray("luaArrays/arr3");
  auto arr4 = inlet->getGlobalTable()->addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};

  EXPECT_TRUE(arr1->getArray(intMap));
  EXPECT_EQ(intMap, expectedInts);

  EXPECT_TRUE(arr2->getArray(boolMap));
  EXPECT_EQ(boolMap, expectedBools);

  EXPECT_TRUE(arr3->getArray(strMap));
  EXPECT_EQ(strMap, expectedStrs);

  EXPECT_TRUE(arr4->getArray(doubleMap));
  EXPECT_EQ(doubleMap, expectedDoubles);
}

// Checks the underlying Sidre representation of the arrays added from Lua
TEST(inletArrays, inletArraysInSidre)
{
  DataStore ds;
  std::string testString =
    "luaArrays = { arr1 = { [1] = 4, [2] = 5, [3] = 6 , [12] = 2.4}, "
    "              arr2 = { [4] = true, [8] = false}, "
    "              arr3 = { [33] = 'hello', [2] = 'bye'}, "
    "              arr4 = { [12] = 2.4 } }";
  auto inlet = createBasicInlet(&ds, testString);

  inlet->addIntArray("luaArrays/arr1");

  auto group = inlet->getGlobalTable()->sidreGroup()->getGroup(
    "luaArrays/arr1/_inlet_array");
  auto idx = group->getGroup("1");
  EXPECT_TRUE(idx);
  int val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 4);

  idx = group->getGroup("12");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 2);

  idx = group->getGroup("2");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 5);

  idx = group->getGroup("3");
  EXPECT_TRUE(idx);
  val = idx->getView("value")->getScalar();
  EXPECT_EQ(val, 6);

  inlet->getGlobalTable()->addBoolArray("luaArrays/arr2");
  group = inlet->getTable("luaArrays/arr2/_inlet_array")->sidreGroup();

  idx = group->getGroup("4");
  EXPECT_TRUE(idx);
  int8_t boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 1);

  idx = group->getGroup("8");
  EXPECT_TRUE(idx);
  boolVal = idx->getView("value")->getScalar();
  EXPECT_EQ(boolVal, 0);

  inlet->getGlobalTable()->addStringArray("luaArrays/arr3");
  group = inlet->getTable("luaArrays/arr3/_inlet_array")->sidreGroup();

  idx = group->getGroup("33");
  EXPECT_TRUE(idx);
  std::string str = idx->getView("value")->getString();
  EXPECT_EQ(str, "hello");

  idx = group->getGroup("2");
  EXPECT_TRUE(idx);
  str = idx->getView("value")->getString();
  EXPECT_EQ(str, "bye");

  inlet->getGlobalTable()->addDoubleArray("luaArrays/arr4");
  group = inlet->getTable("luaArrays/arr4/_inlet_array")->sidreGroup();

  idx = group->getGroup("12");
  EXPECT_TRUE(idx);
  double doubleVal = idx->getView("value")->getScalar();
  EXPECT_EQ(doubleVal, 2.4);
}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
