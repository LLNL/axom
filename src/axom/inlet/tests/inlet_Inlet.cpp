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

using axom::inlet::Field;
using axom::inlet::Inlet;
using axom::inlet::InletType;
using axom::inlet::LuaReader;
using axom::inlet::Proxy;
using axom::inlet::Table;
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

  // Check for existing fields
  auto currVerifiable = inlet->addBool("foo", "foo's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addBool("bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addBool("non/existant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = (*inlet)["foo"];
  EXPECT_TRUE(value);

  value = inlet->get<bool>("foo");
  EXPECT_TRUE(value);

  value = (*inlet)["bar"];
  EXPECT_FALSE(value);

  value = inlet->get<bool>("bar");
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto table = (*inlet)["non/existant"];
  EXPECT_EQ(table.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedBools)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addBool("foo/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = (*inlet)["foo/bar"];
  EXPECT_TRUE(value);

  value = (*inlet)["foo/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto table = (*inlet)["foo/nonexistant"];
  EXPECT_EQ(table.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getDoublyNestedBools)
{
  std::string testString = "foo = { quux = { bar = true; baz = false } }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addBool("foo/quux/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addBool("foo/quux/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addBool("foo/quux/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = (*inlet)["foo/quux/bar"];
  EXPECT_TRUE(value);

  value = (*inlet)["foo/quux/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto table = (*inlet)["foo/quux/nonexistant"];
  EXPECT_EQ(table.type(), InletType::Nothing);
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

  // Check for existing fields
  auto currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  // Check for existing fields
  value = (*inlet)["foo/quux/corge/quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = (*inlet)["foo/quux/corge/quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto table = (*inlet)["foo/quux/corge/quuz/grault/nonexistant"];
  EXPECT_EQ(table.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedBoolsThroughTable)
{
  std::string testString = "foo = { bar = true; baz = false }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addBool("foo/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addBool("foo/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addBool("foo/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  // Grab the subtable
  auto table = inlet->getTable("foo");

  // Check for existing fields
  value = (*table)["bar"];
  EXPECT_TRUE(value);

  value = (*table)["baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*table)["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
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

  // Check for existing fields
  auto currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  auto table = inlet->getTable("foo/quux/corge");

  // Check for existing fields
  value = (*table)["quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = (*table)["quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*table)["quuz/grault/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
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

  // Check for existing fields
  auto currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable =
    inlet->addBool("foo/quux/corge/quuz/grault/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  bool value = false;

  auto table = inlet->getTable("foo/quux/corge");

  // Check for existing fields
  value = (*table)["quuz/grault/bar"];
  EXPECT_TRUE(value);

  value = (*table)["quuz/grault/baz"];
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  auto nonexistant_field = (*table)["quuz/grault/nonexistant"];
  EXPECT_EQ(nonexistant_field.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getTopLevelDoubles)
{
  std::string testString = "foo = 5.05; bar = 15.1";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addDouble("foo", "foo's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addDouble("bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addDouble("nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  double value = -1;

  // Check for existing fields
  value = (*inlet)["foo"];
  EXPECT_EQ(value, 5.05);

  value = (*inlet)["bar"];
  EXPECT_EQ(value, 15.1);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedDoubles)
{
  std::string testString = "foo = { bar = 200.5; baz = 100.987654321 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addDouble("foo/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addDouble("foo/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addDouble("foo/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  double value = -1;

  // Check for existing fields
  value = (*inlet)["foo/bar"];
  EXPECT_EQ(value, 200.5);

  value = (*inlet)["foo/baz"];
  EXPECT_EQ(value, 100.987654321);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getTopLevelInts)
{
  std::string testString = "foo = 5; bar = 15";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addInt("foo", "foo's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addInt("bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addInt("nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  value = (*inlet)["foo"];
  EXPECT_EQ(value, 5);

  value = (*inlet)["bar"];
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedInts)
{
  std::string testString = "foo = { bar = 200; baz = 100 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addInt("foo/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addInt("foo/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addInt("foo/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  int value = -1;

  // Check for existing fields
  value = (*inlet)["foo/bar"];
  EXPECT_EQ(value, 200);

  value = (*inlet)["foo/baz"];
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getTopLevelStrings)
{
  std::string testString = "foo = 'test string'; bar = '15'";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addString("foo", "foo's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addString("bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addString("nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = (*inlet)["foo"].get<std::string>();
  EXPECT_EQ(value, "test string");

  value = (*inlet)["bar"].get<std::string>();
  EXPECT_EQ(value, "15");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedStrings)
{
  std::string testString = "foo = { bar = 'yet another string'; baz = '' }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  // Check for existing fields
  auto currVerifiable = inlet->addString("foo/bar", "bar's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = inlet->addString("foo/baz", "baz's description");
  EXPECT_TRUE(currVerifiable);

  // Check one that doesn't exist and doesn't have a default value
  currVerifiable = inlet->addString("foo/nonexistant", "nothing");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  std::string value = "";

  // Check for existing fields
  value = (*inlet)["foo/bar"].get<std::string>();
  EXPECT_EQ(value, "yet another string");

  value = (*inlet)["foo/baz"].get<std::string>();
  EXPECT_EQ(value, "");

  // Check one that doesn't exist and doesn't have a default value
  auto proxy = (*inlet)["foo/nonexistant"];
  EXPECT_EQ(proxy.type(), InletType::Nothing);
}

TEST(inlet_Inlet_basic, getNestedValuesAddedUsingTable)
{
  std::string strVal = "";
  int intVal = 0;
  double doubleVal = 0;
  bool boolVal = false;

  std::string testString =
    "foo = { bar = 'yet another string'; so = 3.5; re = 9; mi = true }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  // Check for existing fields
  auto table = inlet->addTable("foo", "A table called foo");
  table->required(true);

  auto currVerifiable = table->addString("bar", "bar's description");
  EXPECT_TRUE(currVerifiable);
  currVerifiable->required(true);

  currVerifiable = table->addDouble("so", "so's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = table->addInt("re", "re's description");
  EXPECT_TRUE(currVerifiable);

  currVerifiable = table->addBool("mi", "mi's description");
  EXPECT_TRUE(currVerifiable);

  //
  // Check stored values from get
  //

  strVal = (*inlet)["foo/bar"].get<std::string>();
  EXPECT_EQ(strVal, "yet another string");

  boolVal = (*inlet)["foo/mi"];
  EXPECT_EQ(boolVal, true);

  doubleVal = (*inlet)["foo/so"];
  EXPECT_EQ(doubleVal, 3.5);

  intVal = (*inlet)["foo/re"];
  EXPECT_EQ(intVal, 9);
}

TEST(inlet_Inlet_views, NestedTableViewCheck1)
{
  std::string testString =
    "field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  auto currVerifiable =
    inlet->addBool("field1", "this is field #1, a boolean value");
  currVerifiable->required(true);
  currVerifiable = inlet->addInt("field2", "this is field #2, an integer");
  currVerifiable->required(false);
  auto t = inlet->addTable("NewTable", "It's blue");
  t->required(false);
  currVerifiable = t->addString("str", "str's description");
  currVerifiable->required(true);
  currVerifiable = t->addInt("integer", "a whole number");
  currVerifiable->required(false);

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  EXPECT_EQ((*inlet)["field1"].type(), InletType::Bool);
  EXPECT_EQ((*inlet)["field2"].type(), InletType::Integer);
  EXPECT_EQ((*inlet)["NewTable/str"].type(), InletType::String);
  EXPECT_EQ((*inlet)["NewTable/integer"].type(), InletType::Integer);

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
  auto currVerifiable = inlet->addBool("foo", "foo's description");
  currVerifiable->required(true);
  currVerifiable = inlet->addBool("bar", "bar's description");
  currVerifiable->required(false);

  auto t = inlet->addTable("Table1", "The first table");
  t->required(false);
  currVerifiable =
    t->addDouble("float1", "floating point number within table 1");
  currVerifiable->required(true);
  t = t->addTable("Table11", "Table within Table 1");
  t = t->addTable("Table111", "Table within Table 11");
  t->addInt("x", "A variable");

  axom::sidre::Group* sidreGroup = inlet->sidreGroup();

  EXPECT_EQ((*inlet)["foo"].type(), InletType::Bool);
  EXPECT_EQ((*inlet)["bar"].type(), InletType::Bool);
  EXPECT_EQ((*inlet)["Table1/float1"].type(), InletType::Double);

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

  EXPECT_EQ((*inlet)["Table1/float1"].type(), InletType::Double);
  EXPECT_EQ((*inlet)["Table2/int1"].type(), InletType::Integer);
  EXPECT_EQ((*inlet)["Table3/bool1"].type(), InletType::Bool);

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
  std::string strVal;
  int intVal;
  double doubleVal;

  // u0
  strVal = (*inlet)["thermal_solver/u0/type"].get<std::string>();
  EXPECT_EQ(strVal, "function");

  strVal = (*inlet)["thermal_solver/u0/func"].get<std::string>();
  EXPECT_EQ(strVal, "BoundaryTemperature");

  // kappa
  strVal = (*inlet)["thermal_solver/kappa/type"].get<std::string>();
  EXPECT_EQ(strVal, "constant");

  doubleVal = (*inlet)["thermal_solver/kappa/constant"];
  EXPECT_EQ(doubleVal, 0.5);

  // solver
  doubleVal = (*inlet)["thermal_solver/solver/rel_tol"];
  EXPECT_EQ(doubleVal, 1.e-6);

  doubleVal = (*inlet)["thermal_solver/solver/abs_tol"];
  EXPECT_EQ(doubleVal, 1.e-12);

  intVal = (*inlet)["thermal_solver/solver/print_level"];
  EXPECT_EQ(intVal, 0);

  intVal = (*inlet)["thermal_solver/solver/max_iter"];
  EXPECT_EQ(intVal, 100);

  doubleVal = (*inlet)["thermal_solver/solver/dt"];
  EXPECT_EQ(doubleVal, 1.0);

  intVal = (*inlet)["thermal_solver/solver/steps"];
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
  field1->defaultValue(5).range(0, 60);  // ints will get casted to double
  EXPECT_TRUE(inlet1->verify());

  field1 = inlet1->addDouble("field3");
  field1->defaultValue(-12).range(-10.5, 13.23);
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
  field->range(0, 56).defaultValue(32);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("field3");
  field->range(-12, 13).defaultValue(-3);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field4");
  field->range(1, 23).defaultValue(24);
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
  field->validValues({1, 2, 3, 56, 57, 58}).defaultValue(2);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("field3");
  field->validValues(nums).defaultValue(21);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field4");
  field->validValues({44, 40, 48, 22}).defaultValue(48);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addInt("NewTable/field5");
  field->validValues(nums).defaultValue(90);
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
  field->validValues({1, 2, 3, 56, 57, 58}).defaultValue(2.);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("field3");
  field->validValues(nums).defaultValue(21);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("NewTable/field4");
  field->validValues({44.05, 40.13, 48.28, 22.}).defaultValue(48.28);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addDouble("NewTable/field5");
  field->validValues(nums).defaultValue(90.1);
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
  field->validValues({"abc", "defg", "hijk", "lm"}).defaultValue("defg");
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("field3");
  field->validValues(strs).defaultValue("wx");
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("NewTable/field4");
  field->validValues(strs);
  EXPECT_TRUE(inlet1->verify());

  field = inlet1->addString("NewTable/field5");
  field->validValues(strs).defaultValue("zyx");
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

  field2->registerVerifier([](const Field& field) -> bool {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet->verify());

  field3->registerVerifier([](const Field& field) -> bool {
    auto str = field.get<std::string>();
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

  field2->registerVerifier([](const Field& field) -> bool {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'a');
  });
  EXPECT_TRUE(inlet->verify());

  field3->registerVerifier([](const Field& field) -> bool {
    auto str = field.get<std::string>();
    return (str.size() >= 1 && str[0] == 'x');
  });
  EXPECT_TRUE(inlet->verify());

  EXPECT_TRUE(table1->hasField("field3"));

  table1->registerVerifier(
    [](const Table& table) -> bool { return table.contains("field3"); });
  EXPECT_TRUE(inlet->verify());

  table1->registerVerifier(
    [](const Table& table) -> bool { return table.contains("field22"); });
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

  globalTable->registerVerifier([](const Table& table) -> bool {
    bool verifySuccess = true;
    if(table.contains("thermal_solver") &&
       !table["material"].contains("thermalview"))
    {
      verifySuccess = false;
    }
    if(table.contains("solid_solver") && !table["material"].contains("solidview"))
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
  auto v = myInlet->addTable("vector");
  v->required(true);
  v->addInt("x");

  v->registerVerifier([&myInlet](const Table& table) -> bool {
    int dim = (*myInlet)["dimensions"];
    bool x_present =
      table.contains("x") && (table["x"].type() == InletType::Integer);
    bool y_present =
      table.contains("y") && (table["y"].type() == InletType::Integer);
    bool z_present =
      table.contains("z") && (table["z"].type() == InletType::Integer);
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

  auto arr1 = inlet->getGlobalTable()->addIntArray("luaArrays/arr1");
  auto arr2 = inlet->getGlobalTable()->addBoolArray("luaArrays/arr2");
  auto arr3 = inlet->getGlobalTable()->addStringArray("luaArrays/arr3");
  auto arr4 = inlet->getGlobalTable()->addDoubleArray("luaArrays/arr4");

  std::unordered_map<int, int> expectedInts {{1, 4}};
  std::unordered_map<int, bool> expectedBools {{4, true}, {8, false}};
  std::unordered_map<int, double> expectedDoubles {{12, 2.4}};
  std::unordered_map<int, std::string> expectedStrs {{33, "hello"}, {2, "bye"}};

  std::unordered_map<int, int> intMap = (*inlet)["luaArrays/arr1"];
  std::unordered_map<int, bool> boolMap = (*inlet)["luaArrays/arr2"];
  std::unordered_map<int, std::string> strMap = (*inlet)["luaArrays/arr3"];
  std::unordered_map<int, double> doubleMap = (*inlet)["luaArrays/arr4"];

  EXPECT_EQ(intMap, expectedInts);

  EXPECT_EQ(boolMap, expectedBools);

  EXPECT_EQ(strMap, expectedStrs);

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
