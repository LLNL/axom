// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include <iostream>

#include "axom/sidre.hpp"

#include "axom/inlet/LuaReader.hpp"
#include "axom/inlet/Inlet.hpp"


using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;


std::shared_ptr<Inlet> createBasicInlet(DataStore* ds,
                                        const std::string& luaString)
{
  auto lr = std::make_shared<LuaReader>();
  lr->parseString(luaString);

  return std::make_shared<Inlet>(lr, ds->getRoot());
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
  currField = inlet->addBool("nonexistant", "nothing");
  EXPECT_TRUE(currField);

  //
  // Check stored values from get
  //

  bool value = false;
  bool found = false;

  // Check for existing fields
  found = inlet->get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get("nonexistant", value);
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
  found = inlet->get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_TRUE(value);

  found = inlet->get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_FALSE(value);

  // Check one that doesn't exist and doesn't have a default value
  value = true;
  found = inlet->get("foo/nonexistant", value);
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
  found = inlet->get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5.05);

  found = inlet->get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15.1);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get("nonexistant", value);
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
  found = inlet->get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200.5);

  found = inlet->get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100.987654321);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get("foo/nonexistant", value);
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
  found = inlet->get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 5);

  found = inlet->get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 15);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get("nonexistant", value);
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
  found = inlet->get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 200);

  found = inlet->get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, 100);

  // Check one that doesn't exist and doesn't have a default value
  value = 1000;
  found = inlet->get("foo/nonexistant", value);
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
  found = inlet->get("foo", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "test string");

  found = inlet->get("bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "15");

  // Check one that doesn't exist and doesn't have a default value
  value = "don't change";
  found = inlet->get("nonexistant", value);
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
  found = inlet->get("foo/bar", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "yet another string");

  found = inlet->get("foo/baz", value);
  EXPECT_TRUE(found);
  EXPECT_EQ(value, "");

  // Check one that doesn't exist and doesn't have a default value
  value = "1000";
  found = inlet->get("foo/nonexistant", value);
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

  std::string testString = "foo = { bar = 'yet another string'; so = 3.5; re = 9; mi = true }";
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

  found = inlet->get("foo/bar", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "yet another string");

  found = inlet->get("foo/mi", boolVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(boolVal, true);

  found = inlet->get("foo/so", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 3.5);

  found = inlet->get("foo/re", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 9);

}

TEST(inlet_Inlet_views, NestedTableViewCheck1)
{
  std::string testString = "field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }";
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

  found = inlet->get("field1", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get("field2", intVal);
  EXPECT_TRUE(found);
  found = inlet->get("NewTable/str", strVal);
  EXPECT_TRUE(found);  
  found = inlet->get("NewTable/integer", intVal);
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
  std::string testString = "foo = false; bar = true; Table1 = { float1 = 3.14; Table11 = { Table111 = { x = 4 } } }";
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

  found = inlet->get("foo", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get("bar", boolVal);
  EXPECT_TRUE(found);
  found = inlet->get("Table1/float1", doubleVal);
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
  std::string testString = "Table1 = { float1 = 5.6 }; Table2 = { int1 = 95 }; Table3 = { bool1 = true }";
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

  found = inlet->get("Table1/float1", doubleVal);
  ASSERT_TRUE(found);
  found = inlet->get("Table2/int1", intVal);
  ASSERT_TRUE(found);  
  found = inlet->get("Table3/bool1", boolVal);
  ASSERT_TRUE(found);

  ASSERT_TRUE(sidreGroup->hasView("Table1/float1/description"));
  ASSERT_TRUE(sidreGroup->hasView("Table2/int1/description"));
  ASSERT_TRUE(sidreGroup->hasView("Table3/bool1/description"));
}

TEST(inlet_Inlet, mixLevelTables)
{
  // regression test for the lua stack being at the wrong
  // level with nested tables with non-nested values

  axom::sidre::DataStore dataStore;
  auto luaReader = std::make_shared<axom::inlet::LuaReader>();
  luaReader->parseString(
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
    "}");

  auto inlet = std::make_shared<axom::inlet::Inlet>(luaReader, dataStore.getRoot());

  //
  // Define input deck schema
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
  //  Verify values found in input deck
  //
  bool found = false;
  std::string strVal;
  int intVal;
  double doubleVal;

  // u0
  found = inlet->get("thermal_solver/u0/type", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "function");

  found = inlet->get("thermal_solver/u0/func", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "BoundaryTemperature");

  // kappa
  found = inlet->get("thermal_solver/kappa/type", strVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(strVal, "constant");

  found = inlet->get("thermal_solver/kappa/constant", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 0.5);

  // solver
  found = inlet->get("thermal_solver/solver/rel_tol", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.e-6);

  found = inlet->get("thermal_solver/solver/abs_tol", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.e-12);

  found = inlet->get("thermal_solver/solver/print_level", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 0);

  found = inlet->get("thermal_solver/solver/max_iter", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 100);

  found = inlet->get("thermal_solver/solver/dt", doubleVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(doubleVal, 1.0);

  found = inlet->get("thermal_solver/solver/steps", intVal);
  EXPECT_TRUE(found);
  EXPECT_EQ(intVal, 1);
}

TEST(inlet_Inlet_verify, checkVerify) {
  std::string testString = "field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  std::shared_ptr<axom::inlet::Field> currField;
  currField = inlet->addBool("NewTable/field1");
  ds.getRoot()->getGroup("NewTable")->createViewString("required", "string");
  EXPECT_FALSE(inlet->verify());
  
  ds.getRoot()->getView("NewTable/required")->rename("renamedView");
  currField = inlet->addString("NewTable/str");
  currField->required(true);
  currField = inlet->addString("NewTable/str");
  EXPECT_TRUE(inlet->verify());
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
