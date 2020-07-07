// usage : ./document_generation_example_example --enableDocs --deck lua_file --verbose --example [example #]
// example 1: field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }
// example 2: foo = false; bar = true; Table1 = { float1 = 3.14; Table11 = { Table111 = { x = 4 } } }
// example 3: Table1 = { float1 = 5.6 }; Table2 = { int1 = 95 }; Table3 = { bool1 = true }

#include "axom/inlet/DocWriter.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/LuaReader.hpp"

#include "CLI11/CLI11.hpp"
#include <iostream>

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;
using axom::inlet::DocWriter;

void createInletEx1(std::shared_ptr<Inlet> inlet)
{
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
  SLIC_ASSERT_MSG(found == true, "field1 not found");
  found = inlet->get("field2", intVal);
  SLIC_ASSERT_MSG(found == true, "field2 not found");
  found = inlet->get("NewTable/str", strVal);
  SLIC_ASSERT_MSG(found == true, "NewTable/str not found");
  found = inlet->get("NewTable/integer", intVal);
  SLIC_ASSERT_MSG(found == true, "NewTable/integer not found");

  SLIC_ASSERT_MSG(sidreGroup->hasView("field1/required"),"field1/required not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("field2/required"),"field2/required not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("NewTable/str/required"),"NewTable/str/required not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("NewTable/integer/required"),"NewTable/integer/required not found");

  SLIC_ASSERT_MSG(sidreGroup->hasView("field1/description"), "field1/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("field2/description"), "field2/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("NewTable/str/description"), "NewTable/str/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("NewTable/integer/description"), "NewTable/integer/description not found");
}

void createInletEx2(std::shared_ptr<Inlet> inlet)
{
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
  int intVal;
  double doubleVal;  

  found = inlet->get("foo", boolVal);
  SLIC_ASSERT_MSG(found == true, "foo not found");
  found = inlet->get("bar", boolVal);
  SLIC_ASSERT_MSG(found == true, "bar not found");
  found = inlet->get("Table1/float1", doubleVal);
  SLIC_ASSERT_MSG(found == true, "Table1/float1 not found");

  SLIC_ASSERT_MSG(sidreGroup->hasView("foo/required"),"foo/required not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("bar/required"),"bar/required not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("Table1/float1/required"),"Table1float1/required not found");

  SLIC_ASSERT_MSG(sidreGroup->hasView("foo/description"), "foo/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("bar/description"), "bar/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("Table1/float1/description"), "Table1/float1/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("Table1/Table11/Table111/x/description"), "Table1/Table11/Table111/x/description not found");
}

void createInletEx3(std::shared_ptr<Inlet> inlet)
{
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
  SLIC_ASSERT_MSG(found == true, "Table1/float1 not found");
  found = inlet->get("Table2/int1", intVal);
  SLIC_ASSERT_MSG(found == true, "Table2/int1 not found");
  found = inlet->get("Table3/bool1", boolVal);
  SLIC_ASSERT_MSG(found == true, "Table3/bool1 not found");

  SLIC_ASSERT_MSG(sidreGroup->hasView("Table1/float1/description"),"Table1/float1/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("Table2/int1/description"), "Table2/int1/description not found");
  SLIC_ASSERT_MSG(sidreGroup->hasView("Table3/bool1/description"), "Table3/bool1/description not found");
}

std::shared_ptr<Inlet> createExampleInlet(DataStore* ds,
                                        const std::string& luaFile, int exNum)
{
  auto lr = std::make_shared<LuaReader>();
  // lr->parseString(luaString);
  lr->parseFile(luaFile);
  auto inlet = std::make_shared<Inlet>(lr, ds->getRoot());
  if (exNum == 1) {
    createInletEx1(inlet);
  } else if (exNum == 2) {
    createInletEx2(inlet);
  } else {
    createInletEx3(inlet);
  }
  return inlet;
}

int main(int argc, char** argv) {
  CLI::App app {"Description here"};
  bool docs_enabled{false};
  app.add_flag("--enableDocs", docs_enabled, "Enables documentation generation");

  std::string strOption;
  auto opt = app.add_option("--deck", strOption, "Path to input deck file");
  opt->check(CLI::ExistingFile);

  int exampleNum{1};
  app.add_option("--example", exampleNum, "Example number to be run");

  bool isVerbose{false};
  app.add_flag("--verbose", isVerbose, "Enables output of inlet information during run-time");

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  std::shared_ptr<Inlet> inlet = createExampleInlet(&ds, strOption, exampleNum);

  std::cout << "____________________________________________" << std::endl;
  if (docs_enabled) {
    DocWriter doc("example_docs.rst", inlet->sidreGroup(), isVerbose);
    std::cout << "Documentation was written to example_docs.rst" << std::endl;
  } else {
    std::cout << "Documentation generation was not enabled; "
              << "if you would like documentation generated, please rerun using the --enableDocs flag.\n";
  }

  return 0;
}