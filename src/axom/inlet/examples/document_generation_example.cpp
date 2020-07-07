// usage : ./document_generation_example_example --enableDocs --deck lua_file --verbose --example [example #]
// example 1: field1 = true; field2 = 5632; NewTable = { str = 'hello'; integer = 32 }
// example 2: foo = false; bar = true; Table1 = { float1 = 3.14; Table1.1 = { Table 1.1.1 = { x = 4 } } }
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
  t = t->addTable("Table1.1", "Table within Table 1");
  t = t->addTable("Table1.1.1", "Table within Table 1.1");
  t->addInt("x", "A variable");
}

void createInletEx3(std::shared_ptr<Inlet> inlet)
{
  auto t = inlet->addTable("Table1", "The first table");
  t->addDouble("float1", " A floating point number in Table 1");
  t = inlet->addTable("Table2", "The second table");
  t->addInt("int1", "An integer in Table 2");
  t = inlet->addTable("Table3", "The third table");
  t->addBool("bool1", "A boolean value in Table 3");
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
  app.add_flag("--verbose",isVerbose, "Enables output of inlet information during run-time");

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  std::shared_ptr<Inlet> inlet = createExampleInlet(&ds, strOption, exampleNum);

  if (docs_enabled) {
    DocWriter doc("example_docs.rst", inlet->sidreGroup(), isVerbose);
  } 

  return 0;
}