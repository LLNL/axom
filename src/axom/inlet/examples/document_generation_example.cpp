#include "axom/inlet/DocWriter.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/LuaReader.hpp"

#include "CLI11/CLI11.hpp"
#include <iostream>

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;
using axom::inlet::DocWriter;


std::shared_ptr<Inlet> createBasicInlet(DataStore* ds,
                                        const std::string& luaString)
{
  auto lr = std::make_shared<LuaReader>();
  lr->parseString(luaString);

  return std::make_shared<Inlet>(lr, ds->getRoot());
}

int main(int argc, char** argv) {
  // CLI::App app;
  // CLI11_PARSE(app, argc, argv);
  // bool docs_enabled;
  // app.add_flag("-f", docs_enabled, "Enables documentation generation");
  
  std::string testString = "foo = true; bar = false";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);
  
  std::shared_ptr<axom::inlet::Field> currField;

  currField = inlet->addBool("foo", "foo's description");
  assert(currField);
  currField->required(true);
  currField = inlet->addBool("bar", "bar's description");
  assert(currField);
  currField->required(false);

  DocWriter doc("example_docs.rst", inlet->sidreGroup()); 

  return 0;
}