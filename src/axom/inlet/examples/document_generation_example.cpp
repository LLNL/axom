#include "axom/inlet/DocWriter.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/LuaReader.hpp"

#include "CLI11/CLI11.hpp"
#include <iostream>

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

int main(int argc, char** argv) {
  CLI::App app;
  CLI11_PARSE(app, argc, argv);
  bool docs_enabled;
  app.add_flag("-genDocs", docs_enabled, "Enables documentation generation");
  
  std::string testString = "foo = true; bar = false";
  DataStore ds;
  auto inlet = createBasicInlet(&ds, testString);

  //
  // Define schema
  //

  std::shared_ptr<axom::inlet::Field> currField;

  // Check for existing fields
  currField = inlet->addBool("foo", "foo's description");
  assert(currField);
  currField->required(true);
  currField = inlet->addBool("bar", "bar's description");
  assert(currField);
  currField->required(false);

  return 0;
}