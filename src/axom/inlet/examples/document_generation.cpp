// usage : ./document_generation_example --enableDocs --deck lua_file 

#include "axom/inlet/SphinxDocWriter.hpp"
#include "axom/inlet/Inlet.hpp"
#include "axom/inlet/LuaReader.hpp"

#include "CLI11/CLI11.hpp"
#include <iostream>

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;
using axom::inlet::SphinxDocWriter;

void findStr(std::string path, std::shared_ptr<Inlet> inlet) {
  std::string strVal;
  bool found = inlet->get(path, strVal);
  std::cout << path << ": ";
  found ? std::cout << "found " << strVal : std::cout << "not found ";
  std::cout << std::endl;
} 

void findInt(std::string path, std::shared_ptr<Inlet> inlet) {
  int val;
  bool found = inlet->get(path, val);
  std::cout << path << ": ";
  found ? std::cout << "found " << val : std::cout << "not found ";
  std::cout << std::endl;
} 

void findDouble(std::string path, std::shared_ptr<Inlet> inlet) {
  double val;
  bool found = inlet->get(path, val);
  std::cout << path << ": ";
  found ? std::cout << "found " << val : std::cout << "not found ";
  std::cout << std::endl;
} 

void exampleInlet(DataStore* ds, const std::string& luaFile, bool docsEnabled)
{
  auto lr = std::make_shared<LuaReader>();
  lr->parseFile(luaFile);

  auto inlet = std::make_shared<Inlet>(lr, ds->getRoot(), docsEnabled);

  auto docWriter = std::make_shared<SphinxDocWriter>("example_docs.rst", inlet->sidreGroup());
  inlet->registerDocWriter(docWriter);
  
  std::shared_ptr<axom::inlet::Field> currField;

  currField = inlet->addString("thermal_solver/mesh/filename", "file for thermal solver");
  currField->required(true);
  
  currField = inlet->addInt("thermal_solver/mesh/serial", "serial value");  
  currField->required(false);
  
  currField = inlet->addInt("thermal_solver/mesh/parallel", "parallel value");
  currField->required(true);
  
  currField = inlet->addInt("thermal_solver/order", "thermal solver order");
  currField->required(true);
  
  currField = inlet->addString("thermal_solver/timestepper", "thermal solver timestepper");
  currField->required(true);
  
  currField = inlet->addString("thermal_solver/u0/type", "description for u0 type");
  currField->required(true);

  currField = inlet->addString("thermal_solver/u0/func", "description for u0 func"); 
  currField->required(true);
  
  currField = inlet->addString("thermal_solver/kappa/type", "description for kappa type");
  currField->required(true);

  currField = inlet->addDouble("thermal_solver/kappa/constant", "description for kappa constant");
  currField->required(true);
  
  currField = inlet->addDouble("thermal_solver/solver/rel_tol", "description for solver rel tol");
  currField->required(true);
  
  currField = inlet->addDouble("thermal_solver/solver/abs_tol", "description for solver abs tol");
  currField->required(true);  

  currField = inlet->addInt("thermal_solver/solver/print_level", "description for solver print level");
  currField->required(true);  

  currField = inlet->addInt("thermal_solver/solver/max_iter", "description for solver max iter");
  currField->required(true);  
  
  currField = inlet->addDouble("thermal_solver/solver/dt", "description for solver dt");
  currField->required(true); 

  currField = inlet->addInt("thermal_solver/solver/steps", "description for solver steps");
  currField->required(true);  

  // Checking contents of inlet

  findStr("thermal_solver/mesh/filename", inlet);
  findStr("thermal_solver/timestepper", inlet);
  findStr("thermal_solver/u0/type", inlet);
  findStr("thermal_solver/u0/func", inlet);
  findStr("thermal_solver/kappa/type", inlet);

  findInt("thermal_solver/mesh/serial", inlet);
  findInt("thermal_solver/mesh/parallel", inlet);
  findInt("thermal_solver/order", inlet);
  findInt("thermal_solver/solver/print_level", inlet);
  findInt("thermal_solver/solver/max_iter", inlet);
  findInt("thermal_solver/solver/steps", inlet);

  findDouble("thermal_solver/solver/dt", inlet);
  findDouble("thermal_solver/solver/abs_tol", inlet);
  findDouble("thermal_solver/solver/rel_tol", inlet);
  findDouble("thermal_solver/kappa/constant", inlet);
  
  inlet->writeDocs();
}

int main(int argc, char** argv) {
  CLI::App app {"Description here"};
  bool docsEnabled{false};
  app.add_flag("--enableDocs", docsEnabled, "Enables documentation generation");

  std::string strOption;
  auto opt = app.add_option("--deck", strOption, "Path to input deck file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  exampleInlet(&ds, strOption, docsEnabled);

  if (docsEnabled) {
    std::cout << "Documentation was written to example.rst" << std::endl;
  }

  return 0;
}