// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// usage : ./document_generation_example --enableDocs --deck lua_file 

#include "axom/inlet.hpp"

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

void defineSchema(std::shared_ptr<Inlet> inlet)
{
  std::shared_ptr<axom::inlet::Field> currField;

  // Add the description to the thermal_solver/mesh/filename Field
  currField = inlet->addString("thermal_solver/mesh/filename", "file for thermal solver");
  // Set the field's required property to true
  currField->required(true);
  
  currField = inlet->addInt("thermal_solver/mesh/serial", "serial value");  
  // The property called required is left unspecified here

  // The description for thermal_solver/mesh/parallel is left unspecified
  currField = inlet->addInt("thermal_solver/mesh/parallel");
  currField->required(false);

  currField = inlet->addInt("thermal_solver/order", "thermal solver order");
  currField->required(true);
  
  currField = inlet->addString("thermal_solver/timestepper", "thermal solver timestepper");
  currField->required(false);
  currField->defaultValue("this is default");

  currField = inlet->addString("thermal_solver/u0/type", "description for u0 type");
  currField->required(true);

  currField = inlet->addString("thermal_solver/u0/func", "description for u0 func"); 
  currField->required(true);
  
  currField = inlet->addString("thermal_solver/kappa/type", "description for kappa type");
  currField->required(true);

  currField = inlet->addDouble("thermal_solver/kappa/constant", "description for kappa constant");
  currField->required(true);
  currField->defaultValue(0.0);

  // Add description to solver table by using the addTable function
  auto table = inlet->addTable("thermal_solver/solver", "This is the solver sub-table in the thermal_solver table");

  // You can also add fields through a table

  currField = table->addDouble("rel_tol", "description for solver rel tol");
  currField->required(false);
  currField->range(0.5, 100.7);
  
  currField = table->addDouble("abs_tol", "description for solver abs tol");
  currField->required(true);  

  currField = table->addInt("print_level", "description for solver print level");
  currField->required(true); 
  currField->validValues({1, 3, 5, 7});

  currField = table->addInt("max_iter", "description for solver max iter");
  currField->required(false);
  currField->defaultValue(10);
  
  currField = table->addDouble("dt", "description for solver dt");
  currField->required(true); 

  currField = table->addInt("steps", "description for solver steps");
  currField->required(true);
}

// Checking the contents of the passed inlet
void checkValues(std::shared_ptr<Inlet> inlet) {
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
}

int main(int argc, char** argv) {
  CLI::App app {"Basic example of Axom's Inlet component"};
  bool docsEnabled{false};
  app.add_flag("--enableDocs", docsEnabled, "Enables documentation generation");

  std::string inputFileName;
  auto opt = app.add_option("--deck", inputFileName, "Path to input deck file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  // Create inlet and parse input file data into the inlet

  DataStore ds;
  auto lr = std::make_shared<LuaReader>();
  lr->parseFile(inputFileName);
  auto inlet = std::make_shared<Inlet>(lr, ds.getRoot(), docsEnabled);
  auto docWriter = std::make_shared<SphinxDocWriter>("example_doc.rst", inlet->sidreGroup());
  inlet->registerDocWriter(docWriter);

  defineSchema(inlet);
  checkValues(inlet);
  
  // Generate the documentation
  inlet->writeDoc();

  if (docsEnabled) {
    std::cout << "Documentation was written to example_doc.rst" << std::endl;
  }

  return 0;
}
