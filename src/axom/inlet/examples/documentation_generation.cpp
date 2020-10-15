// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// usage : ./inlet_documentation_generation_example --enableDocs --fil lua_file.lua

#include "axom/inlet.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

#include "CLI11/CLI11.hpp"
#include <iostream>

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::inlet::SphinxDocWriter;
using axom::sidre::DataStore;

void findStr(std::string path, std::shared_ptr<Inlet> inlet)
{
  auto proxy = (*inlet)[path];
  if(proxy.type() == axom::inlet::InletType::String)
  {
    std::cout << "found " << proxy.get<std::string>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void findInt(std::string path, std::shared_ptr<Inlet> inlet)
{
  auto proxy = (*inlet)[path];
  if(proxy.type() == axom::inlet::InletType::Integer)
  {
    std::cout << "found " << proxy.get<int>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void findDouble(std::string path, std::shared_ptr<Inlet> inlet)
{
  auto proxy = (*inlet)[path];
  if(proxy.type() == axom::inlet::InletType::Double)
  {
    std::cout << "found " << proxy.get<double>();
  }
  else
  {
    std::cout << "not found ";
  }
  std::cout << std::endl;
}

void defineSchema(std::shared_ptr<Inlet> inlet)
{
  // Add the description to the thermal_solver/mesh/filename Field
  auto currVerifiable =
    inlet->addString("thermal_solver/mesh/filename", "file for thermal solver");
  // Set the field's required property to true
  currVerifiable->required(true);

  currVerifiable = inlet->addInt("thermal_solver/mesh/serial", "serial value");
  currVerifiable->range(0, INT_MAX);
  currVerifiable->defaultValue(1);

  // The description for thermal_solver/mesh/parallel is left unspecified
  currVerifiable = inlet->addInt("thermal_solver/mesh/parallel");
  currVerifiable->range(1, INT_MAX);
  currVerifiable->defaultValue(1);

  currVerifiable = inlet->addInt("thermal_solver/order", "thermal solver order");
  currVerifiable->required(true);
  currVerifiable->range(1, INT_MAX);

  currVerifiable =
    inlet->addString("thermal_solver/timestepper", "thermal solver timestepper");
  currVerifiable->defaultValue("quasistatic");
  currVerifiable->validValues({"quasistatic", "forwardeuler", "backwardeuler"});

  currVerifiable =
    inlet->addString("thermal_solver/u0/type", "description for u0 type");
  currVerifiable->defaultValue("constant");
  currVerifiable->validValues({"constant", "function"});

  currVerifiable =
    inlet->addString("thermal_solver/u0/func", "description for u0 func");
  currVerifiable->required(true);

  currVerifiable =
    inlet->addString("thermal_solver/kappa/type", "description for kappa type");
  currVerifiable->required(true);
  currVerifiable->validValues({"constant", "function"});

  currVerifiable = inlet->addDouble("thermal_solver/kappa/constant",
                                    "description for kappa constant");
  currVerifiable->required(true);

  // Add description to solver table by using the addTable function
  auto table =
    inlet->addTable("thermal_solver/solver",
                    "This is the solver sub-table in the thermal_solver table");

  // You can also add fields through a table

  currVerifiable = table->addDouble("rel_tol", "description for solver rel tol");
  currVerifiable->required(false);
  currVerifiable->defaultValue(1.e-6);
  currVerifiable->range(0.0, __DBL_MAX__);

  currVerifiable = table->addDouble("abs_tol", "description for solver abs tol");
  currVerifiable->required(true);
  currVerifiable->defaultValue(1.e-12);
  currVerifiable->range(0.0, __DBL_MAX__);

  currVerifiable =
    table->addInt("print_level", "description for solver print level");
  currVerifiable->required(true);
  currVerifiable->defaultValue(0);
  currVerifiable->range(0, 3);

  currVerifiable = table->addInt("max_iter", "description for solver max iter");
  currVerifiable->required(false);
  currVerifiable->defaultValue(100);
  currVerifiable->range(1, INT_MAX);

  currVerifiable = table->addDouble("dt", "description for solver dt");
  currVerifiable->required(true);
  currVerifiable->defaultValue(1);
  currVerifiable->range(0.0, __DBL_MAX__);

  currVerifiable = table->addInt("steps", "description for solver steps");
  currVerifiable->required(true);
  currVerifiable->defaultValue(1);
  currVerifiable->range(1, INT_MAX);
}

// Checking the contents of the passed inlet
void checkValues(std::shared_ptr<Inlet> inlet)
{
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

  // Verify that contents of Inlet meet the requirements of the specified schema
  if(inlet->verify())
  {
    SLIC_INFO("Inlet verify successful.");
  }
  else
  {
    SLIC_INFO("Inlet verify failed.");
  }
}

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  // This is a generic basic SLIC logger
  axom::slic::UnitTestLogger logger;

  CLI::App app {"Basic example of Axom's Inlet component"};
  bool docsEnabled {false};
  app.add_flag("--enableDocs", docsEnabled, "Enables documentation generation");

  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  // Create inlet and parse input file data into the inlet

  DataStore ds;
  auto lr = std::make_shared<LuaReader>();
  lr->parseFile(inputFileName);
  auto inlet = std::make_shared<Inlet>(lr, ds.getRoot(), docsEnabled);

  // _inlet_documentation_generation_start
  auto docWriter =
    std::make_shared<SphinxDocWriter>("example_doc.rst", inlet->sidreGroup());
  inlet->registerDocWriter(docWriter);
  // _inlet_documentation_generation_end

  defineSchema(inlet);
  checkValues(inlet);

  // Generate the documentation
  inlet->writeDoc();

  if(docsEnabled)
  {
    SLIC_INFO("Documentation was written to example_doc.rst\n");
  }

  return 0;
}
