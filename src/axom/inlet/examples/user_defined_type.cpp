// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include "CLI11/CLI11.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::inlet::SchemaCreator;
using axom::sidre::DataStore;

struct Mesh
{
  std::string filename;
  int serial_ref_iter;
  int par_ref_iter;

  // Each class should define a static method that adds the fields it
  // will grab from inlet
  static void defineSchema(SchemaCreator& schema)
  {
    schema.addString("filename", "Path to mesh file");
    schema.addInt("serial", "Number of serial refinement iterations");
    schema.addInt("parallel", "Number of parallel refinement iterations");
  }
};

// Additionally, each class should specialize this struct as follows
// in the global namespace so that Inlet can access it
// Alternatively, if a class is default-constructible, an equivalent
// by-reference function can be implemented (see LinearSolver below)
template <>
struct FromInlet<Mesh>
{
  Mesh operator()(axom::inlet::Table& base)
  {
    return {base["filename"], base["serial"], base["parallel"]};
  }
};

struct LinearSolver
{
  double rel_tol;
  double abs_tol;
  int print_level;
  int max_iter;
  double dt;
  int steps;
  static void defineSchema(SchemaCreator& schema)
  {
    schema.addDouble("rel_tol", "Relative convergence criterion");
    schema.addDouble("abs_tol", "Relative convergence criterion");
    schema.addInt("print_level", "Logging level for iterative solver");
    schema.addInt("max_iter", "Maximum iteration count");
    schema.addDouble("dt", "Time step");
    schema.addInt("steps", "Number of simulation iterations/frames");
  }
};

// Example definition of a "deserializer" using a free function
// that takes a reference param
void from_inlet(axom::inlet::Table& base, LinearSolver& lin_solve)
{
  lin_solve.rel_tol = base["rel_tol"];
  lin_solve.abs_tol = base["abs_tol"];
  lin_solve.print_level = base["print_level"];
  lin_solve.max_iter = base["max_iter"];
  lin_solve.dt = base["dt"];
  lin_solve.steps = base["steps"];
}

struct ThermalSolver
{
  Mesh mesh;
  LinearSolver solver;
  // defineSchema is intended to be used recursively
  // Tables are created for subobjects and passed to
  // subobject defineSchema implementations
  static void defineSchema(SchemaCreator& schema)
  {
    auto mesh_table = schema.addTable("mesh", "Information about the mesh");
    Mesh::defineSchema(*mesh_table);
    auto solver_table =
      schema.addTable("solver",
                      "Information about the iterative solver used for Ku = f");
    LinearSolver::defineSchema(*solver_table);
  }
};

template <>
struct FromInlet<ThermalSolver>
{
  // This is also implicitly recursive - will call the from_inlet
  // functions defined for the subobjects
  ThermalSolver operator()(axom::inlet::Table& base)
  {
    return {base["mesh"], base["solver"]};
  }
};

int main(int argc, char** argv)
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::UnitTestLogger logger;

  CLI::App app {"Example of Axom's Inlet component with user-defined types"};
  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  auto lr = std::make_shared<LuaReader>();
  lr->parseFile(inputFileName);
  auto inlet = std::make_shared<Inlet>(lr, ds.getRoot());

  // Create a table off the global table for the thermal_solver object
  // then define its schema
  auto thermal_solver_table =
    inlet->addTable("thermal_solver",
                    "Configuration for a thermal conduction module");
  ThermalSolver::defineSchema(*thermal_solver_table);

  if(!inlet->verify())
  {
    SLIC_INFO("Inlet failed to verify against provided schema");
  }

  // Read all the data into a thermal solver object
  ThermalSolver thermal_solver = (*inlet)["thermal_solver"];
}
