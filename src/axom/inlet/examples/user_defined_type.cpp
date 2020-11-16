// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include <unordered_map>
#include "CLI11/CLI11.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

namespace inlet = axom::inlet;

struct Mesh
{
  std::string filename;
  int serial_ref_iter;
  int par_ref_iter;

  // Each class should define a static method that adds the fields it
  // will grab from inlet
  static void defineSchema(inlet::Table& schema)
  {
    schema.addString("filename", "Path to mesh file");
    schema.addInt("serial", "Number of serial refinement iterations");
    schema.addInt("parallel", "Number of parallel refinement iterations");
  }
};

// Additionally, each class should specialize this struct as follows
// in the global namespace so that Inlet can access it
/**
 * Example Lua definition:
 * \code{.lua}
 * mesh = {
 *    filename = 'data/square.mesh',
 *    serial = 12,
 *    parallel = 7
 * }
 * \endcode
 */
template <>
struct FromInlet<Mesh>
{
  Mesh operator()(const inlet::Table& base)
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
  static void defineSchema(inlet::Table& schema)
  {
    schema.addDouble("rel_tol", "Relative convergence criterion");
    schema.addDouble("abs_tol", "Relative convergence criterion");
    schema.addInt("print_level", "Logging level for iterative solver");
    schema.addInt("max_iter", "Maximum iteration count");
    schema.addDouble("dt", "Time step");
    schema.addInt("steps", "Number of simulation iterations/frames");
  }
};

/**
 * Example Lua definition:
 * \code{.lua}
 * solver = {
 *    rel_tol = 1.e-6,
 *    abs_tol = 1.e-12,
 *    print_level = 0,
 *    max_iter = 100,
 *    dt = 1.0,
 *    steps = 1 
 * }
 * \endcode
 */
template <>
struct FromInlet<LinearSolver>
{
  LinearSolver operator()(const inlet::Table& base)
  {
    LinearSolver lin_solve;
    lin_solve.rel_tol = base["rel_tol"];
    lin_solve.abs_tol = base["abs_tol"];
    lin_solve.print_level = base["print_level"];
    lin_solve.max_iter = base["max_iter"];
    lin_solve.dt = base["dt"];
    lin_solve.steps = base["steps"];
    return lin_solve;
  }
};

struct BoundaryCondition
{
  std::unordered_map<int, int> attrs;
  double constant;
  static void defineSchema(inlet::Table& schema)
  {
    schema.addIntArray("attrs", "List of boundary attributes");
    schema.addDouble("constant",
                     "The scalar to fix the value of the solution to");
  }
};

/**
 * Example Lua definition:
 * \code{.lua}
 * bc = {
 *   attrs = {
 *      3, 4, 6, 9
 *   }
 *   constant = 12.55
 * }
 * \endcode
 */
template <>
struct FromInlet<BoundaryCondition>
{
  BoundaryCondition operator()(const inlet::Table& base)
  {
    BoundaryCondition bc;
    bc.attrs = base["attrs"];
    bc.constant = base["constant"];
    return bc;
  }
};

struct ThermalSolver
{
  Mesh mesh;
  LinearSolver solver;
  std::unordered_map<std::string, BoundaryCondition> bcs;
  // defineSchema is intended to be used recursively
  // Tables are created for subobjects and passed to
  // subobject defineSchema implementations
  static void defineSchema(inlet::Table& schema)
  {
    auto& mesh_table = schema.addTable("mesh", "Information about the mesh");
    Mesh::defineSchema(mesh_table);
    auto& solver_table =
      schema.addTable("solver",
                      "Information about the iterative solver used for Ku = f");
    LinearSolver::defineSchema(solver_table);

    // Schema only needs to be defined once, will propagate through to each
    // element of the array, namely, the subtable at each found index in the input file
    auto& bc_table =
      schema.addGenericDict("bcs", "List of boundary conditions");
    BoundaryCondition::defineSchema(bc_table);
  }
};

/**
 * Example Lua definition:
 * \code{.lua}
 * thermal_solver = {
 *    mesh = {
 *      -- see above FromInlet<Mesh>
 *    },
 *    solver = {
 *      -- see above FromInlet<LinearSolver>
 *    },
 *    bcs = {
 *        ["temperature"] = {
 *          -- see above FromInlet<BoundaryCondition>
 *        },
 *        ["flux"] = {
 *          -- see above FromInlet<BoundaryCondition>
 *        },
 *    }
 * }
 * \endcode
 */
template <>
struct FromInlet<ThermalSolver>
{
  // This is also implicitly recursive - will call the FromInlet
  // functions defined for the subobjects
  ThermalSolver operator()(const inlet::Table& base)
  {
    return {base["mesh"].get<Mesh>(),
            base["solver"].get<LinearSolver>(),
            base["bcs"].get<std::unordered_map<std::string, BoundaryCondition>>()};
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
  auto lr = std::make_unique<LuaReader>();
  lr->parseFile(inputFileName);
  Inlet inlet(std::move(lr), ds.getRoot());

  // Create a table off the global table for the thermal_solver object
  // then define its schema
  auto& thermal_solver_table =
    inlet.addTable("thermal_solver",
                   "Configuration for a thermal conduction module");
  ThermalSolver::defineSchema(thermal_solver_table);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Read all the data into a thermal solver object
  auto thermal_solver = inlet["thermal_solver"].get<ThermalSolver>();
}
