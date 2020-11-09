// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include <unordered_map>
#include "CLI11/CLI11.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"
#include "mfem.hpp"

using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

namespace inlet = axom::inlet;

// This example relies on mfem's support for std::function-based coefficients
#if MFEM_VERSION_MAJOR >= 4 && MFEM_VERSION_MINOR >= 2
  #define MFEM_STDFUNCTION_COEF
#endif

#ifdef MFEM_STDFUNCTION_COEF

struct BoundaryCondition
{
  std::unordered_map<int, int> attrs;
  // These act as a sort of variant without C++17
  std::unique_ptr<mfem::Coefficient> coef;
  std::unique_ptr<mfem::VectorCoefficient> vec_coef;

  static void defineSchema(inlet::Table& schema)
  {
    schema.addIntArray("attrs", "List of boundary attributes");
    // Inlet does not support sum types, so both options are added to the schema
    // Supported function parameter/return types are Double, Vec2D, and Vec3D
    // Only single-argument functions are currently supported
    schema.addFunction("vec_coef",
                       inlet::InletFunctionType::Vec3D,  // Return type
                       inlet::InletFunctionType::Vec3D,  // Argument type
                       "The function representing the BC coefficient");

    schema.addFunction("coef",
                       inlet::InletFunctionType::Double,  // Return type
                       inlet::InletFunctionType::Vec3D,   // Argument type
                       "The function representing the BC coefficient");
  }
};

/**
 * Example Lua definition:
 * \code{.lua}
 * bc = {
 *   attrs = {
 *      3, 4, 6, 9
 *   }
 *   coef = function (x, y, z)
 *     return x * 0.12
 *   end
 * }
 * -- or, for vector coefficients:
 * [8] = {
 *   attrs = { [4] = 14, [8] = 62, [6] = 11},
 *   vec_coef = function (x, y, z)
 *     scale = 0.12
 *     return x * scale, y * scale, z * scale
 *   end
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
    if(base.contains("vec_coef"))
    {
      // We assume the dimension is 3
      bc.vec_coef = std::make_unique<mfem::VectorFunctionCoefficient>(
        3,
        [&base](const mfem::Vector& input, mfem::Vector& output) {
          // One way of accessing the function is with "call" - need
          // to specify the desired return and argument types explicitly
          auto ret = base["vec_coef"].call<axom::primal::Vector3D>(
            axom::primal::Vector3D {input.GetData()});
          // Copy from the primal vector into the MFEM vector
          std::copy(ret.data(), ret.data() + ret.dimension(), output.GetData());
        });
      SLIC_INFO("Created an mfem::VectorCoefficient with dimension "
                << bc.vec_coef->GetVDim());
    }
    else
    {
      // Another way of accessing the function is by extracting the std::function
      auto func =
        base["coef"].get<std::function<double(axom::primal::Vector3D)>>();
      bc.coef = std::make_unique<mfem::FunctionCoefficient>(
        [func(std::move(func))](const mfem::Vector& vec) {
          // func is a concrete function type so calls can leverage conversions
          return func(vec.GetData());
        });
      SLIC_INFO("Created an mfem::Coefficient");
    }
    return bc;
  }
};

#endif  // MFEM_STDFUNCTION_COEF

int main(int argc, char** argv)
{
#ifdef MFEM_STDFUNCTION_COEF
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

  // We only need the boundary condition sub-table
  auto& bc_table =
    inlet.addGenericArray("thermal_solver/bcs", "List of boundary conditions");
  BoundaryCondition::defineSchema(bc_table);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Read all the data into a thermal solver object
  auto bcs =
    inlet["thermal_solver/bcs"].get<std::unordered_map<int, BoundaryCondition>>();

#endif  // MFEM_STDFUNCTION_COEF
}
