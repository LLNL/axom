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
    // Inlet does not support sum types, so options are added to the schema
    // for vector/scalar coefficients and for time-dependent versions of each
    // Supported function parameter/return types are Double and Vec3D
    schema.addFunction("vec_coef",
                       inlet::FunctionType::Vec3D,    // Return type
                       {inlet::FunctionType::Vec3D},  // Argument type
                       "The function representing the BC coefficient");

    schema.addFunction("coef",
                       inlet::FunctionType::Double,
                       {inlet::FunctionType::Vec3D},
                       "The function representing the BC coefficient");

    schema.addFunction("vec_coef_t",
                       inlet::FunctionType::Vec3D,
                       {inlet::FunctionType::Vec3D,
                        inlet::FunctionType::Double},  // Multiple argument types
                       "The function representing the BC coefficient");

    schema.addFunction("coef_t",
                       inlet::FunctionType::Double,
                       {inlet::FunctionType::Vec3D, inlet::FunctionType::Double},
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
      // We assume the dimension is 3 - if the mfem::Vector is 2D,
      // the z-component will be set to zero due to how primal::Vector's
      // (T*, size) constructor works
      bc.vec_coef = std::make_unique<mfem::VectorFunctionCoefficient>(
        3,
        [&base](const mfem::Vector& input, mfem::Vector& output) {
          // One way of accessing the function is with "call" - need
          // to specify the desired return and argument types explicitly
          auto ret = base["vec_coef"].call<axom::primal::Vector3D>(
            axom::primal::Vector3D {input.GetData(), input.Size()});
          // Copy from the primal vector into the MFEM vector
          std::copy(ret.data(), ret.data() + ret.dimension(), output.GetData());
        });
      SLIC_INFO("Created an mfem::VectorCoefficient with dimension "
                << bc.vec_coef->GetVDim());
    }
    else if(base.contains("coef"))
    {
      // Another way of accessing the function is by extracting the std::function
      auto func =
        base["coef"].get<std::function<double(axom::primal::Vector3D)>>();
      bc.coef = std::make_unique<mfem::FunctionCoefficient>(
        [func(std::move(func))](const mfem::Vector& vec) {
          // func is a concrete function type so calls can leverage conversions
          return func({vec.GetData(), vec.Size()});
        });
      SLIC_INFO("Created an mfem::Coefficient");
    }
    else if(base.contains("coef_t"))
    {
      auto func =
        base["coef_t"].get<std::function<double(axom::primal::Vector3D, double)>>();
      bc.coef = std::make_unique<mfem::FunctionCoefficient>(
        [func(std::move(func))](const mfem::Vector& vec, double t) {
          return func({vec.GetData(), vec.Size()}, t);
        });
      SLIC_INFO("Created a time-dependent mfem::Coefficient");
    }
    else if(base.contains("vec_coef_t"))
    {
      auto func =
        base["vec_coef_t"]
          .get<std::function<axom::primal::Vector3D(axom::primal::Vector3D, double)>>();
      bc.vec_coef = std::make_unique<mfem::VectorFunctionCoefficient>(
        3,
        [func(std::move(
          func))](const mfem::Vector& input, double t, mfem::Vector& output) {
          auto ret = func({input.GetData(), input.Size()}, t);
          // Copy from the primal vector into the MFEM vector
          std::copy(ret.data(), ret.data() + ret.dimension(), output.GetData());
        });
      SLIC_INFO(
        "Created a time-dependent mfem::VectorCoefficient with dimension "
        << bc.vec_coef->GetVDim());
    }
    else
    {
      SLIC_ERROR("Table did not contain a coefficient function: " << base.name());
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
  // Intended to be used with mfem_coef.lua
  std::string inputFileName;
  auto opt = app.add_option("--file", inputFileName, "Path to input file");
  opt->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  DataStore ds;
  auto lr = std::make_unique<LuaReader>();
  lr->parseFile(inputFileName);
  Inlet inlet(std::move(lr), ds.getRoot());

  // We only need the boundary condition sub-table
  auto& bc_table = inlet.addGenericArray("bcs", "List of boundary conditions");
  BoundaryCondition::defineSchema(bc_table);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Read all the data into a thermal solver object
  auto bcs = inlet["bcs"].get<std::unordered_map<int, BoundaryCondition>>();

#endif  // MFEM_STDFUNCTION_COEF
}
