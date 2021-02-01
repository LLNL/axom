// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include <unordered_map>
#include "axom/slic/core/SimpleLogger.hpp"

using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::LuaReader;
using axom::sidre::DataStore;

namespace inlet = axom::inlet;
// _inlet_userdef_simple_start
struct Mesh
{
  std::string filename;
  int serial_ref_iter;
  int par_ref_iter;

  // A function can be used to define a schema for a particular struct
  // For convenience, it can be implemented as a static method of the struct
  static void defineSchema(inlet::Table& schema)
  {
    schema.addString("filename", "Path to mesh file");
    schema.addInt("serial", "Number of serial refinement iterations");
    schema.addInt("parallel", "Number of parallel refinement iterations");
  }
};
// _inlet_userdef_simple_end

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
// _inlet_userdef_simple_frominlet_start
template <>
struct FromInlet<Mesh>
{
  Mesh operator()(const inlet::Table& base)
  {
    return {base["filename"], base["serial"], base["parallel"]};
  }
};
// _inlet_userdef_simple_frominlet_end

const std::string input = R"(
mesh = {
  filename = 'data/square.mesh',
  serial = 12,
  parallel = 7
}
)";

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  DataStore ds;
  auto lr = std::make_unique<LuaReader>();
  lr->parseString(input);
  Inlet inlet(std::move(lr), ds.getRoot());

  // Create a table off the global table for the mesh object
  // then define its schema
  auto& mesh_schema =
    inlet.addStruct("mesh", "Information used to read in/process a mesh");
  Mesh::defineSchema(mesh_schema);

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Read all the data into a thermal solver object
  Mesh mesh = inlet["mesh"].get<Mesh>();
  SLIC_INFO(fmt::format("Mesh has filename '{0}'", mesh.filename));
  SLIC_INFO(fmt::format("Mesh has {0} serial refinement iterations",
                        mesh.serial_ref_iter));
  SLIC_INFO(fmt::format("Mesh has {0} parallel refinement iterations",
                        mesh.par_ref_iter));
}
