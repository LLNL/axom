// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

#include <conduit.hpp>
#include <conduit_relay_io_blueprint.hpp>

#include <string>

// namespace aliases
namespace mir = axom::mir;
namespace bputils = axom::mir::utilities::blueprint;

//--------------------------------------------------------------------------------

std::string usageString()
{
  return "Args are <grid size> <number of circles> <output file path>";
}

//--------------------------------------------------------------------------------

void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  //  SLIC_ERROR(axom::fmt::format("Error from Conduit: s1={}, s2={}, i1={}", s1, s2, i1));
  // This is on purpose.
  while(1)
    ;
}

//--------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  axom::slic::SimpleLogger logger;  // create & initialize test logger
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);
  conduit::utils::set_error_handler(conduit_debug_err_handler);
#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper("report");
#endif

  if(argc != 4)
  {
    SLIC_WARNING("Incorrect number of args. " << usageString());
    return -1;
  }

  int retval = 0;
  try
  {
    // Parse the command line arguments
    int gridSize = std::stoi(argv[1]);
    int numCircles = std::stoi(argv[2]);
    std::string outputFilePath = std::string(argv[3]);

    // Initialize a mesh for testing MIR
    auto timer = axom::utilities::Timer(true);
    mir::MeshTester tester;
    conduit::Node mesh;
    tester.initTestCaseFive(gridSize, numCircles, mesh);
    mesh.print();
    timer.stop();
    SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");
    // Output initial mesh.
    conduit::relay::io::blueprint::save_mesh(mesh, "concentric_circles", "hdf5");

    // Make views (we know beforehand which types to make)
    using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 2>;
    CoordsetView coordsetView(
      bputils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
      bputils::make_array_view<float>(mesh["coordsets/coords/values/y"]));

    using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<
      axom::mir::views::QuadShape<int>>;
    TopoView topoView(bputils::make_array_view<int>(
      mesh["topologies/mesh/elements/connectivity"]));

    constexpr int MAXMATERIALS = 20;
    auto materialInfo = axom::mir::views::materials(mesh["matsets/mat"]);
    if(materialInfo.size() >= MAXMATERIALS)
    {
      SLIC_WARNING(
        axom::fmt::format("To use more than {} materials, recompile with "
                          "larger MAXMATERIALS value.",
                          MAXMATERIALS));
      return -4;
    }
    using MatsetView =
      axom::mir::views::UnibufferMaterialView<int, float, MAXMATERIALS>;
    MatsetView matsetView;
    matsetView.set(
      bputils::make_array_view<int>(mesh["matsets/mat/material_ids"]),
      bputils::make_array_view<float>(mesh["matsets/mat/volume_fractions"]),
      bputils::make_array_view<int>(mesh["matsets/mat/sizes"]),
      bputils::make_array_view<int>(mesh["matsets/mat/offsets"]),
      bputils::make_array_view<int>(mesh["matsets/mat/indices"]));

    // Begin material interface reconstruction
    using MIR =
      axom::mir::EquiZAlgorithm<axom::SEQ_EXEC, TopoView, CoordsetView, MatsetView>;
    timer.start();
    conduit::Node options, processedMesh;
    MIR m(topoView, coordsetView, matsetView);
    options["matset"] = "mat";
    m.execute(mesh, options, processedMesh);
    timer.stop();
    SLIC_INFO("Material interface reconstruction time: "
              << timer.elapsedTimeInMilliSec() << " ms.");

    // Output results
    conduit::relay::io::blueprint::save_mesh(processedMesh,
                                             outputFilePath,
                                             "hdf5");

    retval = 0;
  }
  catch(std::invalid_argument const &e)
  {
    SLIC_WARNING("Bad input. " << usageString());
    retval = -2;
  }
  catch(std::out_of_range const &e)
  {
    SLIC_WARNING("Integer overflow. " << usageString());
    retval = -3;
  }
  return retval;
}
