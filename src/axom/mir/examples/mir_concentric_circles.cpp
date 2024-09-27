// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

#include <conduit.hpp>
#include <conduit_relay_io_blueprint.hpp>

#include <string>

// namespace aliases
namespace mir = axom::mir;
namespace bputils = axom::mir::utilities::blueprint;

using RuntimePolicy = axom::runtime_policy::Policy;

//--------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  SLIC_ERROR(
    axom::fmt::format("Error from Conduit: s1={}, s2={}, i1={}", s1, s2, i1));
  // This is on purpose.
  //while(1)
  //  ;
}

//--------------------------------------------------------------------------------
void printNode(const conduit::Node &n)
{
  conduit::Node options;
  options["num_children_threshold"] = 10000;
  options["num_elements_threshold"] = 10000;
  n.to_summary_string_stream(std::cout, options);
}

//--------------------------------------------------------------------------------
template <typename ExecSpace>
int runMIR(const conduit::Node &hostMesh,
           const conduit::Node &options,
           conduit::Node &hostResult)
{
  using namespace axom::mir::views;
  SLIC_INFO(axom::fmt::format("Using policy {}",
                              axom::execution_space<ExecSpace>::name()));

  // Check materials.
  constexpr int MAXMATERIALS = 20;
  auto materialInfo = materials(hostMesh["matsets/mat"]);
  if(materialInfo.size() >= MAXMATERIALS)
  {
    SLIC_WARNING(
      axom::fmt::format("To use more than {} materials, recompile with "
                        "larger MAXMATERIALS value.",
                        MAXMATERIALS));
    return -4;
  }

  // host->device
  conduit::Node deviceMesh;
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);

  AXOM_ANNOTATE_BEGIN("runMIR");
  // Make views (we know beforehand which types to make)
  using CoordsetView = ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(
    bputils::make_array_view<float>(deviceMesh["coordsets/coords/values/x"]),
    bputils::make_array_view<float>(deviceMesh["coordsets/coords/values/y"]));

  using TopoView = UnstructuredTopologySingleShapeView<QuadShape<int>>;
  TopoView topoView(bputils::make_array_view<int>(
    deviceMesh["topologies/mesh/elements/connectivity"]));

  using MatsetView = UnibufferMaterialView<int, float, MAXMATERIALS>;
  MatsetView matsetView;
  matsetView.set(
    bputils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
    bputils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
    bputils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));

  using MIR =
    axom::mir::EquiZAlgorithm<ExecSpace, TopoView, CoordsetView, MatsetView>;
  MIR m(topoView, coordsetView, matsetView);
  conduit::Node deviceResult;
  m.execute(deviceMesh, options, deviceResult);

  AXOM_ANNOTATE_END("runMIR");

  // device->host
  bputils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);

  return 0;
}

//--------------------------------------------------------------------------------
int runMIR(RuntimePolicy policy,
           int gridSize,
           int numCircles,
           const std::string &outputFilePath)
{
  // Initialize a mesh for testing MIR
  auto timer = axom::utilities::Timer(true);
  mir::MeshTester tester;
  conduit::Node mesh;
  {
    AXOM_ANNOTATE_SCOPE("generate");
    tester.initTestCaseFive(gridSize, numCircles, mesh);
    // printNode(mesh);
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Output initial mesh.
  {
    AXOM_ANNOTATE_SCOPE("save_input");
    conduit::relay::io::blueprint::save_mesh(mesh, "concentric_circles", "hdf5");
  }

  // Begin material interface reconstruction
  timer.start();
  conduit::Node options, resultMesh;
  options["matset"] = "mat";

  int retval = 0;
std::cout << axom::fmt::format("policy={}", policy) << std::endl;
  if(policy == RuntimePolicy::seq)
  {
#pragma message "SEQ supported"
    retval = runMIR<axom::SEQ_EXEC>(mesh, options, resultMesh);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
  else if(policy == RuntimePolicy::omp)
  {
#pragma message "OMP supported"
    retval = runMIR<axom::OMP_EXEC>(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(policy == RuntimePolicy::cuda)
  {
#pragma message "CUDA supported"
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
    retval = runMIR<cuda_exec>(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(policy == RuntimePolicy::hip)
  {
#pragma message "HIP supported"
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
    retval = runMIR<hip_exec>(mesh, options, resultMesh);
  }
  #endif
#endif
  else
  {
    retval = -1;
    SLIC_ERROR("Unhandled policy.");
  }
  timer.stop();
  SLIC_INFO("Material interface reconstruction time: "
            << timer.elapsedTimeInMilliSec() << " ms.");

  // Output results
  {
    AXOM_ANNOTATE_SCOPE("save_output");
    conduit::relay::io::blueprint::save_mesh(resultMesh, outputFilePath, "hdf5");
  }

  return retval;
}

//--------------------------------------------------------------------------------
int runMIROld(int gridSize, int numCircles, const std::string &outputFilePath)
{
  // Initialize a mesh for testing MIR
  auto timer = axom::utilities::Timer(true);
  mir::MeshTester tester;
  mir::MIRMesh mesh;
  {
    AXOM_ANNOTATE_SCOPE("generate");
    mesh = tester.initTestCaseFive(gridSize, numCircles);
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Begin material interface reconstruction
  timer.start();
  mir::MIRMesh outputMesh;
  {
    AXOM_ANNOTATE_SCOPE("runMIR");
    mir::InterfaceReconstructor m;
    m.computeReconstructedInterface(mesh, outputMesh);
  }
  timer.stop();
  SLIC_INFO("Material interface reconstruction time: "
            << timer.elapsedTimeInMilliSec() << " ms.");

  // Output results
  {
    AXOM_ANNOTATE_SCOPE("save_output");
    outputMesh.writeMeshToFile(".", outputFilePath + ".vtk");
  }

  return 0;
}

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  axom::slic::SimpleLogger logger;  // create & initialize test logger
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  // Define command line options.
  bool handler = true, old = false;
  int gridSize = 5;
  int numCircles = 2;
  std::string outputFilePath("output");
  axom::CLI::App app;
  app.add_flag("--handler", handler)
    ->description("Install a custom error handler that loops forever.")
    ->capture_default_str();
  app.add_flag("--old", old)
    ->description("Use the old MIR method.")
    ->capture_default_str();
  app.add_option("--gridsize", gridSize)
    ->description("The number of zones along an axis.");
  app.add_option("--numcircles", numCircles)
    ->description("The number of circles to use for material creation.");
  app.add_option("--output", outputFilePath)
    ->description("The file path for output files");

#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("report");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  RuntimePolicy policy {RuntimePolicy::seq};
  std::stringstream pol_sstr;
  pol_sstr << "Set runtime policy for intersection-based sampling method.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
  pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
  pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
  pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
  #endif
#endif
  app.add_option("-p, --policy", policy, pol_sstr.str())
    ->capture_default_str()
    ->transform(
      axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

  // Parse command line options.
  app.parse(argc, argv);

  if(handler)
  {
    conduit::utils::set_error_handler(conduit_debug_err_handler);
  }
#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    annotationMode);
#endif

  int retval = 0;
  try
  {
    if(old)
    {
      retval = runMIROld(gridSize, numCircles, outputFilePath);
    }
    else
    {
      retval = runMIR(policy, gridSize, numCircles, outputFilePath);
    }
  }
  catch(std::invalid_argument const &e)
  {
    SLIC_WARNING("Bad input. " << e.what());
    retval = -2;
  }
  catch(std::out_of_range const &e)
  {
    SLIC_WARNING("Integer overflow. " << e.what());
    retval = -3;
  }
  return retval;
}