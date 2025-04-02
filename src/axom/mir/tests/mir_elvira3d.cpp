// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

//------------------------------------------------------------------------------

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression", "mir_elvira3d");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_Elvira3D
{
  static const int gridSize = 10;
  static const int numSpheres = 2;

  static void initialize(conduit::Node &n_mesh)
  {
    AXOM_ANNOTATE_SCOPE("initialize");
    axom::mir::MeshTester M;
    M.setStructured(true);
    M.initTestCaseSix(gridSize, numSpheres, n_mesh);
  }

  // Select a chunk of zones.
  static int selectZones(conduit::Node &n_options)
  {
    std::vector<axom::IndexType> selected;
    for(int k = 0; k < gridSize; k++)
    {
      for(int j = 0; j < gridSize; j++)
      {
        for(int i = 0; i < gridSize; i++)
        {
          // Save all but an octant so we can see inside.
          if(i > gridSize / 3 && j > gridSize / 3 && k > gridSize / 3)
          {
            continue;
          }
          auto idx = k * gridSize * gridSize + j * gridSize + i;
          selected.push_back(idx);
        }
      }
    }
    n_options["selectedZones"].set(selected);
    return static_cast<int>(selected.size());
  }

  static void test(const std::string &AXOM_UNUSED_PARAM(name), bool selectedZones = false)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    namespace views = axom::mir::views;

    const double expectedVolume = gridSize * gridSize * gridSize;
    double mirExpectedVolume = expectedVolume;

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(hostMesh);
    {
      AXOM_ANNOTATE_SCOPE("host_to_device");
      axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
    }
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    {
      AXOM_ANNOTATE_SCOPE("save_original");
  #if defined(AXOM_USE_HDF5)
      conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  #endif
      conduit::relay::io::save(hostMesh, name + "_orig.yaml", "yaml");
    }
#endif

    //--------------------------------------------------------------------------
    const conduit::Node &n_coordset =
      deviceMesh.fetch_existing("coordsets/coords");
    const conduit::Node &n_topology =
      deviceMesh.fetch_existing("topologies/mesh");
    const conduit::Node &n_matset = deviceMesh.fetch_existing("matsets/mat");

    // Make views.
    auto coordsetView = views::make_explicit_coordset<float, 3>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

    auto topologyView = views::make_structured<3>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView = views::make_unibuffer_matset<int, float, 3>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Do MIR
    using MIR =
      axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node deviceMIRMesh;
    conduit::Node options;
    options["matset"] = "mat";
    // Be more lenient in how away far points are in order to combine them.
    options["point_tolerance"] = 1.e-4;
    if(selectedZones)
    {
      // The MIR volume changes when the set of selected zones changes.
      mirExpectedVolume = selectZones(options);
    }
    m.execute(deviceMesh, options, deviceMIRMesh);

    //--------------------------------------------------------------------------
    // Compute volumes for original mesh as a field.
    AXOM_ANNOTATE_BEGIN("volume");
    bputils::MakeZoneVolumes<ExecSpace, TopologyView, CoordsetView> origZV(
      topologyView,
      coordsetView);
    origZV.execute(n_topology, n_coordset, deviceMesh["fields/volume"]);

    //--------------------------------------------------------------------------
    // Compute volumes for MIR mesh as a field.
    conduit::Node &n_mir_coordset = deviceMIRMesh["coordsets/coords"];
    auto mirCoordsetView =
      views::make_explicit_coordset<float, 3>::view(n_mir_coordset);
    using MirCoordsetView = decltype(mirCoordsetView);

    // Make polyhedral topology view.
    using MirTopologyView = views::UnstructuredTopologyPolyhedralView<int>;
    const conduit::Node &n_mir_topology = deviceMIRMesh["topologies/mesh"];
    MirTopologyView mirTopoView(
      bputils::make_array_view<int>(n_mir_topology["subelements/connectivity"]),
      bputils::make_array_view<int>(n_mir_topology["subelements/sizes"]),
      bputils::make_array_view<int>(n_mir_topology["subelements/offsets"]),
      bputils::make_array_view<int>(n_mir_topology["elements/connectivity"]),
      bputils::make_array_view<int>(n_mir_topology["elements/sizes"]),
      bputils::make_array_view<int>(n_mir_topology["elements/offsets"]));

    const conduit::Node &n_mir_matset = deviceMIRMesh["matsets/mat"];
    auto mirMatsetView =
      views::make_unibuffer_matset<int, float, 3>::view(n_mir_matset);

    bputils::MakeZoneVolumes<ExecSpace, MirTopologyView, MirCoordsetView> mirZV(
      mirTopoView,
      mirCoordsetView);
    mirZV.execute(n_mir_topology, n_mir_coordset, deviceMIRMesh["fields/volume"]);
    AXOM_ANNOTATE_END("volume");

    //--------------------------------------------------------------------------
    // device->host
    conduit::Node hostMIRMesh;
    {
      AXOM_ANNOTATE_SCOPE("device_to_host");
      bputils::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);
    }

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    {
      AXOM_ANNOTATE_SCOPE("save_baseline");
  #if defined(AXOM_USE_HDF5)
      conduit::relay::io::blueprint::save_mesh(hostMIRMesh, name, "hdf5");
  #endif
      conduit::relay::io::save(hostMIRMesh, name + ".yaml", "yaml");
    }
#endif

    // Handle baseline comparison.
    // NOTE: Comparing against a baseline is turned off for now because on various
    //       backends, we do not get quite the same answer for new coordinates in
    //       clipped zones. This ultimately results in different coordset and
    //       connectivity values. Rather than compare the Conduit nodes, we skip
    //       it for now and compare material volumes before/after MIR.
#if 0
    {
      AXOM_ANNOTATE_SCOPE("compare_baseline");
      std::string baselineName(yamlRoot(name));
      const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
      saveBaseline(paths, baselineName, hostMIRMesh);
#else
      constexpr double tolerance = 2.6e-06;
      EXPECT_TRUE(compareBaseline(paths, baselineName, hostMIRMesh, tolerance));
#endif
    }
#endif
    //--------------------------------------------------------------------------
    // Compute the total volumes on the original and MIR meshes.
    constexpr double tolerance = 1.5e-5;

    const auto orig_volume =
      bputils::make_array_view<double>(deviceMesh["fields/volume/values"]);
    EXPECT_NEAR(expectedVolume, variableSum(orig_volume), tolerance);
    const auto origMatInfo = views::materials(n_matset);
    const auto origTotalVolumes =
      sumMaterialVolumes(matsetView, orig_volume, origMatInfo);

    const auto mir_volume =
      bputils::make_array_view<double>(deviceMIRMesh["fields/volume/values"]);
    EXPECT_NEAR(mirExpectedVolume, variableSum(mir_volume), tolerance);
    const auto mirMatInfo = views::materials(n_mir_matset);
    const auto mirTotalVolumes =
      sumMaterialVolumes(mirMatsetView, mir_volume, mirMatInfo);

    //--------------------------------------------------------------------------
    // comparisons
    EXPECT_EQ(origTotalVolumes.size(), mirTotalVolumes.size());
    double volumeSums[2] = {0., 0.};
    for(size_t i = 0; i < origTotalVolumes.size(); i++)
    {
      if(!selectedZones)
      {
        SLIC_INFO(axom::fmt::format("Material {}: origVF = {}, mirVF = {}",
                                    i,
                                    origTotalVolumes[i],
                                    mirTotalVolumes[i]));
      }
      volumeSums[0] += origTotalVolumes[i];
      volumeSums[1] += mirTotalVolumes[i];
    }
    EXPECT_NEAR(volumeSums[0], expectedVolume, tolerance);
    EXPECT_NEAR(volumeSums[1], mirExpectedVolume, tolerance);
    if(!selectedZones)
    {
      for(size_t i = 0; i < origTotalVolumes.size(); i++)
      {
        EXPECT_NEAR(origTotalVolumes[i], mirTotalVolumes[i], tolerance);
      }
    }
  }

  static double variableSum(axom::ArrayView<double> var)
  {
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    RAJA::ReduceSum<reduce_policy, double> reduceVar(0.);
    axom::for_all<ExecSpace>(
      var.size(),
      AXOM_LAMBDA(axom::IndexType i) { reduceVar += var[i]; });
    return reduceVar.get();
  }

  /*!
   * \brief Compute the total volumes for the materials in the mesh.
   *
   * \param matsetView The matset view that contains the material data (on device)
   * \param zoneVolumes An array view that contains zone volumes (on device)
   * \param matInfo The material information for the mesh.
   *
   * \return A vector of volumes on the host, 1 element per material.
   */
  template <typename MatsetView>
  static std::vector<double> sumMaterialVolumes(
    MatsetView matsetView,
    axom::ArrayView<double> zoneVolumes,
    const axom::mir::views::MaterialInformation &matInfo)
  {
    using atomic_policy =
      typename axom::execution_space<ExecSpace>::atomic_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    AXOM_ANNOTATE_SCOPE("sumMaterialVolumes");

    // Make a sorted list of material numbers on device.
    const int nmats = static_cast<int>(matInfo.size());
    axom::Array<int> sortedIdsHost(matInfo.size());
    int mi = 0;
    for(const auto &mat : matInfo)
    {
      sortedIdsHost[mi++] = mat.number;
    }
    axom::utilities::Sorting<int>::sort(sortedIdsHost.data(),
                                        sortedIdsHost.size());
    axom::Array<int> sortedIds(nmats, nmats, allocatorID);
    axom::copy(sortedIds.data(), sortedIdsHost.data(), nmats * sizeof(int));
    auto sortedIdsView = sortedIds.view();

    // Compute the total volumes for each material.
    axom::Array<double> totalVolume(nmats, nmats, allocatorID);
    auto totalVolumeView = totalVolume.view();
    axom::for_all<ExecSpace>(
      nmats,
      AXOM_LAMBDA(axom::IndexType i) { totalVolumeView[i] = 0.; });
    axom::for_all<ExecSpace>(
      matsetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zi) {
        // Get the materials for this zone.
        typename MatsetView::IDList ids;
        typename MatsetView::VFList vfs;
        matsetView.zoneMaterials(zi, ids, vfs);

        // Add the material volumes to the total volumes.
        for(axom::IndexType i = 0; i < ids.size(); i++)
        {
          auto index = axom::mir::utilities::bsearch(ids[i], sortedIdsView);
          SLIC_ASSERT(index >= 0 && index < nmats);

          // Use an atomic to sum the value.
          RAJA::atomicAdd<atomic_policy>(totalVolumeView.data() + index,
                                         zoneVolumes[zi] * vfs[i]);
        }
      });

    // Move results back to host.
    std::vector<double> hostVolumes(nmats, 0.);
    axom::copy(hostVolumes.data(), totalVolume.data(), nmats * sizeof(double));

    return hostVolumes;
  }
};

//------------------------------------------------------------------------------
TEST(mir_elvira3d, elvira3d_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_seq");
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_seq");
  constexpr bool selectZones = true;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer_sel", selectZones);
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_elvira3d, elvira3d_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_omp");
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_omp");
  constexpr bool selectZones = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel", selectZones);
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_elvira3d, elvira3d_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_cuda");
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_cuda");
  constexpr bool selectZones = true;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer_sel", selectZones);
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_elvira3d, elvira3d_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_hip");
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_hip");
  constexpr bool selectZones = true;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer_sel", selectZones);
}
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  // Define command line options.
  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  bool handlerEnabled = false;
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");

  // Parse command line options.
  try
  {
    app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
    axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
      annotationMode);
#endif

    axom::slic::SimpleLogger logger;  // create & initialize test logger,
    if(handlerEnabled)
    {
      conduit::utils::set_error_handler(conduit_debug_err_handler);
    }

    result = RUN_ALL_TESTS();
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    result = 0;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    result = app.exit(e);
  }

  return result;
}
