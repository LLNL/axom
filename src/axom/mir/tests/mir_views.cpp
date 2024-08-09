// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <cmath>

// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #else
    using hip_exec = seq_exec;
  #endif
#endif
// clang-format on

//------------------------------------------------------------------------------

// Uncomment to generate baselines
#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
#define AXOM_TESTING_SAVE_VISUALIZATION

// Include after seq_exec is defined.
#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "mir"), "regression"),
               "mir_views");
}
//------------------------------------------------------------------------------

TEST(mir_views, shape2conduitName)
{
  EXPECT_EQ(axom::mir::views::LineShape<int>::name(), "line");
  EXPECT_EQ(axom::mir::views::LineShape<long>::name(), "line");

  EXPECT_EQ(axom::mir::views::TriShape<int>::name(), "tri");
  EXPECT_EQ(axom::mir::views::TriShape<long>::name(), "tri");

  EXPECT_EQ(axom::mir::views::QuadShape<int>::name(), "quad");
  EXPECT_EQ(axom::mir::views::QuadShape<long>::name(), "quad");

  EXPECT_EQ(axom::mir::views::TetShape<int>::name(), "tet");
  EXPECT_EQ(axom::mir::views::TetShape<long>::name(), "tet");

  EXPECT_EQ(axom::mir::views::PyramidShape<int>::name(), "pyramid");
  EXPECT_EQ(axom::mir::views::PyramidShape<long>::name(), "pyramid");

  EXPECT_EQ(axom::mir::views::WedgeShape<int>::name(), "wedge");
  EXPECT_EQ(axom::mir::views::WedgeShape<long>::name(), "wedge");

  EXPECT_EQ(axom::mir::views::HexShape<int>::name(), "hex");
  EXPECT_EQ(axom::mir::views::HexShape<long>::name(), "hex");
}

TEST(mir_views, explicit_coordsetview)
{
  axom::Array<float> x {{0., 1., 2., 3., 4., 5.}};
  axom::Array<float> y {{10., 11., 12., 13., 14., 15.}};
  axom::Array<float> z {{20., 21., 22., 23., 24., 25.}};

  axom::mir::views::ExplicitCoordsetView<float, 2> view2d(x.view(), y.view());
  EXPECT_EQ(view2d.size(), 6);
  for(axom::IndexType i = 0; i < view2d.size(); i++)
  {
    axom::primal::Point<float, 2> P({x[i], y[i]});
    EXPECT_EQ(view2d.getPoint(i), P);
    EXPECT_EQ(view2d[i], P);
  }

  axom::mir::views::ExplicitCoordsetView<float, 3> view3d(x.view(),
                                                          y.view(),
                                                          z.view());
  EXPECT_EQ(view3d.size(), 6);
  for(axom::IndexType i = 0; i < view3d.size(); i++)
  {
    axom::primal::Point<float, 3> P({x[i], y[i], z[i]});
    EXPECT_EQ(view3d.getPoint(i), P);
    EXPECT_EQ(view3d[i], P);
  }
}

TEST(mir_views, strided_structured)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::strided_structured<2>(hostMesh);
  //  hostMesh.print();

  // These are the expected zone ids for this strided structured mesh.
  // clang-format off
  const axom::Array<int> expectedZones {{16, 17, 24, 23,
                                         17, 18, 25, 24,
                                         18, 19, 26, 25,
                                         23, 24, 31, 30,
                                         24, 25, 32, 31,
                                         25, 26, 33, 32}};
  // clang-format on
  auto expectedZonesView = expectedZones.view();

  axom::mir::views::dispatch_explicit_coordset(
    hostMesh["coordsets/coords"],
    [&](auto coordsetView) {
      axom::mir::views::dispatch_structured_topology<
        axom::mir::views::select_dimensions(2)>(
        hostMesh["topologies/mesh"],
        [&](const std::string &shape, auto topoView) {
          // Traverse the zones in the mesh and check the zone ids.
          topoView.template for_all_zones<axom::SEQ_EXEC>(
            AXOM_LAMBDA(auto zoneIndex, const auto &zone) {
              // Check zone ids.
              const auto ids = zone.getIds();
              for(axom::IndexType i = 0; i < ids.size(); i++)
              {
                EXPECT_EQ(expectedZonesView[zoneIndex * 4 + i], ids[i]);
              }

              // Check coordinates
              const auto nodeIndexing = topoView.indexing().expand();
              for(axom::IndexType i = 0; i < ids.size(); i++)
              {
                // Get coordinate from coordsetView.
                const auto pt = coordsetView[ids[i]];

                // Get the logical local id for the id.
                const auto index = nodeIndexing.GlobalToLocal(ids[i]);
                const auto logical = nodeIndexing.IndexToLogicalIndex(index);

                // Expected coordinate
                double x = (3. + 1. / 3.) * static_cast<double>(logical[0] - 1);
                const double yvals[] = {-2, 2, 6};
                double y = yvals[logical[1]];

                const double dx = pt[0] - x;
                const double dy = pt[1] - y;
                double d = sqrt(dx * dx + dy * dy);

                EXPECT_TRUE(d < 1.e-10);
              }
            });
        });
    });
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_mat_test(const std::string &type, const std::string &mattype, const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  axom::mir::testing::data::make_matset(mattype, "mesh", zoneDims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig_yaml", "yaml");
#endif

  using MatsetView = axom::mir::views::UnibufferMaterialView<int, float, 3>;
  MatsetView matsetView;
  matsetView.set(axom::mir::utilities::blueprint::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                 axom::mir::utilities::blueprint::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                 axom::mir::utilities::blueprint::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                 axom::mir::utilities::blueprint::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                 axom::mir::utilities::blueprint::make_array_view<int>(deviceMesh["matsets/mat/indices"]));

  constexpr int MATA = 0;
  constexpr int MATB = 1;
  constexpr int MATC = 2;
  const int matids[] = {0, 36, 40};
  // clang-format off
  const int results[] = {/*contains mat*/ 0, 1, 0, /*mats in zone*/ 1, /*mats in zone*/ MATB, -1, -1,
                         /*contains mat*/ 1, 1, 0, /*mats in zone*/ 2, /*mats in zone*/ MATA, MATB, -1,
                         /*contains mat*/ 1, 1, 1, /*mats in zone*/ 3, /*mats in zone*/ MATA, MATB, MATC};
  // clang-format on
  constexpr int nZones = sizeof(matids) / sizeof(int);

  // Get matids into matidsView for device.
  axom::Array<int> matidsArray(nZones, nZones, allocatorID);
  axom::copy(matidsArray.data(), matids, sizeof(int) * nZones);
  auto matidsView = matidsArray.view();

  // Allocate results array on device.
  constexpr int nResults = sizeof(results) / sizeof(int);
  axom::Array<int> resultsArrayDevice(nResults, nResults, allocatorID);
  auto resultsView = resultsArrayDevice.view();

  // host-safe tests.
  EXPECT_EQ(matsetView.numberOfZones(), zoneDims[0] * zoneDims[1]);
  EXPECT_EQ(matsetView.numberOfMaterials(), 3);

  // Fill in resultsView on the device.
  constexpr int nResultsPerZone = nResults / nZones;
  axom::for_all<ExecSpace>(3, AXOM_LAMBDA(auto index)
  {
    resultsView[nResultsPerZone * index + 0] = matsetView.zoneContainsMaterial(matidsView[index], MATA) ? 1 : 0;
    resultsView[nResultsPerZone * index + 1] = matsetView.zoneContainsMaterial(matidsView[index], MATB) ? 1 : 0;
    resultsView[nResultsPerZone * index + 2] = matsetView.zoneContainsMaterial(matidsView[index], MATC) ? 1 : 0;
    resultsView[nResultsPerZone * index + 3] = matsetView.numberOfMaterials(matidsView[index]);

    MatsetView::IDList ids;
    MatsetView::VFList vfs;
    matsetView.zoneMaterials(matidsView[index], ids, vfs);
    for(int i = 0; i < 3; i++)
    {
      resultsView[nResultsPerZone * index + 4 + i] = (i < ids.size()) ? ids[i] : -1;
    }
  });
  // Get containsView data to the host and compare results
  std::vector<int> resultsHost(nResults);
  axom::copy(resultsHost.data(), resultsView.data(), sizeof(int) * nResults);
  for(int i = 0; i < nResults; i++)
  {
    EXPECT_EQ(results[i], resultsHost[i]);
  }
}

TEST(mir_views, matset_unibuffer)
{
  braid2d_mat_test<seq_exec>("uniform", "unibuffer", "uniform2d_unibuffer");

#if defined(AXOM_USE_OPENMP)
  braid2d_mat_test<omp_exec>("uniform", "unibuffer", "uniform2d_unibuffer");
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid2d_mat_test<cuda_exec>("uniform", "unibuffer", "uniform2d_unibuffer");
#endif

#if defined(AXOM_USE_HIP)
  braid2d_mat_test<hip_exec>("uniform", "unibuffer", "uniform2d_unibuffer");
#endif
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
