// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

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
               "mir_clipfield");
}
//------------------------------------------------------------------------------

TEST(mir_clipfield, options)
{
  int nzones = 6;

  conduit::Node options;
  axom::mir::clipping::ClipOptions<seq_exec> opts(nzones, options);

  options["clipField"] = "distance";
  EXPECT_EQ(opts.clipField(), options["clipField"].as_string());

  EXPECT_EQ(opts.clipValue(), 0.);
  options["clipValue"] = 2.5f;
  EXPECT_EQ(opts.clipValue(), 2.5f);

  EXPECT_EQ(opts.topologyName("default"), "default");
  options["topologyName"] = "topo";
  EXPECT_EQ(opts.topologyName("default"), "topo");

  EXPECT_EQ(opts.coordsetName("default"), "default");
  options["coordsetName"] = "coords";
  EXPECT_EQ(opts.coordsetName("default"), "coords");

  EXPECT_EQ(opts.colorField(), "color");
  options["colorField"] = "custom_color";
  EXPECT_EQ(opts.colorField(), "custom_color");

  EXPECT_TRUE(opts.inside());
  options["inside"] = 1;
  EXPECT_TRUE(opts.inside());
  options["inside"] = 0;
  EXPECT_FALSE(opts.inside());

  EXPECT_FALSE(opts.outside());
  options["outside"] = 1;
  EXPECT_TRUE(opts.outside());
  options["outside"] = 0;
  EXPECT_FALSE(opts.outside());

  // The clip field has to be present
  conduit::Node n_fields;
  n_fields["distance/topology"] = "topo";
  n_fields["distance/association"] = "vertex";
  n_fields["distance/values"].set(std::vector<float> {0., 1., 2., 3.});

  // There are currently no fields in the options. fields should just return the clip field.
  auto fields = opts.fields(n_fields);
  EXPECT_EQ(fields.size(), 1);
  EXPECT_EQ(fields.begin()->first, "distance");
  EXPECT_EQ(fields.begin()->second, "distance");

  // Add an empty fields node so we select NO fields.
  (void)options["fields"];
  fields = opts.fields(n_fields);
  EXPECT_EQ(fields.size(), 0);

  // Add some fields
  options["fields/distance"] = "distance";
  options["fields/source"] = "destination";
  options["fields/same"] = 1;
  fields = opts.fields(n_fields);
  EXPECT_EQ(fields.size(), 3);
  int i = 0;
  for(auto it = fields.begin(); it != fields.end(); it++, i++)
  {
    if(i == 0)
    {
      EXPECT_EQ(it->first, "distance");
      EXPECT_EQ(it->second, "distance");
    }
    else if(i == 1)
    {
      EXPECT_EQ(it->first, "same");
      EXPECT_EQ(it->second, "same");
    }
    else if(i == 2)
    {
      EXPECT_EQ(it->first, "source");
      EXPECT_EQ(it->second, "destination");
    }
  }

  // There are no "selectedZones" in the options. We should get nzones values from 0 onward.
  auto selectedZonesView = opts.selectedZonesView();
  EXPECT_EQ(selectedZonesView.size(), 6);
  EXPECT_EQ(selectedZonesView[0], 0);
  EXPECT_EQ(selectedZonesView[1], 1);
  EXPECT_EQ(selectedZonesView[2], 2);
  EXPECT_EQ(selectedZonesView[3], 3);
  EXPECT_EQ(selectedZonesView[4], 4);
  EXPECT_EQ(selectedZonesView[5], 5);

  // Put some "selectedZones" in the options.
  opts.invalidateSelectedZones();
  options["selectedZones"].set(std::vector<axom::IndexType> {5, 4, 3});
  selectedZonesView = opts.selectedZonesView();
  EXPECT_EQ(selectedZonesView.size(), 3);
  EXPECT_EQ(selectedZonesView[0], 3);
  EXPECT_EQ(selectedZonesView[1], 4);
  EXPECT_EQ(selectedZonesView[2], 5);
}

TEST(mir_clipfield, blend_group_builder)
{
  using IndexType = axom::IndexType;
  using KeyType = std::uint64_t;

  /*

  We'll make 2 quads

  3      4      5
  *--8---*------*
  |      |      |
  |  6   9      |
  |      |      |
  *--7---*------*
  0      1      2

  */
  axom::Array<IndexType> blendGroups {{8, 5}};
  axom::Array<IndexType> blendGroupsLen {
    {/*zone 0*/ 4 + 1 + 1 + 1 + 1 + 2 + 2 + 2, /*zone 1*/ 1 + 1 + 1 + 1 + 2}};
  axom::Array<IndexType> blendGroupOffsets {{0, 8}};
  axom::Array<IndexType> blendOffsets {{0, blendGroupsLen[0]}};

  axom::Array<KeyType> blendNames {
    {/*zone 0*/ 6, 0, 1, 3, 4, 7, 8, 9, /*zone 1*/ 1, 2, 4, 5, 9}};
  axom::Array<IndexType> blendGroupSizes {
    {/*zone 0*/ 4, 1, 1, 1, 1, 2, 2, 2, /*zone 1*/ 1, 1, 1, 1, 2}};
  axom::Array<IndexType> blendGroupStart {
    {/*zone 0*/ 0, 4, 5, 6, 7, 8, 10, 12, /*zone 1*/ 13, 14, 15, 16, 18}};
  axom::Array<IndexType> blendIds {{
    /*zone 0*/
    0,
    1,
    2,
    3,  // 6 (bgname) // 0 (bgindex)
    0,  // 0          // 1
    1,  // 1          // 2
    3,  // 3          // 3
    4,  // 4          // 4
    0,
    1,  // 7          // 5
    3,
    4,  // 8          // 6
    1,
    4,  // 9          // 7
    /*zone 1*/
    1,  // 1          // 8
    2,  // 2          // 9
    4,  // 4          // 10
    5,  // 5          // 11
    1,
    4  // 9          // 12
  }};
  axom::Array<float> blendCoeff {{/*zone 0*/
                                  0.25,
                                  0.25,
                                  0.25,
                                  0.25,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  0.5,
                                  /*zone 1*/
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  0.5,
                                  0.5}};
  axom::Array<KeyType> blendUniqueNames {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};
  axom::Array<KeyType> blendUniqueIndices {{1, 2, 9, 3, 4, 11, 0, 5, 6, 7}};

  axom::mir::clipping::BlendGroupBuilder<seq_exec> builder;
  builder.setBlendGroupSizes(blendGroups.view(), blendGroupsLen.view());
  builder.setBlendGroupOffsets(blendOffsets.view(), blendGroupOffsets.view());
  builder.setBlendViews(blendNames.view(),
                        blendGroupSizes.view(),
                        blendGroupStart.view(),
                        blendIds.view(),
                        blendCoeff.view());

  //std::cout << "-------- zone 0 --------" << std::endl;
  auto z0 = builder.blendGroupsForZone(0);
  EXPECT_EQ(z0.size(), 8);
  IndexType index = 0;
  for(IndexType i = 0; i < z0.size(); i++, index++)
  {
    //z0.print(std::cout);
    EXPECT_EQ(z0.ids().size(), blendGroupSizes[index]);

    z0++;
  }

  //std::cout << "-------- zone 1 --------" << std::endl;
  auto z1 = builder.blendGroupsForZone(1);
  EXPECT_EQ(z1.size(), 5);
  for(IndexType i = 0; i < z1.size(); i++, index++)
  {
    //z1.print(std::cout);
    EXPECT_EQ(z1.ids().size(), blendGroupSizes[index]);
    z1++;
  }
}

//------------------------------------------------------------------------------
template <typename ArrayType>
bool increasing(const ArrayType &arr)
{
  bool retval = true;
  for(size_t i = 1; i < arr.size(); i++) retval &= (arr[i] >= arr[i - 1]);
  return retval;
}

std::vector<int> permute(const std::vector<int> &input)
{
  std::vector<int> values, indices;
  std::vector<double> order;

  values.resize(input.size());
  indices.resize(input.size());
  order.resize(input.size());

  std::iota(indices.begin(), indices.end(), 0);
  for(size_t i = 0; i < input.size(); i++)
  {
    order[i] = drand48();
  }
  std::sort(indices.begin(), indices.end(), [&](int a, int b) {
    return order[a] < order[b];
  });
  for(size_t i = 0; i < input.size(); i++) values[i] = input[indices[i]];
  return values;
}

std::vector<int> makeUnsortedArray(int n)
{
  std::vector<int> values;
  values.resize(n);
  std::iota(values.begin(), values.end(), 0);
  return permute(values);
}

std::vector<int> makeRandomArray(int n)
{
  constexpr int largestId = 1 << 28;
  std::vector<int> values;
  values.resize(n);
  for(int i = 0; i < n; i++)
  {
    values[i] = static_cast<int>(largestId * drand48());
  }
  return permute(values);
}

//------------------------------------------------------------------------------
TEST(mir_clipfield, sort_values)
{
  for(int n = 1; n < 15; n++)
  {
    for(int trial = 1; trial <= n; trial++)
    {
      auto values = makeUnsortedArray(n);
      axom::mir::utilities::sort_values(values.data(), values.size());
      EXPECT_TRUE(increasing(values));
    }
  }
}
//------------------------------------------------------------------------------
TEST(mir_clipfield, unique)
{
  /*
  8---9---10--11
  |   |   |   |
  4---5---6---7
  |   |   |   |
  0---1---2---3

 */
  axom::Array<int> ids {{0, 1, 5, 4, 1, 2, 6,  5, 2, 3, 7,  6,
                         4, 5, 9, 8, 5, 6, 10, 9, 6, 7, 11, 10}};
  axom::Array<int> uIds, uIndices;

  axom::mir::utilities::unique<seq_exec>(ids.view(), uIds, uIndices);
  EXPECT_EQ(uIds.size(), 12);
  EXPECT_EQ(uIndices.size(), 12);
  for(axom::IndexType i = 0; i < uIds.size(); i++)
  {
    EXPECT_EQ(uIds[i], i);
    EXPECT_EQ(uIds[i], ids[uIndices[i]]);
  }
}

//------------------------------------------------------------------------------
TEST(mir_clipfield, make_name)
{
  for(int n = 1; n < 14; n++)
  {
    // Make a set of scrambled ids.
    auto values = makeRandomArray(n);
    std::uint64_t name = 0, name2 = 0;
    // Compute the name for that list of ids.
    if(n == 1)
      name = axom::mir::utilities::make_name_1(values[0]);
    else if(n == 2)
      name = axom::mir::utilities::make_name_2(values[0], values[1]);
    else
      name = axom::mir::utilities::make_name_n(values.data(), values.size());

    for(int trial = 0; trial < 1000; trial++)
    {
      // Scramble the id list.
      auto values2 = permute(values);
      // Compute the name for that list of ids.
      if(n == 1)
        name2 = axom::mir::utilities::make_name_1(values2[0]);
      else if(n == 2)
        name2 = axom::mir::utilities::make_name_2(values2[0], values2[1]);
      else
        name2 = axom::mir::utilities::make_name_n(values2.data(), values2.size());

      // The names for the 2 scrambled lists of numbers should be the same.
      EXPECT_EQ(name, name2);
    }
  }
}

template <typename ExecSpace, typename ShapeType>
void test_one_shape(const conduit::Node &hostMesh, const std::string &name)
{
  using TopoView =
    axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<float, 3>;

  // Copy mesh to device
  conduit::Node deviceMesh;
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);

  // Make views for the device mesh.
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  axom::ArrayView<float> xView(static_cast<float *>(n_x.data_ptr()),
                               n_x.dtype().number_of_elements());
  axom::ArrayView<float> yView(static_cast<float *>(n_y.data_ptr()),
                               n_y.dtype().number_of_elements());
  axom::ArrayView<float> zView(static_cast<float *>(n_z.data_ptr()),
                               n_z.dtype().number_of_elements());
  CoordsetView coordsetView(xView, yView, zView);

  conduit::Node &n_conn =
    deviceMesh.fetch_existing("topologies/topo/elements/connectivity");
  axom::ArrayView<int> connView(static_cast<int *>(n_conn.data_ptr()),
                                n_conn.dtype().number_of_elements());
  TopoView topoView(connView);

  // Clip the data
  conduit::Node deviceClipMesh, options;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
    topoView,
    coordsetView);
  options["clipField"] = "distance";
  options["clipValue"] = 0.;
  options["inside"] = 1;
  options["outside"] = 1;
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

template <typename ShapeType>
void test_one_shape_exec(const conduit::Node &hostMesh, const std::string &name)
{
  test_one_shape<seq_exec, ShapeType>(hostMesh, name);

#if defined(AXOM_USE_OPENMP)
  test_one_shape<omp_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_one_shape<cuda_exec, ShapeType>(hostMesh, name);
#endif

#if defined(AXOM_USE_HIP)
  test_one_shape<hip_exec, ShapeType>(hostMesh, name);
#endif
}

TEST(mir_clipfield, onetet)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_tet(hostMesh);
  test_one_shape_exec<axom::mir::views::TetShape<int>>(hostMesh, "one_tet");
}

TEST(mir_clipfield, onepyr)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_pyr(hostMesh);
  test_one_shape_exec<axom::mir::views::PyramidShape<int>>(hostMesh, "one_pyr");
}

TEST(mir_clipfield, onewdg)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_wdg(hostMesh);
  test_one_shape_exec<axom::mir::views::WedgeShape<int>>(hostMesh, "one_wdg");
}

TEST(mir_clipfield, onehex)
{
  conduit::Node hostMesh;
  axom::mir::testing::data::make_one_hex(hostMesh);
  test_one_shape_exec<axom::mir::views::HexShape<int>>(hostMesh, "one_hex");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_clip_test(const std::string &type, const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Create views
  axom::StackArray<double, 2> origin {0., 0.}, spacing {1., 1.};
  CoordsetView coordsetView(dims, origin, spacing);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  const std::string clipTopoName("cliptopo");
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;
  options["topologyName"] = clipTopoName;
  options["coordsetName"] = "clipcoords";
  options["fields/braid"] = "new_braid";
  options["fields/radial"] = "new_radial";

#if 0
  // Select the left half of zones
  std::vector<int> sel;
  for(int i = 0; i < topoView.numberOfZones(); i++)
  {
    if(i % zoneDims[0] < 5)
      sel.push_back(i);
  }
  options["selectedZones"].set(sel);
#endif

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
    topoView,
    coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
  }

  // Now, take the clipped mesh and clip it again using a mixed topology view.
  using MixedTopoView =
    axom::mir::views::UnstructuredTopologyMixedShapeView<axom::IndexType>;
  using ExpCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
  conduit::Node &n_x =
    deviceClipMesh.fetch_existing("coordsets/clipcoords/values/x");
  conduit::Node &n_y =
    deviceClipMesh.fetch_existing("coordsets/clipcoords/values/y");
  axom::ArrayView<double> xView(static_cast<double *>(n_x.data_ptr()),
                                n_x.dtype().number_of_elements());
  axom::ArrayView<double> yView(static_cast<double *>(n_y.data_ptr()),
                                n_y.dtype().number_of_elements());
  ExpCoordsetView mixedCoordsetView(xView, yView);

  conduit::Node &n_device_topo =
    deviceClipMesh.fetch_existing("topologies/" + clipTopoName);
  conduit::Node &n_conn = n_device_topo.fetch_existing("elements/connectivity");
  conduit::Node &n_shapes = n_device_topo.fetch_existing("elements/shapes");
  conduit::Node &n_sizes = n_device_topo.fetch_existing("elements/sizes");
  conduit::Node &n_offsets = n_device_topo.fetch_existing("elements/offsets");
  axom::ArrayView<axom::IndexType> connView(
    static_cast<axom::IndexType *>(n_conn.data_ptr()),
    n_conn.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> shapesView(
    static_cast<axom::IndexType *>(n_shapes.data_ptr()),
    n_shapes.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> sizesView(
    static_cast<axom::IndexType *>(n_sizes.data_ptr()),
    n_sizes.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> offsetsView(
    static_cast<axom::IndexType *>(n_offsets.data_ptr()),
    n_offsets.dtype().number_of_elements());
  MixedTopoView mixedTopoView(n_device_topo,
                              connView,
                              shapesView,
                              sizesView,
                              offsetsView);

  // Clip the data
  conduit::Node deviceClipMixedMesh;
  axom::mir::clipping::ClipField<ExecSpace, MixedTopoView, ExpCoordsetView>
    mixedClipper(mixedTopoView, mixedCoordsetView);
  options["clipField"] = "new_braid";
  options["clipValue"] = 1.;
  options["fields"].reset();
  options["fields/new_braid"] = "new_braid2";
  options["fields/color"] = "new_color";
  options["fields/new_radial"] = "new_radial2";

  mixedClipper.execute(deviceClipMesh, options, deviceClipMixedMesh);

  // Copy device->host
  conduit::Node hostClipMixedMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMixedMesh,
                                                  deviceClipMixedMesh);

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name + "_mixed"));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMixedMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMixedMesh));
#endif
  }
}

TEST(mir_clipfield, uniform2d)
{
  braid2d_clip_test<seq_exec>("uniform", "uniform2d");

#if defined(AXOM_USE_OPENMP)
  braid2d_clip_test<omp_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid2d_clip_test<cuda_exec>("uniform", "uniform2d");
#endif

#if defined(AXOM_USE_HIP)
  braid2d_clip_test<hip_exec>("uniform", "uniform2d");
#endif
}

template <typename ExecSpace, int NDIMS>
void braid_rectilinear_clip_test(const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, NDIMS>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;

  axom::StackArray<axom::IndexType, NDIMS> dims, zoneDims;
  for(int i = 0; i < NDIMS; i++)
  {
    dims[i] = 10;
    zoneDims[i] = dims[i] - 1;
  }

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid("rectilinear", dims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Create views
  auto coordsetView =
    axom::mir::views::make_rectilinear_coordset<double, NDIMS>::view(
      hostMesh["coordsets/coords"]);
  using CoordsetView = decltype(coordsetView);
  TopoView topoView(Indexing {zoneDims});

  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
    topoView,
    coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
  }
}

TEST(mir_clipfield, rectilinear2d)
{
  braid_rectilinear_clip_test<seq_exec, 2>("rectilinear2d");

#if defined(AXOM_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid_rectilinear_clip_test<cuda_exec, 2>("rectilinear2d");
#endif

#if defined(AXOM_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 2>("rectilinear2d");
#endif
}

TEST(mir_clipfield, rectilinear3d)
{
  braid_rectilinear_clip_test<seq_exec, 3>("rectilinear3d");

#if defined(AXOM_USE_OPENMP)
  braid_rectilinear_clip_test<omp_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid_rectilinear_clip_test<cuda_exec, 3>("rectilinear3d");
#endif

#if defined(AXOM_USE_HIP)
  braid_rectilinear_clip_test<hip_exec, 3>("rectilinear3d");
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace, int NDIMS>
void strided_structured_clip_test(const std::string &name)
{
  using Indexing =
    axom::mir::views::StridedStructuredIndexing<axom::IndexType, NDIMS>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<double, NDIMS>;

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::strided_structured<NDIMS>(hostMesh);
  hostMesh.print();

  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig_yaml", "yaml");
#endif

  conduit::Node options, deviceClipMesh, hostClipMesh;

  // Create options to control the clipping.
  options["clipField"] = "vert_vals";
  options["clipValue"] = 6.5;
  options["inside"] = 1;
  options["outside"] = 1;

  // Create views
  axom::mir::views::dispatch_explicit_coordset(deviceMesh["coordsets/coords"], [&](auto coordsetView)
  {
    auto topoView = axom::mir::views::make_strided_structured<2>::view(deviceMesh["topologies/mesh"]);

    using CoordsetView = decltype(coordsetView);
    using TopoView = decltype(topoView);

    // Clip the data
    axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
      topoView,
      coordsetView);
    clipper.execute(deviceMesh, options, deviceClipMesh);

    // Copy device->host
    axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);
  });

  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
  #if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostClipMesh);
  #else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
  #endif
  }
}

TEST(mir_clipfield, strided_structured_2d)
{
  strided_structured_clip_test<seq_exec, 2>("strided_structured_2d");

  //#if defined(AXOM_USE_OPENMP)
  //  strided_structured_clip_test<omp_exec, 2>("strided_structured_2d");
  //#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  strided_structured_clip_test<cuda_exec, 2>("strided_structured_2d");
#endif

#if defined(AXOM_USE_HIP)
  strided_structured_clip_test<hip_exec, 2>("strided_structured_2d");
#endif
}

template <typename ExecSpace, typename ShapeType>
void braid3d_clip_test(const std::string &type, const std::string &name)
{
  using TopoView =
    axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<double, 3>;

  // Create the data
  const axom::StackArray<axom::IndexType, 3> dims {10, 10, 10};
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Create views
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  const axom::ArrayView<double> x(static_cast<double *>(n_x.data_ptr()),
                                  n_x.dtype().number_of_elements());
  const axom::ArrayView<double> y(static_cast<double *>(n_y.data_ptr()),
                                  n_y.dtype().number_of_elements());
  const axom::ArrayView<double> z(static_cast<double *>(n_z.data_ptr()),
                                  n_z.dtype().number_of_elements());
  CoordsetView coordsetView(x, y, z);

  conduit::Node &n_conn =
    deviceMesh.fetch_existing("topologies/mesh/elements/connectivity");
  const axom::ArrayView<int> conn(static_cast<int *>(n_conn.data_ptr()),
                                  n_conn.dtype().number_of_elements());
  TopoView topoView(conn);

  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
    topoView,
    coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

/// Execute the braid3d test for a single shape on multiple ExecSpaces
template <typename ShapeType>
void braid3d_clip_test_exec(const std::string &type, const std::string &name)
{
  braid3d_clip_test<seq_exec, ShapeType>(type, name);

#if defined(AXOM_USE_OPENMP)
  braid3d_clip_test<omp_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid3d_clip_test<cuda_exec, ShapeType>(type, name);
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec, ShapeType>(type, name);
#endif
}

TEST(mir_clipfield, tet)
{
  braid3d_clip_test_exec<axom::mir::views::TetShape<int>>("tets", "tet");
}

TEST(mir_clipfield, pyramid)
{
  braid3d_clip_test_exec<axom::mir::views::PyramidShape<int>>("pyramids", "pyr");
}

TEST(mir_clipfield, wedge)
{
  braid3d_clip_test_exec<axom::mir::views::WedgeShape<int>>("wedges", "wdg");
}

TEST(mir_clipfield, hex)
{
  braid3d_clip_test_exec<axom::mir::views::HexShape<int>>("hexs", "hex");
}

template <typename ExecSpace>
void braid3d_mixed_clip_test(const std::string &name)
{
  using CoordType = float;
  using ConnType = int;
  using TopoView = axom::mir::views::UnstructuredTopologyMixedShapeView<ConnType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<CoordType, 3>;

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::mixed3d(hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Create views
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  const axom::ArrayView<CoordType> x(static_cast<CoordType *>(n_x.data_ptr()),
                                     n_x.dtype().number_of_elements());
  const axom::ArrayView<CoordType> y(static_cast<CoordType *>(n_y.data_ptr()),
                                     n_y.dtype().number_of_elements());
  const axom::ArrayView<CoordType> z(static_cast<CoordType *>(n_z.data_ptr()),
                                     n_z.dtype().number_of_elements());
  CoordsetView coordsetView(x, y, z);

  conduit::Node &n_device_topo = deviceMesh.fetch_existing("topologies/mesh");
  conduit::Node &n_conn = n_device_topo.fetch_existing("elements/connectivity");
  conduit::Node &n_shapes = n_device_topo.fetch_existing("elements/shapes");
  conduit::Node &n_sizes = n_device_topo.fetch_existing("elements/sizes");
  conduit::Node &n_offsets = n_device_topo.fetch_existing("elements/offsets");
  axom::ArrayView<ConnType> connView(static_cast<ConnType *>(n_conn.data_ptr()),
                                     n_conn.dtype().number_of_elements());
  axom::ArrayView<ConnType> shapesView(
    static_cast<ConnType *>(n_shapes.data_ptr()),
    n_shapes.dtype().number_of_elements());
  axom::ArrayView<ConnType> sizesView(static_cast<ConnType *>(n_sizes.data_ptr()),
                                      n_sizes.dtype().number_of_elements());
  axom::ArrayView<ConnType> offsetsView(
    static_cast<ConnType *>(n_offsets.data_ptr()),
    n_offsets.dtype().number_of_elements());
  TopoView topoView(n_device_topo, connView, shapesView, sizesView, offsetsView);

  // Create options to control the clipping.
  conduit::Node options;
  options["clipField"] = "distance";
  options["clipValue"] = 12.f;
  options["inside"] = 1;
  options["outside"] = 0;

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(
    topoView,
    coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Handle baseline comparison.
  std::string baselineName(yamlRoot(name));
  const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
  saveBaseline(paths, baselineName, hostClipMesh);
#else
  EXPECT_TRUE(compareBaseline(paths, baselineName, hostClipMesh));
#endif
}

TEST(mir_clipfield, mixed)
{
  const std::string name("mixed");
  braid3d_mixed_clip_test<seq_exec>(name);

#if defined(AXOM_USE_OPENMP)
  braid3d_mixed_clip_test<omp_exec>(name);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  braid3d_mixed_clip_test<cuda_exec>(name);
#endif

#if defined(AXOM_USE_HIP)
  braid3d_mixed_clip_test<hip_exec>(name);
#endif
}

//------------------------------------------------------------------------------
#if defined(DEBUGGING_TEST_CASES)
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
#if defined(DEBUGGING_TEST_CASES)
  conduit::utils::set_error_handler(conduit_debug_err_handler);
#endif
  result = RUN_ALL_TESTS();
  return result;
}
