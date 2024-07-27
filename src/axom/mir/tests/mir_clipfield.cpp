// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>

// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA)
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

template <typename Dimensions>
void braid(const std::string &type, const Dimensions &dims, conduit::Node &mesh)
{
  int d[3] = {0, 0, 0};
  for(int i = 0; i < dims.size(); i++)
    d[i] = dims[i];
  conduit::blueprint::mesh::examples::braid(type, d[0], d[1], d[2], mesh);

  // Make a new distance field.
  const float dist = 6.5f;
  const conduit::Node &n_coordset = mesh["coordsets"][0];
  axom::mir::views::dispatch_coordset(n_coordset, [&](auto coordsetView)
  {
    mesh["fields/distance/topology"] = "mesh";
    mesh["fields/distance/association"] = "vertex";
    conduit::Node &n_values = mesh["fields/distance/values"];
    const auto nnodes = coordsetView.size();
    n_values.set(conduit::DataType::float32(nnodes));
    float *valuesPtr = static_cast<float *>(n_values.data_ptr());
    for(int index = 0; index < nnodes; index++)
    {
      const auto pt = coordsetView[index];
      float norm2 = 0.f;
      for(int i = 0; i < pt.DIMENSION; i++)
        norm2 += pt[i] * pt[i];
      valuesPtr[index] = sqrt(norm2) - dist;
    }
  });
}

int conduit_to_vtk_cell(int shape_value)
{
  int vtktype = 0;
  if(shape_value == axom::mir::views::Tri_ShapeID)
      vtktype = 5; // VTK_TRIANGLE
  else if(shape_value == axom::mir::views::Quad_ShapeID)
      vtktype = 9; // VTK_QUAD
  else if(shape_value == axom::mir::views::Polygon_ShapeID)
      vtktype = 7; // VTK_POLYGON
  else if(shape_value == axom::mir::views::Tet_ShapeID)
      vtktype = 10; // VTK_TETRA
  else if(shape_value == axom::mir::views::Pyramid_ShapeID)
      vtktype = 14; // VTK_PYRAMID
  else if(shape_value == axom::mir::views::Wedge_ShapeID)
      vtktype = 13; // VTK_WEDGE
  else if(shape_value == axom::mir::views::Hex_ShapeID)
      vtktype = 12; // VTK_HEXAHEDRON

  return vtktype;
}

void
conduit_save_vtk(const conduit::Node &node, const std::string &path)
{
    FILE *file = fopen(path.c_str(), "wt");

    // Write the VTK file header
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Unstructured Grid Example\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    // Write the points
    const conduit::Node &points = node["coordsets/coords/values"];
    const conduit::Node &x = points["x"];
    const conduit::Node &y = points["y"];
    const conduit::Node &z = points["z"];
    size_t num_points = x.dtype().number_of_elements();

    fprintf(file, "POINTS %zu float\n", num_points);
    for (size_t i = 0; i < num_points; ++i) {
        fprintf(file, "%f %f %f\n", x.as_float64_array()[i], y.as_float64_array()[i], z.as_float64_array()[i]);
    }

    // Write the cells
    const conduit::Node &topologies = node["topologies"];
    const conduit::Node &topo = topologies[0];
    const conduit::Node &elements = topo["elements"];
    const conduit::Node &connectivity = elements["connectivity"];
    size_t num_cells = elements["sizes"].dtype().number_of_elements();
    size_t total_num_indices = connectivity.dtype().number_of_elements();

    fprintf(file, "CELLS %zu %zu\n", num_cells, total_num_indices + num_cells);
    size_t index = 0;
    for (size_t i = 0; i < num_cells; ++i) {
        size_t cell_size = elements["sizes"].as_int32_array()[i];
        fprintf(file, "%zu", cell_size);
        for (size_t j = 0; j < cell_size; ++j) {
            fprintf(file, " %d", connectivity.as_int32_array()[index++]);
        }
        fprintf(file, "\n");
    }

    // Write the cell types
    const conduit::Node &shapes = elements["shapes"];
    fprintf(file, "CELL_TYPES %zu\n", num_cells);
    for (size_t i = 0; i < num_cells; ++i) {
        const auto type = conduit_to_vtk_cell(shapes.as_int32_array()[i]);
        fprintf(file, "%d\n", type);
    }

    // Close the file
    fclose(file);
}



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
  n_fields["distance/values"].set(std::vector<float>{0.,1.,2.,3.});

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
  options["selectedZones"].set(std::vector<axom::IndexType>{5,4,3});
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
  axom::Array<IndexType> blendGroups{{8, 5}};
  axom::Array<IndexType> blendGroupsLen{{/*zone 0*/ 4+1+1+1+1+2+2+2, /*zone 1*/ 1+1+1+1+2}};
  axom::Array<IndexType> blendGroupOffsets{{0, 8}};
  axom::Array<IndexType> blendOffsets{{0, blendGroupsLen[0]}};

  axom::Array<KeyType>   blendNames{{/*zone 0*/ 6,0,1,3,4,7,8,9, /*zone 1*/ 1,2,4,5,9}};
  axom::Array<IndexType> blendGroupSizes{{/*zone 0*/4,1,1,1,1,2,2,2, /*zone 1*/1,1,1,1,2}};
  axom::Array<IndexType> blendGroupStart{{/*zone 0*/0,4,5,6,7,8,10,12, /*zone 1*/ 13, 14, 15, 16, 18}};
  axom::Array<IndexType> blendIds{{/*zone 0*/
                                   0,1,2,3,   // 6 (bgname) // 0 (bgindex)
                                   0,         // 0          // 1
                                   1,         // 1          // 2
                                   3,         // 3          // 3
                                   4,         // 4          // 4
                                   0,1,       // 7          // 5
                                   3,4,       // 8          // 6
                                   1,4,       // 9          // 7
                                   /*zone 1*/
                                   1,         // 1          // 8
                                   2,         // 2          // 9
                                   4,         // 4          // 10
                                   5,         // 5          // 11
                                   1,4        // 9          // 12
                                 }};
  axom::Array<float> blendCoeff{{/*zone 0*/
                                   0.25, 0.25, 0.25, 0.25,
                                   1.,
                                   1.,
                                   1.,
                                   1.,
                                   0.5, 0.5,
                                   0.5, 0.5,
                                   0.5, 0.5,
                                   /*zone 1*/
                                   1.,
                                   1.,
                                   1.,
                                   1.,
                                   0.5, 0.5
                                 }};
  axom::Array<KeyType>   blendUniqueNames{{0,1,2,3,4,5,6,7,8,9}};
  axom::Array<KeyType>   blendUniqueIndices{{1,2,9,3,4,11,0,5,6,7}};

  axom::mir::clipping::BlendGroupBuilder<seq_exec> builder;
  builder.setBlendGroupSizes(blendGroups.view(), blendGroupsLen.view());
  builder.setBlendGroupOffsets(blendOffsets.view(), blendGroupOffsets.view());
  builder.setBlendViews(blendNames.view(), blendGroupSizes.view(), blendGroupStart.view(), blendIds.view(), blendCoeff.view());

  std::cout << "-------- zone 0 --------" << std::endl;
  auto z0 = builder.blendGroupsForZone(0);
  EXPECT_EQ(z0.size(), 8);
  IndexType index = 0;
  for(IndexType i = 0; i < z0.size(); i++, index++)
  {
    z0.print(std::cout);
    EXPECT_EQ(z0.ids().size(), blendGroupSizes[index]);

    z0++;
  }

  std::cout << "-------- zone 1 --------" << std::endl;
  auto z1 = builder.blendGroupsForZone(1);
  EXPECT_EQ(z1.size(), 5);
  for(IndexType i = 0; i < z1.size(); i++, index++)
  {
    z1.print(std::cout);
    EXPECT_EQ(z1.ids().size(), blendGroupSizes[index]);
    z1++;
  }
}

template <typename ExecSpace>
void braid2d_clip_test(const std::string &type, const std::string &name)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

  axom::StackArray<axom::IndexType, 2> dims{10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims{dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  braid(type, dims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");

  // Create views
  axom::StackArray<double, 2> origin{0., 0.}, spacing{1., 1.};
  CoordsetView coordsetView(dims, origin, spacing);
  TopoView topoView(Indexing{zoneDims});

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
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Save data.
  conduit::relay::io::blueprint::save_mesh(hostClipMesh, name, "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostClipMesh, name + "_yaml", "yaml");

  // Now, take the clipped mesh and clip it again using a mixed topology view.
  using MixedTopoView = axom::mir::views::UnstructuredTopologyMixedShapeView<axom::IndexType>;
  using ExpCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
  conduit::Node &n_x = deviceClipMesh.fetch_existing("coordsets/clipcoords/values/x");
  conduit::Node &n_y = deviceClipMesh.fetch_existing("coordsets/clipcoords/values/y");
  axom::ArrayView<double> xView(static_cast<double *>(n_x.data_ptr()), n_x.dtype().number_of_elements());
  axom::ArrayView<double> yView(static_cast<double *>(n_y.data_ptr()), n_y.dtype().number_of_elements());
  ExpCoordsetView mixedCoordsetView(xView, yView);

  conduit::Node &n_device_topo = deviceClipMesh.fetch_existing("topologies/" + clipTopoName);
  conduit::Node &n_conn = n_device_topo.fetch_existing("elements/connectivity");
  conduit::Node &n_shapes = n_device_topo.fetch_existing("elements/shapes");
  conduit::Node &n_sizes = n_device_topo.fetch_existing("elements/sizes");
  conduit::Node &n_offsets = n_device_topo.fetch_existing("elements/offsets");
  axom::ArrayView<axom::IndexType> connView(static_cast<axom::IndexType *>(n_conn.data_ptr()), n_conn.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> shapesView(static_cast<axom::IndexType *>(n_shapes.data_ptr()), n_shapes.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> sizesView(static_cast<axom::IndexType *>(n_sizes.data_ptr()), n_sizes.dtype().number_of_elements());
  axom::ArrayView<axom::IndexType> offsetsView(static_cast<axom::IndexType *>(n_offsets.data_ptr()), n_offsets.dtype().number_of_elements());
  MixedTopoView mixedTopoView(n_device_topo, connView, shapesView, sizesView, offsetsView);

  // Clip the data
  conduit::Node deviceClipMixedMesh; 
  axom::mir::clipping::ClipField<ExecSpace, MixedTopoView, ExpCoordsetView> mixedClipper(mixedTopoView, mixedCoordsetView);
  options["clipField"] = "new_braid";
  options["clipValue"] = 1.;
  options["fields"].reset();
  options["fields/new_braid"] = "new_braid2";
  options["fields/color"] = "new_color";
  options["fields/new_radial"] = "new_radial2";

  mixedClipper.execute(deviceClipMesh, options, deviceClipMixedMesh);

  // Copy device->host
  conduit::Node hostClipMixedMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMixedMesh, deviceClipMixedMesh);

  // Save data.
  conduit::relay::io::blueprint::save_mesh(hostClipMixedMesh, name + "_mixed", "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostClipMixedMesh, name + "_mixed_yaml", "yaml");

  // Load a clipped baseline file & compare.
}

template <typename ExecSpace, typename ShapeType>
void braid3d_clip_test(const std::string &type, const std::string &name)
{
  using TopoView = axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  using CoordsetView = axom::mir::views::ExplicitCoordsetView<double, 3>;

  axom::StackArray<axom::IndexType, 3> dims{8,8,8};//10, 10, 10};
  axom::StackArray<axom::IndexType, 3> zoneDims{dims[0] - 1, dims[1] - 1, dims[2] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  braid(type, dims, hostMesh);
  hostMesh.print();
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig_yaml", "yaml");

  // Create views
  conduit::Node &n_x = deviceMesh.fetch_existing("coordsets/coords/values/x");
  conduit::Node &n_y = deviceMesh.fetch_existing("coordsets/coords/values/y");
  conduit::Node &n_z = deviceMesh.fetch_existing("coordsets/coords/values/z");
  const axom::ArrayView<double> x(static_cast<double *>(n_x.data_ptr()), n_x.dtype().number_of_elements());
  const axom::ArrayView<double> y(static_cast<double *>(n_y.data_ptr()), n_y.dtype().number_of_elements());
  const axom::ArrayView<double> z(static_cast<double *>(n_z.data_ptr()), n_z.dtype().number_of_elements());
  CoordsetView coordsetView(x, y, z);

  conduit::Node &n_conn = deviceMesh.fetch_existing("topologies/mesh/elements/connectivity");
  //conduit::Node &n_sizes = deviceMesh.fetch_existing("topologies/mesh/elements/sizes");
  //conduit::Node &n_offsets = deviceMesh.fetch_existing("topologies/mesh/elements/offsets");
  const axom::ArrayView<int> conn(static_cast<int *>(n_conn.data_ptr()), n_conn.dtype().number_of_elements());
  //const axom::ArrayView<int> sizes(static_cast<int *>(n_sizes.data_ptr()), n_sizes.dtype().number_of_elements());
  //const axom::ArrayView<int> offsets(static_cast<int *>(n_offsets.data_ptr()), n_offsets.dtype().number_of_elements()); 
  TopoView topoView(conn); //, sizes, offsets);

  // Create options to control the clipping.
  const std::string clipTopoName("cliptopo");
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 0;
  //options["topologyName"] = clipTopoName;
  //options["coordsetName"] = "clipcoords";
  //options["fields/braid"] = "new_braid";
  //options["fields/radial"] = "new_radial";

  // Clip the data
  conduit::Node deviceClipMesh;
  axom::mir::clipping::ClipField<ExecSpace, TopoView, CoordsetView> clipper(topoView, coordsetView);
  clipper.execute(deviceMesh, options, deviceClipMesh);

  // Copy device->host
  conduit::Node hostClipMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostClipMesh, deviceClipMesh);

  // Save data.
  conduit::relay::io::blueprint::save_mesh(hostClipMesh, name, "hdf5");
  conduit::relay::io::blueprint::save_mesh(hostClipMesh, name + "_yaml", "yaml");
  conduit_save_vtk(hostClipMesh, name + ".vtk");
}

#if 0
TEST(mir_clipfield, uniform2d)
{
  braid2d_clip_test<seq_exec>("uniform", "uniform2d");

//#if defined(AXOM_USE_OPENMP)
//  braid2d_clip_test<omp_exec>("uniform", "uniform2d_omp");
//#endif

#if defined(AXOM_USE_CUDA)
  braid2d_clip_test<cuda_exec>("uniform", "uniform2d_cuda");
#endif

#if defined(AXOM_USE_HIP)
  braid2d_clip_test<hip_exec>("uniform", "uniform2d_hip");
#endif
}
#endif

#if 0
TEST(mir_clipfield, hex)
{
  using ShapeType = axom::mir::views::HexShape<int>;

  const std::string type("hexs");
  braid3d_clip_test<seq_exec, ShapeType>(type, "hex");

//#if defined(AXOM_USE_OPENMP)
//  braid3d_clip_test<omp_exec>(type, "hex_omp");
//#endif

#if defined(AXOM_USE_CUDA)
  braid3d_clip_test<cuda_exec>(type, "hex_cuda");
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec>(type, "hex_hip");
#endif
}
#endif

#if 0
TEST(mir_clipfield, tet)
{
  using ShapeType = axom::mir::views::TetShape<int>;

  const std::string type("tets");
  braid3d_clip_test<seq_exec, ShapeType>(type, "tet");

//#if defined(AXOM_USE_OPENMP)
//  braid3d_clip_test<omp_exec>(type, "tet_omp");
//#endif

#if defined(AXOM_USE_CUDA)
  braid3d_clip_test<cuda_exec>(type, "tet_cuda");
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec>(type, "tet_hip");
#endif
}
#endif

#if 0
TEST(mir_clipfield, pyramid)
{
  using ShapeType = axom::mir::views::PyramidShape<int>;

  const std::string type("pyramids");
  braid3d_clip_test<seq_exec, ShapeType>(type, "pyr");

//#if defined(AXOM_USE_OPENMP)
//  braid3d_clip_test<omp_exec>(type, "pyr_omp");
//#endif

#if defined(AXOM_USE_CUDA)
  braid3d_clip_test<cuda_exec>(type, "pyr_cuda");
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec>(type, "pyr_hip");
#endif
}
#endif

#if 1
TEST(mir_clipfield, wedge)
{
  using ShapeType = axom::mir::views::WedgeShape<int>;

  const std::string type("wedges");
  braid3d_clip_test<seq_exec, ShapeType>(type, "wdg");

//#if defined(AXOM_USE_OPENMP)
//  braid3d_clip_test<omp_exec>(type, "wdg_omp");
//#endif

#if defined(AXOM_USE_CUDA)
  braid3d_clip_test<cuda_exec>(type, "wdg_cuda");
#endif

#if defined(AXOM_USE_HIP)
  braid3d_clip_test<hip_exec>(type, "wdg_hip");
#endif
}
#endif
//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
