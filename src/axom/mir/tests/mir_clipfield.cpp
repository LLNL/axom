// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>

using seq_exec = axom::SEQ_EXEC;

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

TEST(mir_clipfield, uniform2d)
{
  using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
  using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
  using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

  axom::StackArray<axom::IndexType, 2> dims{10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims{dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node mesh;
  braid("uniform", dims, mesh);
  conduit::relay::io::blueprint::save_mesh(mesh, "uniform2d_orig", "hdf5");

  // Create views
  axom::StackArray<double, 2> origin{0., 0.}, spacing{1., 1.};
  CoordsetView coordsetView(dims, origin, spacing);
  TopoView topoView(Indexing{zoneDims});

  // Clip the data
  conduit::Node clipmesh;
  axom::mir::clipping::ClipField<seq_exec, TopoView, CoordsetView> clipper(topoView, coordsetView);
  conduit::Node options;
  options["clipField"] = "distance";
  options["inside"] = 1;
  options["outside"] = 1;
  options["topologyName"] = "cliptopo";
  options["coordsetName"] = "clipcoords";
  options["fields/braid"] = "new_braid";
  options["fields/radial"] = "new_radial";

  clipper.execute(mesh, options, clipmesh);
  conduit::relay::io::blueprint::save_mesh(clipmesh, "uniform2d", "hdf5");
  conduit::relay::io::blueprint::save_mesh(clipmesh, "uniform2d_yaml", "yaml");

  // Load a clipped baseline file & compare.
}


//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
