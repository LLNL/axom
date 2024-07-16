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


TEST(mir_clipfield, uniform2d)
{
   using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
   using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
   using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

   axom::StackArray<axom::IndexType, 2> dims{10, 10};

   // Create the data
   conduit::Node mesh;
   braid("uniform", dims, mesh);
   conduit::relay::io::blueprint::save_mesh(mesh, "uniform2d_orig", "hdf5");
   mesh.print();
   const conduit::Node &topo = mesh["topologies"][0];
   const conduit::Node &coordset = mesh["coordsets"][0];

   // Create views
   axom::StackArray<double, 2> origin{0., 0.}, spacing{1., 1.};
   CoordsetView coordsetView(dims, origin, spacing);
   Indexing zoneIndexing(dims);
   TopoView topoView(zoneIndexing);

   // Clip the data
   conduit::Node clipmesh;
   axom::mir::clipping::ClipField<seq_exec, TopoView, CoordsetView> clipper(topoView, coordsetView);
   clipper.execute(mesh, "distance", clipmesh);
   conduit::relay::io::blueprint::save_mesh(clipmesh, "uniform2d", "hdf5");

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
