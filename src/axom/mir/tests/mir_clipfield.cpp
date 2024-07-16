// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>


using seq_exec = axom::SEQ_EXEC;

TEST(mir_clipfield, uniform2d)
{
   using Indexing = axom::mir::views::StructuredIndexing<axom::IndexType, 2>;
   using TopoView = axom::mir::views::StructuredTopologyView<Indexing>;
   using CoordsetView = axom::mir::views::UniformCoordsetView<double, 2>;

   axom::StackArray<axom::IndexType, 2> dims{10, 10};

   // Create the data
   conduit::Node mesh;
   conduit::blueprint::mesh::examples::braid("uniform", dims[0], dims[1], 0, mesh);
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
   clipper.execute(mesh, "radial", clipmesh);
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
