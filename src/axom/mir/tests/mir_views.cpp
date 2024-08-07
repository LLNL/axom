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
  const axom::Array<int> expectedZones {{16, 17, 24, 23, 17, 18, 25, 24,
                                         18, 19, 26, 25, 23, 24, 31, 30,
                                         24, 25, 32, 31, 25, 26, 33, 32}};
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

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
