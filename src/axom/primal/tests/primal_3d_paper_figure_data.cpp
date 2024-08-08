// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! 
 * \file primal_bezier_curve.cpp
 * \brief This file tests primal's Bezier curve functionality
 */

#include "gtest/gtest.h"

#include "axom/core/FlatMap.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"
#include "axom/primal/operators/winding_number.hpp"
#include "axom/primal/operators/detail/winding_number_impl.hpp"

#include "axom/primal/operators/printers.hpp"

#include <chrono>
#include <string>
#include <iomanip>

namespace primal = axom::primal;
using Point3D = primal::Point<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using BPatch = primal::BezierPatch<double, 3>;
using BBox = primal::BoundingBox<double, 3>;

TEST(primal_3d_paper_figure_data, plotting_demo)
{
  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\plotting_demo";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  double rt2 = sqrt(2), rt3 = sqrt(3), rt6 = sqrt(6);

  // clang-format off
  axom::Array<Point3D> node_data = {
    Point3D {4*(1-rt3),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(rt3-4),            -rt2, rt2*(rt3-4)}, Point3D {4*(1-2*rt3)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(rt3-4),           rt2,   rt2*(rt3-4)}, Point3D {4*(1-rt3),     4*(rt3-1),     4*(1-rt3)},
    Point3D {     -rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(2-3*rt3)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(2*rt3-7)/3, 0,      -5*rt6/3}, Point3D {(2-3*rt3)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {     -rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
    Point3D {        0, 4*(1-2*rt3)/3, 4*(1-2*rt3)/3}, Point3D {          0, rt2*(2*rt3-7)/3,    -5*rt6/3}, Point3D {0,               0,   4*(rt3-5)/3}, Point3D {          0, rt2*(7-2*rt3)/3,    -5*rt6/3}, Point3D {        0, 4*(2*rt3-1)/3, 4*(1-2*rt3)/3},
    Point3D {      rt2, rt2*(rt3 - 4), rt2*(rt3 - 4)}, Point3D {(3*rt3-2)/2,     (2-3*rt3)/2,  -(rt3+6)/2}, Point3D {rt2*(7-2*rt3)/3, 0,      -5*rt6/3}, Point3D {(3*rt3-2)/2,   (3*rt3-2)/2,    -(rt3+6)/2}, Point3D {      rt2,   rt2*(4-rt3),   rt2*(rt3-4)},
    Point3D {4*(rt3-1),     4*(1-rt3),     4*(1-rt3)}, Point3D {rt2*(4-rt3),            -rt2, rt2*(rt3-4)}, Point3D {4*(2*rt3-1)/3,   0, 4*(1-2*rt3)/3}, Point3D {rt2*(4-rt3),           rt2,   rt2*(rt3-4)}, Point3D {4*(rt3-1),     4*(rt3-1),     4*(1-rt3)}};

  axom::Array<double> weight_data = {
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3),
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
       4*(5-rt3)/3, rt2*(rt3+6)/3, 4*(5*rt3-1)/9, rt2*(rt3+6)/3,   4*(5-rt3)/3,
     rt2*(3*rt3-2),     (rt3+6)/2, rt2*(rt3+6)/3,     (rt3+6)/2, rt2*(3*rt3-2),
         4*(3-rt3), rt2*(3*rt3-2),   4*(5-rt3)/3, rt2*(3*rt3-2),     4*(3-rt3)};
  // clang-format on

  BPatch sphere_faces[6];
  for(int n = 0; n < 6; ++n)
  {
    sphere_faces[n].setOrder(4, 4);
    sphere_faces[n].makeRational();
  }

  sphere_faces[0].setOrder(4, 4);
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      for(int n = 0; n < 6; ++n)
      {
        sphere_faces[n].setWeight(i, j, weight_data[idx]);
      }

      // Set up each face by rotating one of the patch faces
      sphere_faces[0](i, j)[0] = node_data[idx][1];
      sphere_faces[0](i, j)[1] = node_data[idx][0];
      sphere_faces[0](i, j)[2] = node_data[idx][2];
      sphere_faces[0](i, j).array() /= weight_data[idx];

      sphere_faces[1](i, j)[0] = -node_data[idx][0];
      sphere_faces[1](i, j)[1] = -node_data[idx][1];
      sphere_faces[1](i, j)[2] = -node_data[idx][2];
      sphere_faces[1](i, j).array() /= weight_data[idx];

      sphere_faces[2](i, j)[0] = node_data[idx][2];
      sphere_faces[2](i, j)[1] = node_data[idx][1];
      sphere_faces[2](i, j)[2] = node_data[idx][0];
      sphere_faces[2](i, j).array() /= weight_data[idx];

      sphere_faces[3](i, j)[0] = -node_data[idx][1];
      sphere_faces[3](i, j)[1] = -node_data[idx][2];
      sphere_faces[3](i, j)[2] = -node_data[idx][0];
      sphere_faces[3](i, j).array() /= weight_data[idx];

      sphere_faces[4](i, j)[0] = node_data[idx][0];
      sphere_faces[4](i, j)[1] = node_data[idx][2];
      sphere_faces[4](i, j)[2] = node_data[idx][1];
      sphere_faces[4](i, j).array() /= weight_data[idx];

      sphere_faces[5](i, j)[0] = -node_data[idx][2];
      sphere_faces[5](i, j)[1] = -node_data[idx][0];
      sphere_faces[5](i, j)[2] = -node_data[idx][1];
      sphere_faces[5](i, j).array() /= weight_data[idx];
    }
  }

  BBox bbox;
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      for(int n = 0; n < 6; ++n)
      {
        bbox.addPoint(sphere_faces[n](i, j));
      }
    }
  }

  axom::Array<BPatch> patches;
  patches.push_back(sphere_faces[0]);

  auto wn_accurate = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    double wn = 0.0;
    for(const auto& patch : patches)
      wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
    return wn;
  };

  auto wn_stokes = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    double wn = 0.0;
    for(const auto& patch : patches)
      wn += just_the_stokes(query, patch, edge_tol, quad_tol, EPS);
    return wn;
  };

  // clang-format off
  exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_accurate, bbox, 50, 50, 50);
  exportScalarFieldToVTK(data_dir + "/sphere_field_stokes.vtk", wn_stokes, bbox, 50, 50, 50);
  // clang-format on

  exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
