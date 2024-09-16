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
using Point2D = primal::Point<double, 2>;
using Point3D = primal::Point<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using BCurve2D = primal::BezierCurve<double, 2>;
using BPatch = primal::BezierPatch<double, 3>;
using BBox = primal::BoundingBox<double, 3>;
using Tri = primal::Triangle<double, 3>;

TEST(primal_3d_paper_figure_data, plotting_demo)
{
  return;

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

  double the_wn = winding_number_casting(Point3D {0.1, 0.1},
                                         sphere_faces[1],
                                         edge_tol,
                                         quad_tol,
                                         EPS);
  std::cout << the_wn << std::endl;

  return;

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
      wn += winding_number_casting(query, patch, edge_tol, quad_tol, EPS);
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
  // exportScalarFieldToVTK(data_dir + "/sphere_field_stokes.vtk", wn_stokes, bbox, 50, 50, 50);
  // clang-format on

  exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches);
}

TEST(primal_3d_paper_figure_data, rotating_patch)
{
  //   return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\rotating_patch";

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

  BPatch sphere_face;
  sphere_face.setOrder(4, 4);
  sphere_face.makeRational();

  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      const int idx = 5 * i + j;
      sphere_face.setWeight(i, j, weight_data[idx]);

      // Set up each face by rotating one of the patch faces
      sphere_face(i, j)[0] = -node_data[idx][0];
      sphere_face(i, j)[1] = -node_data[idx][1];
      sphere_face(i, j)[2] = -node_data[idx][2];
      sphere_face(i, j).array() /= weight_data[idx];
    }
  }

  auto the_wn =
    winding_number_casting_split(Point3D {0.63494712, -0.427279, 1.5401849},
                                 sphere_face,
                                 edge_tol,
                                 quad_tol,
                                 EPS);

  std::cout << "----------------" << std::endl;
  std::cout << the_wn.first << std::endl;
  std::cout << the_wn.second << std::endl;
  std::cout << the_wn.first + the_wn.second << std::endl;
  std::cout << "----------------" << std::endl;

  auto true_wn = winding_number(Point3D {0.63494712, -0.427279, 1.5401849},
                                sphere_face,
                                edge_tol,
                                quad_tol,
                                EPS);

  std::cout << true_wn << std::endl;
  std::cout << true_wn << " - " << the_wn.first + the_wn.second << " = " << true_wn - (the_wn.first + the_wn.second) << std::endl;
  return;

  BBox bbox;
  for(int i = 0; i < 5; ++i)
  {
    for(int j = 0; j < 5; ++j)
    {
      bbox.addPoint(sphere_face(i, j));
    }
  }

  // Make the bounding box a cube with the same centroid
  int max_dim = bbox.getLongestDimension();
  double max_len = bbox.getMax()[max_dim] - bbox.getMin()[max_dim];

  primal::Point<double, 3> centroid = bbox.getCentroid();

  primal::Point<double, 3> new_min {centroid[0] - max_len / 2.0,
                                    centroid[1] - max_len / 2.0,
                                    centroid[2] - max_len / 2.0};
  primal::Point<double, 3> new_max {centroid[0] + max_len / 2.0,
                                    centroid[1] + max_len / 2.0,
                                    centroid[2] + max_len / 2.0};

  bbox.addPoint(new_min);
  bbox.addPoint(new_max);
  bbox.scale(1.1);

  std::cout << bbox << std::endl;
  axom::Array<BPatch> patches;
  patches.push_back(sphere_face);

  auto wn_ground_truth =
    [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : patches)
        wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto wn_casting = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    std::pair<double, double> wn_split = {0.0, 0.0};
    for(const auto& patch : patches)
    {
      auto val =
        winding_number_casting_split(query, patch, edge_tol, quad_tol, EPS);
      wn_split.first += val.first;
      wn_split.second += val.second;
    }
    return wn_split;
  };

  // clang-format off
  exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_ground_truth, bbox, 100, 100, 100);
  exportSplitScalarFieldToVTK(data_dir + "/sphere_field_casting_z.vtk", wn_casting, bbox, 100, 100, 100);
  // clang-format on

  exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches);
}

TEST(primal_3d_paper_figure_data, trimmed_sphere)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\trimmed_sphere";

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

  // Create a patch that has the inner 90% of sphere_faces[5]
  BPatch p1, p2, p3, p4, p5;
  sphere_faces[5].split(0.9, 0.9, p1, p2, p3, p4);
  p1.split(0.09, 0.09, p2, p3, p4, p5);
  sphere_faces[5] = p5;

  axom::Array<BPatch> patches;
  for(int n = 0; n < 6; ++n) patches.push_back(sphere_faces[n]);

  auto wn_accurate = [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
    double wn = 0.0;
    for(const auto& patch : patches)
      wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
    return wn;
  };

  exportSurfaceToSTL(data_dir + "/sphere_face.stl", patches, 50, 50);

  // clang-format off
  exportScalarFieldToVTK(data_dir + "/sphere_field.vtk", wn_accurate, bbox, 300, 300, 300);
  // clang-format on
}

TEST(primal_3d_paper_figure_data, vase_shape)
{
  return;
  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\vase_shape";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // lambda function that rotates a given BezierCurve around the z-axis
  auto rotate_curve = [](const BCurve2D& curve) -> axom::Array<BPatch> {
    const int ord = curve.getOrder();
    axom::Array<BPatch> rs(4);
    for(int i = 0; i < 4; ++i)
    {
      rs[i].setOrder(ord, 2);
      rs[i].makeRational();
    }

    for(int i = 0; i <= ord; ++i)
    {
      rs[0](i, 0) = Point3D {curve[i][0], 0.0, curve[i][1]};
      rs[0](i, 1) = Point3D {curve[i][0], curve[i][0], curve[i][1]};
      rs[0](i, 2) = Point3D {0.0, curve[i][0], curve[i][1]};

      rs[1](i, 0) = Point3D {0.0, curve[i][0], curve[i][1]};
      rs[1](i, 1) = Point3D {-curve[i][0], curve[i][0], curve[i][1]};
      rs[1](i, 2) = Point3D {-curve[i][0], 0.0, curve[i][1]};

      rs[2](i, 0) = Point3D {-curve[i][0], 0.0, curve[i][1]};
      rs[2](i, 1) = Point3D {-curve[i][0], -curve[i][0], curve[i][1]};
      rs[2](i, 2) = Point3D {0.0, -curve[i][0], curve[i][1]};

      rs[3](i, 0) = Point3D {0.0, -curve[i][0], curve[i][1]};
      rs[3](i, 1) = Point3D {curve[i][0], -curve[i][0], curve[i][1]};
      rs[3](i, 2) = Point3D {curve[i][0], 0.0, curve[i][1]};

      for(int j = 0; j < 4; ++j)
      {
        rs[j].setWeight(i, 1, 1.0 / std::sqrt(2));
      }
    }

    return rs;
  };

  BCurve2D vase1(5);
  vase1[0] = Point2D {0.0, -1.0};
  vase1[1] = Point2D {1.0, -1.0};
  vase1[2] = Point2D {0.5, 0.5};
  vase1[3] = Point2D {0.0, 0.8};
  vase1[4] = Point2D {0.8, 0.5};
  vase1[5] = Point2D {0.8, 1.0};

  BCurve2D vase2(5);
  vase2[0] = Point2D {0.8, 1.0};
  vase2[1] = Point2D {0.8, 1.5};
  vase2[2] = Point2D {0.0, 0.8};
  vase2[3] = Point2D {0.3, 1.5};
  vase2[4] = Point2D {0.5, -0.8};
  vase2[5] = Point2D {0.0, -0.8};

  auto vase1_patches = rotate_curve(vase1);
  auto vase2_patches = rotate_curve(vase2);

  BCurve2D trim(2);
  trim[0] = Point2D {0.55, 0.2};
  trim[1] = Point2D {1.0, 0.4};
  trim[2] = Point2D {0.4, 0.4};

  auto trim_patches = rotate_curve(trim);

  BBox bbox;
  bbox.addPoint(Point3D {1.25, 1.25, -1.25});
  bbox.addPoint(Point3D {-1.25, -1.25, 1.75});

  axom::Array<BPatch> patches;
  for(int n = 0; n < 4; ++n)
  {
    patches.push_back(vase1_patches[n]);
    patches.push_back(vase2_patches[n]);
  }

  auto vase_wn_accurate =
    [&patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : patches)
        wn += winding_number(query, patch, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto trim_wn_accurate =
    [&trim_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& trim_patches : trim_patches)
        wn += winding_number(query, trim_patches, edge_tol, quad_tol, EPS);
      return wn;
    };

  auto bad_query = Point3D {-0.38793103448275867429,
                            -0.56034482758620685061,
                            -0.93965517241379314939};

  for(int i = 0; i < patches.size(); ++i)
  {
    axom::Array<BPatch> single_patch;
    single_patch.push_back(patches[i]);
    exportSurfaceToSTL(data_dir + "/vase_shape_" + std::to_string(i) + ".stl",
                       single_patch,
                       50,
                       50);
  }

  for(int i = 0; i < trim_patches.size(); ++i)
  {
    axom::Array<BPatch> single_patch;
    single_patch.push_back(trim_patches[i]);
    exportSurfaceToSTL(data_dir + "/trim_shape_" + std::to_string(i) + ".stl",
                       single_patch,
                       50,
                       50);
  }

  exportSurfaceToSTL(data_dir + "/vase_shape.stl", patches, 50, 50);
  exportSurfaceToSTL(data_dir + "/trim_shape.stl", trim_patches, 50, 50);

  // clang-format off
  //exportScalarFieldToVTK(data_dir + "/vase_field.vtk", vase_wn_accurate, bbox, 300, 300, 300);
  //exportScalarFieldToVTK(data_dir + "/trim_field.vtk", trim_wn_accurate, bbox, 300, 300, 300);
  // clang-format on
}

TEST(primal_3d_paper_figure_data, two_teapots)
{
  return;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_3d\\teapot_shape\\";

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  std::ifstream file(data_dir + "teapot.stl");

  BBox bbox;
  bbox.addPoint(Point3D {0, 0, -55});
  bbox.addPoint(Point3D {160, 90, 55});

  if(!file.is_open())
  {
    std::cerr << "Unable to open file" << std::endl;
    return;
  }

  std::string line;
  Point3D p1, p2;
  Point3D p3;
  int vertexIndex = 0;

  axom::Array<Tri> half_patches;
  axom::Array<Tri> mirror_patches;

  while(std::getline(file, line))
  {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;

    if(keyword == "vertex")
    {
      Point3D the_point;
      iss >> the_point[0] >> the_point[1] >> the_point[2];

      switch(vertexIndex)
      {
      case 0:
        p1 = the_point;
        break;
      case 1:
        p2 = the_point;
        break;
      case 2:
        p3 = the_point;
        break;
      }

      vertexIndex++;
    }
    else if(keyword == "endloop")
    {
      half_patches.push_back(Tri {p1, p2, p3});
      p1[2] = -p1[2];
      p2[2] = -p2[2];
      p3[2] = -p3[2];
      mirror_patches.push_back(Tri {p1, p2, p3});
      vertexIndex = 0;
    }
  }

  auto half_teapot_wn =
    [&half_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : half_patches) wn += winding_number(query, patch);
      return wn;
    };

  auto mirror_teapot_wn =
    [&mirror_patches, &edge_tol, &quad_tol, &EPS](const Point3D& query) {
      double wn = 0.0;
      for(const auto& patch : mirror_patches)
        wn += winding_number(query, patch);
      return wn;
    };

  exportScalarFieldToVTK(data_dir + "/teapot_half_field_high.vtk",
                         half_teapot_wn,
                         bbox,
                         400,
                         200,
                         200);
  exportScalarFieldToVTK(data_dir + "/teapot_mirror_field_high.vtk",
                         mirror_teapot_wn,
                         bbox,
                         400,
                         200,
                         200);
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
