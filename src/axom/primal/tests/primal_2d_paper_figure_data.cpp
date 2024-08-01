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

TEST(primal_2d_paper_figure_data, generalized_winding)
{
  return;

  std::string data_dir =
    "E:\\Code\\winding_number\\figures\\generalized_winding\\";
  std::string curve_name = "polygon";

  axom::Array<primal::BezierCurve<double, 2>> polygon;

  // Read in the polygon, write out the curves for python
  primal::convert_from_svg(data_dir + curve_name + ".svg", polygon);
  std::ofstream polygon_curve_out(data_dir + curve_name + "_full_curves.txt");
  for(auto& edge : polygon)
  {
    polygon_curve_out << edge << std::endl;
  }

  std::ofstream edge_curve_out(data_dir + curve_name + "_edge_curves.txt");
  edge_curve_out << polygon[3] << std::endl;

  // Set up fines in whcih to write out winding number data
  std::ofstream full_winding_out(data_dir + curve_name + "_full_wn.csv");
  std::ofstream edge_winding_out(data_dir + curve_name + "_edge_wn.csv");

  // Find maximum and minimum x and y values
  double min_x = polygon[0][0][0], min_y = polygon[0][0][1];
  double max_x = polygon[0][0][0], max_y = polygon[0][0][1];

  for(auto& edge : polygon)
  {
    for(int p = 0; p <= edge.getOrder(); ++p)
    {
      min_x = axom::utilities::min(min_x, edge[p][0]);
      max_x = axom::utilities::max(max_x, edge[p][0]);

      min_y = axom::utilities::min(min_y, edge[p][1]);
      max_y = axom::utilities::max(max_y, edge[p][1]);
    }
  }

  // Put those in a bounding box, then expand it
  primal::Point<double, 2> bounds[2] = {primal::Point<double, 2> {min_x, min_y},
                                        primal::Point<double, 2> {max_x, max_y}};
  primal::BoundingBox<double, 2> bb(bounds, 2);
  bb.scale(1.25);

  // Iterate over 250 points in the x and y direction between max and min
  double npts_x = 250;
  double npts_y = 250;

  for(double x = bb.getMin()[0]; x <= bb.getMax()[0];
      x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double edge_winding_number =
        primal::winding_number(query, polygon[3], nevals);
      double full_winding_number = edge_winding_number;

      for(int i = 0; i < polygon.size(); ++i)
      {
        if(i == 3) continue;

        full_winding_number += primal::winding_number(query, polygon[i], nevals);
      }

      full_winding_out << x << "," << y << "," << full_winding_number
                       << std::endl;
      edge_winding_out << x << "," << y << "," << edge_winding_number
                       << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, closure_example)
{
  return;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\closure_example\\";
  std::string curve_name = "closure_example";

  axom::Array<primal::BezierCurve<double, 2>> closed_curve;

  // Read in the polygon, write out the curves for python
  primal::convert_from_svg(data_dir + curve_name + ".svg", closed_curve);
  for(auto& curve : closed_curve) curve.reverseOrientation();

  std::ofstream curve_out(data_dir + curve_name + "_curve_out.txt");
  std::ofstream curve_winding_out(data_dir + curve_name + "_curve_wn.csv");
  for(int i = 0; i < closed_curve.size() - 1; ++i)
  {
    curve_out << closed_curve[i] << std::endl;
  }

  std::ofstream closure_out(data_dir + curve_name + "_closure_out.txt");
  std::ofstream closure_winding_out(data_dir + curve_name + "_closure_wn.csv");
  closure_out << closed_curve[closed_curve.size() - 1] << std::endl;

  // Find maximum and minimum x and y values
  double min_x = closed_curve[0][0][0], min_y = closed_curve[0][0][1];
  double max_x = closed_curve[0][0][0], max_y = closed_curve[0][0][1];

  for(auto& edge : closed_curve)
  {
    for(int p = 0; p <= edge.getOrder(); ++p)
    {
      min_x = axom::utilities::min(min_x, edge[p][0]);
      max_x = axom::utilities::max(max_x, edge[p][0]);

      min_y = axom::utilities::min(min_y, edge[p][1]);
      max_y = axom::utilities::max(max_y, edge[p][1]);
    }
  }

  // Put those in a bounding box, then expand it
  primal::Point<double, 2> bounds[2] = {primal::Point<double, 2> {min_x, min_y},
                                        primal::Point<double, 2> {max_x, max_y}};
  primal::BoundingBox<double, 2> bb(bounds, 2);
  bb.scale(1.25);

  // Iterate over 250 points in the x and y direction between max and min
  double npts_x = 250;
  double npts_y = 250;

  for(double x = bb.getMin()[0]; x <= bb.getMax()[0];
      x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double curve_wn = 0.0;
      for(int i = 0; i < closed_curve.size() - 1; ++i)
      {
        curve_wn += primal::winding_number(query, closed_curve[i], nevals);
      }

      double closure_wn =
        primal::winding_number(query,
                               closed_curve[closed_curve.size() - 1],
                               nevals);

      curve_winding_out << x << "," << y << "," << curve_wn << std::endl;
      closure_winding_out << x << "," << y << "," << closure_wn << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, generalized_winding_curves)
{
  return;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\generalized_winding_curves\\";

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "heart";
  CurveArray shape;
  primal::convert_from_svg(data_dir + name + ".svg", shape);
  std::ofstream shape_curve_out(data_dir + name + ".txt");
  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_curve_out << shape[i] << std::endl;
  }
  BoundingBox shape_bbox = curves_bbox(shape, 1.25);
  std::ofstream shape_wn_out(data_dir + name + "_wn.csv");
  simple_grid_test(shape, shape_bbox, npts_x, npts_y, shape_wn_out);

  name = "exploded_heart2";
  CurveArray exploded_shape;
  primal::convert_from_svg(data_dir + name + ".svg", exploded_shape);
  std::ofstream exploded_shape_out(data_dir + name + ".txt");
  std::ofstream exploded_shape_wn_out(data_dir + name + "_wn.csv");
  for(int i = 0; i < exploded_shape.size(); ++i)
  {
    exploded_shape[i].reverseOrientation();
    exploded_shape_out << exploded_shape[i] << std::endl;
  }
  //BoundingBox exploded_shape_bbox = curves_bbox(exploded_shape, 1.25);
  simple_grid_test(exploded_shape, shape_bbox,
                   npts_x,
                   npts_y,
                   exploded_shape_wn_out);

  std::cout << exploded_shape << std::endl;

  name = "heart_curve_1";
  CurveArray curve_1;
  curve_1.push_back(shape[0]);
  curve_1.push_back(shape[5]);
  std::ofstream curve_1_out(data_dir + name + ".txt");
  std::ofstream curve_1_wn_out(data_dir + name + "_wn.csv");
  curve_1_out << curve_1[0] << std::endl;
  curve_1_out << curve_1[1] << std::endl;
  BoundingBox curve_1_bbox = curves_bbox(curve_1, 1.5, true);
  simple_grid_test(curve_1, curve_1_bbox, npts_x, npts_y, curve_1_wn_out);

  name = "heart_curve_2";
  CurveArray curve_2; 
  curve_2.push_back(shape[3]);
  curve_2.push_back(shape[4]);
  std::ofstream curve_2_out(data_dir + name + ".txt");
  std::ofstream curve_2_wn_out(data_dir + name + "_wn.csv");
  curve_2_out << curve_2[0] << std::endl;
  curve_2_out << curve_2[1] << std::endl;
  BoundingBox curve_2_bbox = curves_bbox(curve_2, 1.5, true);
  //curve_2_bbox.shift( (curve_1_bbox.getCentroid() - curve_2_bbox.getCentroid()) );
  simple_grid_test(curve_2, curve_2_bbox, npts_x, npts_y, curve_2_wn_out);

  name = "heart_curve_3";
  CurveArray curve_3;
  curve_3.push_back(shape[1]);
  std::ofstream curve_3_out(data_dir + name + ".txt");
  std::ofstream curve_3_wn_out(data_dir + name + "_wn.csv");
  curve_3_out << curve_3[0] << std::endl;
  BoundingBox curve_3_bbox = curves_bbox(curve_3, 1.5, true);
  //curve_3_bbox.shift((curve_1_bbox.getCentroid() - curve_3_bbox.getCentroid()));
  simple_grid_test(curve_3, curve_3_bbox, npts_x, npts_y, curve_3_wn_out);

  name = "heart_curve_4";
  CurveArray curve_4;
  curve_4.push_back(shape[2]);
  std::ofstream curve_4_out(data_dir + name + ".txt");
  std::ofstream curve_4_wn_out(data_dir + name + "_wn.csv");
  curve_4_out << curve_4[0] << std::endl;
  BoundingBox curve_4_bbox = curves_bbox(curve_4, 1.5, true);
  //curve_4_bbox.shift((curve_1_bbox.getCentroid() - curve_4_bbox.getCentroid()));
  simple_grid_test(curve_4, curve_4_bbox, npts_x, npts_y, curve_4_wn_out);
}

TEST(primal_2d_paper_figure_data, linearize_example)
{
  return;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "E:\\Code\\winding_number\\figures\\linearize_example\\";

  int npts_x = 250;
  int npts_y = 250;

  std::string name = "curve_1";
  CurveArray curve_1;
  std::ofstream curve_1_out(data_dir + name + ".txt");
  std::ofstream curve_1_wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", curve_1);
  for(int i = 0; i < curve_1.size(); ++i)
  {
    curve_1_out << curve_1[i] << std::endl;
  }
  BoundingBox largest_bbox = curves_bbox(curve_1, 1.5, true);
  simple_grid_test(curve_1, largest_bbox, npts_x, npts_y, curve_1_wn_out);

  name = "curve_2";
  CurveArray curve_2;
  std::ofstream curve_2_out(data_dir + name + ".txt");
  std::ofstream curve_2_wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", curve_2);
  for(int i = 0; i < curve_2.size(); ++i)
  {
    curve_2_out << curve_2[i] << std::endl;
  }
  simple_grid_test(curve_2, largest_bbox, npts_x, npts_y, curve_2_wn_out);

  name = "curve_3";
  CurveArray curve_3;
  std::ofstream curve_3_out(data_dir + name + ".txt");
  std::ofstream curve_3_wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", curve_3);
  for(int i = 0; i < curve_3.size(); ++i)
  {
    curve_3_out << curve_3[i] << std::endl;
  }
  BoundingBox curve_3_bbox = curves_bbox(curve_3, 1.5, true);
  simple_grid_test(curve_3, largest_bbox, npts_x, npts_y, curve_3_wn_out);

  name = "curve_4";
  CurveArray curve_4;
  std::ofstream curve_4_out(data_dir + name + ".txt");
  std::ofstream curve_4_wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", curve_4);
  for(int i = 0; i < curve_4.size(); ++i)
  {
    curve_4[i].reverseOrientation();
    curve_4_out << curve_4[i] << std::endl;
  }
  BoundingBox curve_4_bbox = curves_bbox(curve_4, 1.5, true);
  simple_grid_test(curve_4, largest_bbox, npts_x, npts_y, curve_4_wn_out);
}
//------------------------------------------------------------------------------

TEST(primal_2d_paper_figure_data, ray_casting_bad)
{
  return;

  std::cout << "Running: \"ray_casting_bad\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\ray_casting_bad\\";

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "butterfly";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream shape_cn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);
  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
  }

  // Replace curve 100 with a point
  //std::cout << shape[100] << std::endl;
  for(int i = 0; i <= shape[100].getOrder(); ++i)
  {
    shape[100][i] = shape[100][0];
  }
  //std::cout << shape[100] << std::endl;

  BoundingBox shape_bbox = curves_bbox(shape, 1.1, false);

  for(double x = shape_bbox.getMin()[0]; x <= shape_bbox.getMax()[0];
      x += (shape_bbox.getMax()[0] - shape_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = shape_bbox.getMin()[1]; y <= shape_bbox.getMax()[1];
        y += (shape_bbox.getMax()[1] - shape_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int crossing_number = primal::crossing_number(query, shape, shape_bbox);

      shape_cn_out << std::setprecision(16) << x << ',' << y << ","
                   << crossing_number << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, winding_number_good)
{
  return;
  std::cout << "Running: \"winding_number_good\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\winding_number_"
    "good\\";

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "butterfly";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream shape_wn_out(data_dir + name + "_wn.csv");
  std::ofstream zoomed_wn_out(data_dir + name + "_zoomed_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);
  BoundingBox shape_bbox = curves_bbox(shape, 1.1, false);

  BoundingBox zoomed_bbox(shape[100].boundingBox());
  double max_dim =
    axom::utilities::max(zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0],
                         zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]);
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_bbox.getMin()[0] + max_dim,
                              zoomed_bbox.getMin()[1] + max_dim});
  zoomed_bbox.scale(2.0);

  // Delete the one curve
  for(int i = 0; i <= shape[100].getOrder(); ++i) shape[100][i] = shape[100][0];

  // Print the rest to the file
  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  simple_grid_test(shape, shape_bbox, npts_x, npts_y, shape_wn_out);
  simple_grid_test(shape, zoomed_bbox, npts_x, npts_y, zoomed_wn_out);
}

TEST(primal_2d_paper_figure_data, fish_figure)
{
  return;
  std::cout << "Running: \"fish_figure\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = 
  "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\fish_figure\\";

  std::string name = "fish";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream shape_wn_out(data_dir + name + "_wn.csv");
  std::ofstream zoomed_wn_out(data_dir + name + "_zoomed_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  // Jiggle the curves, and print them all to the file
  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();

    shape[i][0][0] += (rand() % 30) - 15;
    shape[i][0][1] += (rand() % 30) - 15;

    shape[i][shape[i].getOrder()][0] += (rand() % 30) - 15;
    shape[i][shape[i].getOrder()][0] += (rand() % 30) - 15;

    shape_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox = curves_bbox(shape, 1.01, false);
  BoundingBox zoomed_bbox(shape[240].boundingBox());
  double max_dim =
    axom::utilities::max(zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0],
                         zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]);
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_bbox.getMin()[0] + max_dim,
                              zoomed_bbox.getMin()[1] + max_dim});
  zoomed_bbox.scale(12.0);

  simple_grid_test(shape, shape_bbox, 1000, 1000, shape_wn_out);
  simple_grid_test(shape, zoomed_bbox, 1000, 1000, zoomed_wn_out);
}

TEST(primal_2d_paper_figure_data, proximity_test)
{
  return;
  std::cout << "Running: \"proximity_test\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;
  using Point2D = primal::Point<double, 2>;
  using Vector2D = primal::Vector<double, 2>;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\proximity_test\\";

  std::string name = "parabola";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");

  Point2D ellipse_nodes[] = {Point2D {2.0, 0.0},
                             Point2D {2.0, 1.0},
                             Point2D {0.0, 1.0}};

  double ellipse_weights[] = {1.0, 1.0 / sqrt(2.0), 1.0};
  primal::BezierCurve<double, 2> ellipse(ellipse_nodes, ellipse_weights, 2);
  shape.push_back(ellipse);

  axom::Array<double> ts({0.09, 0.51, 0.91});
  axom::Array<double> dists({-0.1,  -1e-2, -1e-3,  -1e-4, -1e-5, -1e-6, -1e-7,
                             -1e-8, -1e-9, -1e-10, 0.1,   1e-2,  1e-3,  1e-4,
                             1e-5,  1e-6,  1e-7,   1e-8,  1e-9,  1e-10});

  std::ofstream bisection_evals_out(data_dir + "bisection_evals.csv");
  std::ofstream clipping_evals_out(data_dir + "clipping_evals.csv");
  std::ofstream boxes_evals_out(data_dir + "boxes_evals.csv");

  bisection_evals_out << std::setprecision(16) << "0";
  clipping_evals_out << std::setprecision(16) << "0";
  boxes_evals_out << std::setprecision(16) << "0";

  for(auto dist : dists)
  {
    bisection_evals_out << ',' << dist;
    clipping_evals_out << ',' << dist;
    boxes_evals_out << ',' << dist;
  }

  bisection_evals_out << std::endl;
  clipping_evals_out << std::endl;
  boxes_evals_out << std::endl;

  for(auto t0 : ts)
  {
    Point2D on_curve(ellipse.evaluate(t0));
    Vector2D tangent(ellipse.dt(t0));

    bisection_evals_out << t0;
    clipping_evals_out << t0;
    boxes_evals_out << t0;

    for(auto dist : dists)
    {
      Point2D q_int({
        -tangent[1] * dist / tangent.norm() + on_curve[0],
        tangent[0] * dist / tangent.norm() + on_curve[1],
      });

      int bisection_nevals = 0;
      int clipping_nevals = 0;
      int boxes_nevals = 0;

      double bisection_wn = primal::winding_number_bisection(q_int,
                                                             ellipse,
                                                             bisection_nevals,
                                                             1e-16,
                                                             1e-16);
      double clipping_wn = primal::winding_number_clipping(q_int,
                                                           ellipse,
                                                           clipping_nevals,
                                                           1e-16,
                                                           1e-16);
      double boxes_wn =
        primal::winding_number(q_int, ellipse, boxes_nevals, 1e-16, 1e-16);

      bisection_evals_out << ',' << bisection_nevals;
      clipping_evals_out << ',' << clipping_nevals;
      boxes_evals_out << ',' << boxes_nevals;
    }

    bisection_evals_out << std::endl;
    clipping_evals_out << std::endl;
    boxes_evals_out << std::endl;
  }
}

TEST(primal_2d_paper_figure_data, linearization_test)
{
  return;

  std::cout << "Running: \"linearization_test\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\linearization_test\\";

  std::string name = "bubble";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream misclassify_vals_out(data_dir + "misclassify_vals.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox = curves_bbox(shape, 1.05, false);

  std::ofstream shape_wn_out(data_dir + name + "_wnn.csv");
  simple_grid_test(shape, shape_bbox, 1000, 1000, shape_wn_out);

  std::ofstream shape_grid_out(data_dir + name + "_grid.txt");
  simple_grid_test(shape, shape_bbox, 2, 2, shape_grid_out);

  constexpr double tol = 1e-10;
  int npts = 1e5;
  int max_refinement = 25;
  int dummy_int = 0;

  for(int n = 4; n <= 0; ++n)
  {
    std::cout << std::endl << "Doing refinement " << n << std::endl;
    std::ofstream linearized_shape_out(data_dir + name + std::to_string(n) +
                                       ".txt");

    // Construct the linear approximation with the given level of refinement
    CurveArray linearized_shape;
    for(int i = 0; i < shape.size(); ++i)
    {
      for(double j = 0; j < n; ++j)
      {
        primal::BezierCurve<double, 2> line(1);
        line[0] = shape[i].evaluate(j / n);
        line[1] = shape[i].evaluate((j + 1) / n);

        linearized_shape.push_back(line);
      }
    }

    std::ofstream linearized_shape_wn_out(data_dir + name +
                                          "_linearized_wn.csv");
    simple_grid_test(linearized_shape, shape_bbox, 500, 500, linearized_shape_wn_out);

    for(int i = 0; i < linearized_shape.size(); ++i)
    {
      linearized_shape_out << linearized_shape[i] << std::endl;
    }

    std::ofstream misclassifications_out(data_dir + "misclassifications_" +
                                         std::to_string(n) + ".csv");
    int misclassifications = 0;

    for(int m = 0; m < npts; ++m)
    {
      primal::printLoadingBar(m, npts);
      double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                              shape_bbox.getMax()[0]);
      double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                              shape_bbox.getMax()[1]);

      primal::Point<double, 2> query({x, y});

      double linearized_wn = 0;
      for(int i = 0; i < linearized_shape.size(); ++i)
      {
        linearized_wn += primal::winding_number(query,
                                                linearized_shape[i],
                                                dummy_int,
                                                1e-16,
                                                1e-16);
      }
      bool linearized_is_in = std::lround(linearized_wn) != 0;

      double true_wn = 0;
      for(int i = 0; i < shape.size(); ++i)
      {
        true_wn +=
          primal::winding_number(query, shape[i], dummy_int, 1e-16, 1e-16);
      }
      bool true_is_in = std::lround(true_wn) != 0;

      if(linearized_is_in != true_is_in)
      {
        misclassifications++;
      }
    }

    misclassifications_out << misclassifications << std::endl;
  }
}

TEST(primal_2d_paper_figure_data, tolerance_test)
{
  return;

  std::cout << "Running : \"tolerance test\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using SegmentArray = axom::Array<primal::Segment<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\tolerance_test\\";

  std::string name = "tangle";
  CurveArray shape;
  primal::Point<double, 2> tangle_nodes[] = {
    primal::Point<double, 2> {1.000, 0.000},
    primal::Point<double, 2> {-3.000, 0.000},
    primal::Point<double, 2> {-3.225, 3.175},
    primal::Point<double, 2> {1.008, 3.300},
    primal::Point<double, 2> {5.223, -3.836},
    primal::Point<double, 2> {-3.273, -4.321},
    primal::Point<double, 2> {-5.485, 2.728},
    primal::Point<double, 2> {2.681, 3.014},
    primal::Point<double, 2> {3.000, 0.000},
    primal::Point<double, 2> {-2.000, 0.000}};

  double tangle_weights[] = {1.0, 1.1, 0.7, 1.1, 1.0, 1.2, 1.1, 1.0, 1.0, 1.0};
  primal::BezierCurve<double, 2> tangle_curve(tangle_nodes, tangle_weights, 9);
  //shape.push_back(tangle_curve);
  //tangle_curve.reverseOrientation();

  primal::Point<double, 2> closure_nodes[] = {
    primal::Point<double, 2> {-2.0, 0.0},
    primal::Point<double, 2> {-0.5, -3},
    primal::Point<double, 2> {1.0, 0.0}

  };
  primal::BezierCurve<double, 2> closure_curve(closure_nodes, 2);
  shape.push_back(closure_curve);

  primal::Point<double, 2> closure_closure_nodes[] = {
    primal::Point<double, 2> {1.0, 0.0},
    primal::Point<double, 2> {-2.0, 0.0}};

  primal::BezierCurve<double, 2> closure_closure_curve(closure_closure_nodes, 2);
  shape.push_back(closure_closure_curve);

  std::ofstream curves_out(data_dir + name + ".txt");

  //primal::convert_from_svg(data_dir + name + ".svg", shape);

  srand(100);

  double shape_radius = 1;
  BoundingBox shape_bbox;
  primal::Point<double, 2> center({0.2, 0.3});
  shape_bbox.addPoint(primal::Point<double, 2> {center[0] + shape_radius,
                                                center[1] + shape_radius});
  shape_bbox.addPoint(primal::Point<double, 2> {center[0] - shape_radius,
                                                center[1] - shape_radius});

  // Get the largest axis of the bounding box
  double max_dim =
    axom::utilities::max(shape_bbox.getMax()[0] - shape_bbox.getMin()[0],
                         shape_bbox.getMax()[1] - shape_bbox.getMin()[1]);

  for(int i = 0; i < shape.size(); ++i)
  {
    std::cout << shape[i] << std::endl;
    for(int j = 0; j <= shape[i].getOrder(); ++j)
    {
      shape[i][j][0] = 0.25 * shape[i][j][0] + 0.62;
      shape[i][j][1] = 0.25 * shape[i][j][1] + 0.55;
    }

    curves_out << shape[i] << std::endl;

    // Make a fuller version of the shape that has three times as many control nodes
    //  primal::BezierCurve<double, 2> fuller_curve(3 * shape[i].getOrder());
    //  for (int j = 0; j <= 3 * shape[i].getOrder(); ++j)
    //  {
    //	fuller_curve[j] = shape[i][j / 3];
    //}
  }

  shape_bbox.clear();

  shape_bbox.addPoint(primal::Point<double, 2> {-0.1, -0.1});
  shape_bbox.addPoint(primal::Point<double, 2> {1.1, 1.1});

  std::ofstream shape_grid_out(data_dir + name + "_grid.txt");
  simple_grid_test(shape, shape_bbox, 500, 500, shape_grid_out);

  //return;

  std::ofstream tolerance_wn_timing_out(data_dir +
                                        "tolerance_wn_timing_out_simple.csv");

  //constexpr double EPS = 0.0;
  int npts = 1e5;
  const int max_refinement = 25;
  int dummy_int = 0;
  bool dummy_bool = false;
  int idx;

  const int NUM_CURVES = 20;

  // Make an array of duration objects
  double avg_depth[15];
  int misclassifications[15];
  int nevals[15];
  for(int i = 0; i < 15; ++i)
  {
    misclassifications[i] = 0;
    nevals[i] = 0;
    avg_depth[i] = 0;
  }

  //int idx2 = 0;
  //for(int i = 2; i <= 0; ++i)
  //{
  //  double the_tol = std::pow(10, -i);

  //  std::cout << "===================" << the_tol
  //            << "===================" << std::endl;

  //  for(int m = 0; m < npts; ++m)
  //  {
  //    double x = axom::utilities::random_real(shape_bbox.getMin()[0],
  //                                            shape_bbox.getMax()[0]);
  //    double y = axom::utilities::random_real(shape_bbox.getMin()[1],
  //                                            shape_bbox.getMax()[1]);
  //    primal::Point<double, 2> query({x, y});
  //    primal::Polygon<double, 2> temp_approxogon(20);

  //    double adaptive_wn = primal::winding_number_approxogon(query,
  //                                                           shape,
  //                                                           temp_approxogon,
  //                                                           0.0,
  //                                                           the_tol);
  //  }

  //  idx2++;
  //}

  //return;
  // Loop over the number of points
  for(int m = 0; m < npts; ++m)
  {
    primal::printLoadingBar(m, npts);
    double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                            shape_bbox.getMax()[0]);
    double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                            shape_bbox.getMax()[1]);

    primal::Point<double, 2> query({x, y});

    // Loop over the tolerances, 1e-2 to 1e-16
    idx = 0;
    for(int i = 2; i <= 16; ++i)
    {
      double the_tol = std::pow(10, -i);

      //std::cout << the_tol << std::endl;

      double adaptive_wn = 0.0;
      primal::Polygon<double, 2> temp_approxogon(20);

      double true_wn = primal::winding_number(query, shape, 0.0);
      int this_depth = 0;
      adaptive_wn = primal::winding_number_approxogon(query,
                                                      shape,
                                                      temp_approxogon,
                                                      nevals[idx],
                                                      this_depth,
                                                      0.0,
                                                      the_tol);

      if((std::lround(adaptive_wn) != 0) != (std::lround(true_wn) != 0))
      {
        misclassifications[idx]++;
      }
      avg_depth[idx] += static_cast<double>(this_depth) / npts;

      idx++;
    }
  }

  idx = 0;
  for(int i = 2; i <= 16; ++i)
  {
    double the_tol = std::pow(10, -i);

    tolerance_wn_timing_out << the_tol << "," << avg_depth[idx] << ","
                            << misclassifications[idx] << "," << nevals[idx]
                            << std::endl;
    idx++;
  }
}

TEST(primal_2d_paper_figure_data, timing_test)
{
  return;

  std::cout << "Running: \"timing_test\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using SegmentArray = axom::Array<primal::Segment<double, 2>>;
  using PolygonArray = axom::Array<primal::Polygon<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string prefix = "run_3\\";

  std::string data_dir = "E:\\Code\\winding_number\\figures\\timing_test\\";

  std::string name = "bubble";
  CurveArray shape;
  std::ofstream curves_out(data_dir + name + ".txt");
  std::ofstream adaptive_wn_timing_out(data_dir + prefix +
                                       "method0_wn_timing_out.csv");
  std::ofstream memoized_wn_timing_out(data_dir + prefix +
                                       "method00_wn_timing_out.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  srand(10101);

  for(int i = 0; i < shape.size(); ++i)
  {
    //curves_out << shape[i] << std::endl;

    primal::BezierCurve<double, 2> fuller_curve(3 * shape[i].getOrder());
    for(int j = 0; j <= 3 * shape[i].getOrder(); ++j)
    {
      fuller_curve[j] = shape[i][j / 3];
    }
    fuller_curve.makeRational();
    shape[i].makeRational();
    //shape[i] = fuller_curve;
    //shape_fuller.push_back(fuller_curve);

    curves_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox = curves_bbox(shape, 1.05, false);

  std::ofstream shape_grid_out(data_dir + name + "_grid.txt");
  simple_grid_test(shape, shape_bbox, 2, 2, shape_grid_out);

  constexpr double tol = 1e-10;
  int npts = 1e5;
  const int max_refinement = 25;
  int dummy_int = 0;
  bool dummy_bool = false;

  const int NUM_CURVES = 20;

  // Do the preprocessing to get the linearization for each level of refinement
  //SegmentArray bag_of_linearizations[max_refinement];
  //SegmentArray batches_of_linearizations[max_refinement][NUM_CURVES];  // There are 20 curves in the bubble
  PolygonArray bag_of_polygons[max_refinement];
  std::unordered_map<std::pair<int, int>, primal::Segment<double, 2>, primal::detail::PairHash>
    segment_hashes[max_refinement];

  // Array for low order refinement
  std::vector<std::vector<primal::detail::BezierCurveMemo<double>>>
    array_memos[NUM_CURVES];

  // Hash for higher order refinement
  axom::FlatMap<std::pair<int, int>,
                primal::detail::BezierCurveMemo<double>,
                primal::detail::PairHash>
    hash_memos[NUM_CURVES];

  std::ofstream method1_wn_timing_out[max_refinement];  // Big bag of winding numbers. Completely naive
  std::ofstream method2_wn_timing_out[max_refinement];  // Semi-sophisticated. If you're outside the box, do one calculation. If you're inside the box, do linearization
  std::ofstream method3_wn_timing_out[max_refinement];  // Sophisticated. Make a polygon with a fixed linearization, then do polygon_wn - closure
  std::ofstream method4_wn_timing_out[max_refinement];  // Method 3 with no pre-processing

  std::chrono::duration<double> preprocessing_times[max_refinement][4];
  std::ofstream preprocessing_times_out =
    std::ofstream(data_dir + prefix + "preprocessing_times.csv");

  // Preprocessing for method 00
  constexpr int LEVEL = 3;
  auto start = std::chrono::high_resolution_clock::now();
  for(int ci = 0; ci < NUM_CURVES; ++ci)
  {
    array_memos[ci].resize(LEVEL);
    for(int i = 0; i < LEVEL; ++i) array_memos[ci][i].resize(1 << i);

    array_memos[ci][0][0] = {
      //shape[ci].isLinear(),
      primal::is_convex(primal::Polygon<double, 2>(shape[ci].getControlPoints()),
                        primal::PRIMAL_TINY),
      shape[ci]};

    for(int i = 0; i < LEVEL - 1; ++i)
    {
      for(int j = 0; j < (1 << i); ++j)
      {
        primal::BezierCurve<double, 2> c1, c2;
        array_memos[ci][i][j].curve.split(0.5, c1, c2);

        array_memos[ci][i + 1][2 * j] = {
          //c1.isLinear(),
          primal::is_convex(primal::Polygon<double, 2>(c1.getControlPoints()),
                            primal::PRIMAL_TINY),
          c1};

        array_memos[ci][i + 1][2 * j + 1] = {
          //c2.isLinear(),
          primal::is_convex(primal::Polygon<double, 2>(c2.getControlPoints()),
                            primal::PRIMAL_TINY),
          c2};
      }
    }

    //hash_memos.push_back(axom::FlatMap<std::pair<int, int>,
    //                                   primal::detail::BezierCurveMemo<double>,
    //                                   primal::detail::PairHash>());
  }
  auto end = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < max_refinement; ++i)
  {
    preprocessing_times[i][3] = end - start;
  }

  /*std::ofstream shape_memoized_grid_out(data_dir + name + "_memoized_grid.txt");
  simple_grid_test_memoized(shape_marray,
                            shape_bbox,
                            200,
                            200,
                            shape_memoized_grid_out);*/

  // Set up the linearization
  for(int ref = 0; ref < max_refinement; ++ref)
  {
    method1_wn_timing_out[ref] =
      std::ofstream(data_dir + prefix + "method1_wn_timing_out_" +
                    std::to_string(ref + 1) + ".csv");
    method2_wn_timing_out[ref] =
      std::ofstream(data_dir + prefix + "method2_wn_timing_out_" +
                    std::to_string(ref + 1) + ".csv");
    method3_wn_timing_out[ref] =
      std::ofstream(data_dir + prefix + "method3_wn_timing_out_" +
                    std::to_string(ref + 1) + ".csv");
    method4_wn_timing_out[ref] =
      std::ofstream(data_dir + prefix + "method4_wn_timing_out_" +
                    std::to_string(ref + 1) + ".csv");

    // Big boy preprocessing
    start = std::chrono::high_resolution_clock::now();
    for(int ci = 0; ci < shape.size(); ++ci)
    {
      for(int i = 0; i < (ref + 1); ++i)
      {
        primal::BezierCurve<double, 2> curve = shape[ci];
        primal::Point<double, 2> point = curve.evaluate(1.0 * i / (ref + 1.0));
        primal::Point<double, 2> next_point =
          curve.evaluate(1.0 * (i + 1.0) / (ref + 1.0));
        primal::Segment<double, 2> line(point, next_point);
        std::pair<int, int> the_pair = std::make_pair(ci, i);
        segment_hashes[ref][the_pair] = line;
      }
    }
    end = std::chrono::high_resolution_clock::now();
    preprocessing_times[ref][0] = end - start;

    // Preprocessing for method 1
    //start = std::chrono::high_resolution_clock::now();
    //for(int ci = 0; ci < NUM_CURVES; ++ci)
    //{
    //  for(int i = 0; i < (ref + 1); ++i)
    //  {
    //    primal::Segment<double, 2> line(
    //      shape[ci].evaluate(1.0 * i / (ref + 1.0)),
    //      shape[ci].evaluate(1.0 * (i + 1.0) / (ref + 1.0)));
    //    bag_of_linearizations[ref].push_back(line);
    //  }
    //}
    //end = std::chrono::high_resolution_clock::now();
    //preprocessing_times[ref][0] = end - start;

    // Preprocessing for method 2
    start = std::chrono::high_resolution_clock::now();
    //for(int ci = 0; ci < NUM_CURVES; ci++)
    //{
    //  for(int i = 0; i < (ref + 1); ++i)
    //  {
    //    primal::Segment<double, 2> line(
    //      shape[ci].evaluate(1.0 * i / (ref + 1.0)),
    //      shape[ci].evaluate(1.0 * (i + 1.0) / (ref + 1.0)));

    //    batches_of_linearizations[ref][ci].push_back(line);
    //  }
    //}
    end = std::chrono::high_resolution_clock::now();
    preprocessing_times[ref][1] = end - start;

    // Preprocessing for method 3
    start = std::chrono::high_resolution_clock::now();
    //for(int ci = 0; ci < NUM_CURVES; ++ci)
    //{
    //  primal::Polygon<double, 2> poly(ref + 1);

    //  for(int i = 0; i < (ref + 1); ++i)
    //  {
    //    poly.addVertex(shape[ci].evaluate(1.0 * i / (ref + 1.0)));
    //  }

    //  poly.addVertex(shape[ci].evaluate(1.0));
    //  bag_of_polygons[ref].push_back(poly);
    //}
    end = std::chrono::high_resolution_clock::now();
    preprocessing_times[ref][2] = end - start;

    // Construct the linear approximation out of BezierCurves for plotting
    CurveArray linearized_curves;
    for(int ci = 0; ci < shape.size(); ++ci)
    {
      for(int i = 0; i < (ref + 1); ++i)
      {
        primal::BezierCurve<double, 2> curve(1);
        curve[0] = shape[ci].evaluate(1.0 * i / (ref + 1.0));
        curve[1] = shape[ci].evaluate(1.0 * (i + 1.0) / (ref + 1.0));
        linearized_curves.push_back(curve);
      }
    }

    std::ofstream fixed_wn_shape_out(data_dir + "fixed_wn_shape_out_" +
                                     std::to_string(ref + 1) + ".txt");

    for(auto& line : linearized_curves)
    {
      fixed_wn_shape_out << line << std::endl;
    }
  }

  // Write the preprocessing times to the file
  for(int ref = 0; ref < max_refinement; ++ref)
  {
    preprocessing_times_out << ref + 1 << ','
                            << preprocessing_times[ref][0].count() << ','
                            << preprocessing_times[ref][1].count() << ','
                            << preprocessing_times[ref][2].count() << ','
                            << preprocessing_times[ref][3].count() << std::endl;
  }

  // Loop over the number of points
  for(int m = 0; m < npts; ++m)
  {
    primal::printLoadingBar(m, npts);
    double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                            shape_bbox.getMax()[0]);
    double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                            shape_bbox.getMax()[1]);

    primal::Point<double, 2> query({x, y});
    //184.194, 113.01
    //query[0] = 184.194;
    //query[1] = 113.01;

    double adaptive_wn = 0.0, memoized_wn = 0.0;
    primal::Polygon<double, 2> temp_approxogon(20);
    std::stack<std::pair<int, int>> curve_stack;

    auto start = std::chrono::high_resolution_clock::now();
    adaptive_wn = primal::winding_number_approxogon(query,
                                                    shape,
                                                    temp_approxogon,
                                                    dummy_int,
                                                    dummy_int,
                                                    tol);
    auto end = std::chrono::high_resolution_clock::now();

    adaptive_wn_timing_out << x << "," << y << ","
                           << std::chrono::duration<double>(end - start).count()
                           << std::endl;

    for(int ref = 0; ref < max_refinement; ++ref)
    {
      double wn1 = 0.0, wn3 = 0.0, wn4 = 0.0;

      /* Do method 2: bounding box + small bag version */
      start = std::chrono::high_resolution_clock::now();
      double wn2 = 0.0;
      for(int i = 0; i < NUM_CURVES; ++i)
      {
        if(!shape[i].boundingBox().contains(query))
        {
          wn2 -= primal::winding_number(
            query,
            primal::Segment<double, 2>(shape[i][shape[i].getOrder()], shape[i][0]),
            tol);
        }
        else
        {
          for(int j = 0; j < ref + 1; ++j)
          {
            auto the_pair = std::make_pair(i, j);
            wn2 +=
              primal::winding_number(query, segment_hashes[ref][the_pair], tol);
          }
        }
      }
      end = std::chrono::high_resolution_clock::now();

      method2_wn_timing_out[ref]
        << x << ',' << y << ','
        << (std::chrono::duration<double>(end - start)).count() << std::endl;

      /* Do method 1: big bag o winding numbers version */
      auto start = std::chrono::high_resolution_clock::now();
      //for(int i = 0; i < bag_of_linearizations[ref].size(); ++i)
      //{
      //wn1 += primal::winding_number(query, bag_of_linearizations[ref][i], tol);
      //}
      auto end = std::chrono::high_resolution_clock::now();

      method1_wn_timing_out[ref]
        << x << ',' << y << ','
        << (std::chrono::duration<double>(end - start)).count() << std::endl;

      /* Do method 3: polygon version */
      start = std::chrono::high_resolution_clock::now();
      for(int i = 0; i < NUM_CURVES; ++i)
      {
        /*if(!shape[i].boundingBox().contains(query))
        {
          wn3 -= primal::winding_number(
            query,
            primal::Segment<double, 2>(shape[i][shape[i].getOrder()], shape[i][0]),
            tol);
        }
        else
        {
          wn3 +=
            primal::detail::approxogon_winding_number(query,
                                                      bag_of_polygons[ref][i],
                                                      tol);
        }*/
      }
      end = std::chrono::high_resolution_clock::now();

      method3_wn_timing_out[ref]
        << x << ',' << y << ','
        << (std::chrono::duration<double>(end - start)).count() << std::endl;

      ///* Do method 4: polygon with no pre-processing*/
      //primal::Polygon<double, 2> this_poly(ref + 1);

      //start = std::chrono::high_resolution_clock::now();
      //for(int i = 0; i < NUM_CURVES; ++i)
      //{
      //  if(!shape[i].boundingBox().contains(query))
      //  {
      //    wn4 -= primal::winding_number(
      //      query,
      //      primal::Segment<double, 2>(shape[i][shape[i].getOrder()], shape[i][0]),
      //      tol);
      //  }
      //  else
      //  {
      //    this_poly.clear();

      //    for(int j = 0; j < (ref + 1); ++j)
      //    {
      //      this_poly.addVertex(shape[i].evaluate(1.0 * j / (ref + 1.0)));
      //    }
      //    this_poly.addVertex(shape[i].evaluate(1.0));

      //    wn4 += primal::detail::approxogon_winding_number(query, this_poly, tol);
      //  }
      //}
      //end = std::chrono::high_resolution_clock::now();

      //method4_wn_timing_out[ref]
      //  << x << ',' << y << ','
      //  << std::chrono::duration<double>(end - start).count() << std::endl;
    }

    start = std::chrono::high_resolution_clock::now();
    //memoized_wn = primal::winding_number_approxogon_memoized(query,
    //                                                         array_memos,
    //                                                         hash_memos,
    //                                                         temp_approxogon,
    //                                                         curve_stack,
    //                                                         tol);
    for(int i = 0; i < NUM_CURVES; ++i)
    {
      if(!shape[i].boundingBox().contains(query))
      {
        memoized_wn -= primal::winding_number(
          query,
          primal::Segment<double, 2>(shape[i][shape[i].getOrder()], shape[i][0]),
          tol);
      }
      else
      {
        primal::detail::winding_number_adaptive_linear_memoized(query,
                                                                array_memos[i],
                                                                hash_memos[i],
                                                                0,
                                                                tol,
                                                                temp_approxogon,
                                                                curve_stack,
                                                                memoized_wn);
        temp_approxogon.addVertex(shape[i][shape[i].getOrder()]);
        memoized_wn +=
          primal::detail::approxogon_winding_number(query, temp_approxogon, 0);
        temp_approxogon.clear();
      }
    }
    end = std::chrono::high_resolution_clock::now();

    memoized_wn_timing_out << x << "," << y << ","
                           << std::chrono::duration<double>(end - start).count()
                           << std::endl;
  }
}

TEST(primal_2d_paper_figure_data, rose_figure)
{
  return;
  bool collectGridData = false;
  bool collectRandomData = true;

  std::cout << "Running: \"rose_figure\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\rose_figure\\";

  std::string name = "rose";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  // Jiggle the curves, and print them all to the file
  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox = curves_bbox(shape, 1.05, false);

  if(collectGridData)
  {
    std::ofstream grid_boxes_nevals_out(data_dir + "boxes_grid_nevals.csv");
    std::ofstream grid_clipping_nevals_out(data_dir +
                                           "clipping_grid_nevals.csv");
    std::ofstream grid_bisection_nevals_out(data_dir +
                                            "bisection_grid_nevals.csv");

    std::ofstream grid_boxes_wn_out(data_dir + "boxes_grid_wn.csv");
    std::ofstream grid_clipping_wn_out(data_dir + "clipping_grid_wn.csv");
    std::ofstream grid_bisection_wn_out(data_dir + "bisection_grid_wn.csv");

    int npts_x = 250;
    int npts_y = 250;

    for(double x = shape_bbox.getMin()[0]; x <= shape_bbox.getMax()[0];
        x += (shape_bbox.getMax()[0] - shape_bbox.getMin()[0]) / npts_x)
    {
      std::cout << x << std::endl;

      for(double y = shape_bbox.getMin()[1]; y <= shape_bbox.getMax()[1];
          y += (shape_bbox.getMax()[1] - shape_bbox.getMin()[1]) / npts_y)
      {
        primal::Point<double, 2> query({x, y});

        int nevals_boxes = 0;
        int nevals_bisection = 0;
        int nevals_clipping = 0;

        double wn_boxes = 0.0;
        double wn_bisection = 0.0;
        double wn_clipping = 0.0;

        for(int i = 0; i < shape.size(); ++i)
        {
          wn_boxes +=
            primal::winding_number(query, shape[i], nevals_boxes, 1e-16, 1e-16);
          wn_clipping += primal::winding_number_clipping(query,
                                                         shape[i],
                                                         nevals_clipping,
                                                         1e-16,
                                                         1e-16);
          wn_bisection += primal::winding_number_bisection(query,
                                                           shape[i],
                                                           nevals_bisection,
                                                           1e-16,
                                                           1e-16);
        }

        grid_boxes_nevals_out << x << "," << y << "," << nevals_boxes
                              << std::endl;
        grid_clipping_nevals_out << x << "," << y << "," << nevals_clipping
                                 << std::endl;
        grid_bisection_nevals_out << x << "," << y << "," << nevals_bisection
                                  << std::endl;

        grid_boxes_wn_out << x << "," << y << "," << wn_boxes << std::endl;
        grid_clipping_wn_out << x << "," << y << "," << wn_clipping << std::endl;
        grid_bisection_wn_out << x << "," << y << "," << wn_bisection
                              << std::endl;
      }
    }
  }

  if(collectRandomData)
  {
    std::ofstream random_boxes_nevals_out(data_dir + "boxes_random_nevals.csv");
    std::ofstream random_clipping_nevals_out(data_dir +
                                             "clipping_random_nevals.csv");
    std::ofstream random_bisection_nevals_out(data_dir +
                                              "bisection_random_nevals.csv");

    int npts = 250000;
    constexpr double tol = 1e-10;

    for(int n = 0; n < npts; ++n)
    {
      primal::printLoadingBar(n, npts);

      double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                              shape_bbox.getMax()[0]);
      double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                              shape_bbox.getMax()[1]);

      primal::Point<double, 2> query({x, y});

      int nevals_boxes = 0;
      int nevals_bisection = 0;
      int nevals_clipping = 0;

      double wn_boxes = 0.0;
      double wn_bisection = 0.0;
      double wn_clipping = 0.0;
      for(int i = 0; i < shape.size(); ++i)
      {
        //desmos_print(shape[i]);
        //desmos_print(query);
        wn_boxes +=
          primal::winding_number(query, shape[i], nevals_boxes, tol, tol);
        wn_clipping += primal::winding_number_clipping(query,
                                                       shape[i],
                                                       nevals_clipping,
                                                       tol,
                                                       tol);
        wn_bisection += primal::winding_number_bisection(query,
                                                         shape[i],
                                                         nevals_bisection,
                                                         tol,
                                                         tol);
      }

      if(!axom::utilities::isNearlyEqual(wn_bisection, wn_boxes, 1e-5))
      {
        std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                  << std::endl;
        std::cout << "bisection machine broke" << std::endl;
      }

      if(!axom::utilities::isNearlyEqual(wn_boxes, wn_clipping, 1e-5))
      {
        std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                  << std::endl;
        std::cout << "clipping machine broke" << std::endl;
      }

      random_boxes_nevals_out << x << "," << y << "," << nevals_boxes
                              << std::endl;
      random_clipping_nevals_out << x << "," << y << "," << nevals_clipping
                                 << std::endl;
      random_bisection_nevals_out << x << "," << y << "," << nevals_bisection
                                  << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, spiral_figure)
{
  return;

  std::cout << "Running: \"spiral_figure\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;
  using Point2D = primal::Point<double, 2>;

  std::string data_dir = 
  "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\spiral_figure\\";

  std::string name = "spiral";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");

  Point2D spiral_nodes[] = {Point2D {-0.824, -0.927},
                            Point2D {-3.344, -0.802},
                            Point2D {-2.085, 3.242},
                            Point2D {1.794, 4.861},
                            Point2D {2.692, -4.142},
                            Point2D {-2.313, -3.279},
                            Point2D {-6.041, 2.037},
                            Point2D {1.98, 2.444},
                            Point2D {-0.966, -0.939},
                            Point2D {-1.034, 0.113}};

  double spiral_weights[] = {1.0, 1.2, 1.2, 0.9, 1.1, 0.9, 1.0, 1.0, 1.1, 1.0};
  primal::BezierCurve<double, 2> spiral_curve(spiral_nodes, spiral_weights, 9);
  shape.push_back(spiral_curve);

  //Point2D spiral_nodes[] = {Point2D {2.0, 0.0},
  //                          Point2D {2.0, 1.0},
  //                          Point2D {0.0, 1.0}};

  //double spiral_weights[] = {1.0, 1.0 / sqrt(2.0), 1.0};
  //primal::BezierCurve<double, 2> spiral_curve(spiral_nodes, spiral_weights, 2);
  //shape.push_back(spiral_curve);

  // Jiggle the curves, and print them all to the file
  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox;  //  = curves_bbox(shape, 1.05, true);

  shape_bbox.addPoint(spiral_nodes[9]);
  shape_bbox.addPoint(
    Point2D {spiral_nodes[9][0] + 1.5, spiral_nodes[9][1] + 1.5});
  shape_bbox.addPoint(
    Point2D {spiral_nodes[9][0] - 1.5, spiral_nodes[9][1] - 1.5});
  std::cout << shape_bbox << std::endl;

  std::ofstream wn_out(data_dir + name + "_wn.csv");
  simple_grid_test(shape, shape_bbox, 1000, 1000, wn_out);
  return;

  std::ofstream grid_nevals_out(data_dir + "grid_nevals.csv");
  int npts_x = 3;
  int npts_y = 3;

  for(int xi = 0; xi < 2; ++xi)
  {
    for(int yi = 0; yi < 2; ++yi)
    {
      double x = shape_bbox.getMin()[0] +
        (shape_bbox.getMax()[0] - shape_bbox.getMin()[0]) * xi;
      double y = shape_bbox.getMin()[1] +
        (shape_bbox.getMax()[1] - shape_bbox.getMin()[1]) * yi;
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < shape.size(); ++i)
      {
        primal::winding_number(query, shape[i], nevals, 1e-16, 1e-16);
      }

      grid_nevals_out << x << "," << y << "," << nevals << std::endl;
    }
  }

  std::ofstream random_boxes_nevals_out(data_dir + "boxes_random_nevals.csv");
  std::ofstream random_clipping_nevals_out(data_dir +
                                           "clipping_random_nevals.csv");
  std::ofstream random_bisection_nevals_out(data_dir +
                                            "bisection_random_nevals.csv");

  int npts = 250000;
  constexpr double tol = 1e-10;

  for(int n = 0; n < npts; ++n)
  {
    primal::printLoadingBar(n, npts);

    double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                            shape_bbox.getMax()[0]);
    double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                            shape_bbox.getMax()[1]);

    primal::Point<double, 2> query({x, y});

    int nevals_boxes = 0;
    int nevals_bisection = 0;
    int nevals_clipping = 0;

    double wn_boxes = 0.0;
    double wn_bisection = 0.0;
    double wn_clipping = 0.0;
    for(int i = 0; i < shape.size(); ++i)
    {
      //desmos_print(shape[i]);
      //desmos_print(query);
      wn_boxes += primal::winding_number(query, shape[i], nevals_boxes, tol, tol);
      wn_clipping +=
        primal::winding_number_clipping(query, shape[i], nevals_clipping, tol, tol);
      wn_bisection += primal::winding_number_bisection(query,
                                                       shape[i],
                                                       nevals_bisection,
                                                       tol,
                                                       tol);
    }

    if(!axom::utilities::isNearlyEqual(wn_bisection, wn_boxes, 1e-5))
    {
      std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                << std::endl;
      std::cout << "bisection machine broke" << std::endl;
    }

    if(!axom::utilities::isNearlyEqual(wn_boxes, wn_clipping, 1e-5))
    {
      std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                << std::endl;
      std::cout << "clipping machine broke" << std::endl;
    }

    random_boxes_nevals_out << x << "," << y << "," << nevals_boxes << std::endl;
    random_clipping_nevals_out << x << "," << y << "," << nevals_clipping
                               << std::endl;
    random_bisection_nevals_out << x << "," << y << "," << nevals_bisection
                                << std::endl;
  }
}

TEST(primal_2d_paper_figure_data, tangle_figure)
{
  return;

  std::cout << "Running: \"tangle_figure\"" << std::endl;
  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;
  using Point2D = primal::Point<double, 2>;

  std::string data_dir = 
  "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\tangle_figure\\";

  std::string name = "tangle";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");

  Point2D tangle_nodes[] = {Point2D {1.000, 0.000},
                            Point2D {-3.000, 0.000},
                            Point2D {-3.225, 3.175},
                            Point2D {1.008, 3.300},
                            Point2D {5.223, -3.836},
                            Point2D {-3.273, -4.321},
                            Point2D {-5.485, 2.728},
                            Point2D {2.681, 3.014},
                            Point2D {3.000, 0.000},
                            Point2D {-2.000, 0.000}};

  double tangle_weights[] = {1.0, 1.1, 0.7, 1.1, 1.0, 1.2, 1.1, 1.0, 1.0, 1.0};
  primal::BezierCurve<double, 2> tangle_curve(tangle_nodes, tangle_weights, 9);
  shape.push_back(tangle_curve);

  // Jiggle the curves, and print them all to the file
  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  BoundingBox shape_bbox;
  Point2D center({-0.33, 0.2});
  shape_bbox.addPoint(Point2D {center[0] + 2, center[1] + 2});
  shape_bbox.addPoint(Point2D {center[0] - 2, center[1] - 2});

  std::ofstream wn_out(data_dir + name + "_wn.csv");
  simple_grid_test(shape, shape_bbox, 1000, 1000, wn_out);

  return;

  std::ofstream grid_nevals_out(data_dir + "grid_nevals.csv");
  int npts_x = 3;
  int npts_y = 3;

  for(int xi = 0; xi < 2; ++xi)
  {
    for(int yi = 0; yi < 2; ++yi)
    {
      double x = shape_bbox.getMin()[0] +
        (shape_bbox.getMax()[0] - shape_bbox.getMin()[0]) * xi;
      double y = shape_bbox.getMin()[1] +
        (shape_bbox.getMax()[1] - shape_bbox.getMin()[1]) * yi;
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < shape.size(); ++i)
      {
        primal::winding_number(query, shape[i], nevals, 1e-16, 1e-16);
      }

      grid_nevals_out << x << "," << y << "," << nevals << std::endl;
    }
  }

  std::ofstream random_boxes_nevals_out(data_dir + "boxes_random_nevals.csv");
  std::ofstream random_clipping_nevals_out(data_dir +
                                           "clipping_random_nevals.csv");
  std::ofstream random_bisection_nevals_out(data_dir +
                                            "bisection_random_nevals.csv");

  int npts = 1e5 + 1e5 / 2;
  constexpr double tol = 1e-10;

  for(int n = 0; n < npts; ++n)
  {
    primal::printLoadingBar(n, npts);

    double x = axom::utilities::random_real(shape_bbox.getMin()[0],
                                            shape_bbox.getMax()[0]);
    double y = axom::utilities::random_real(shape_bbox.getMin()[1],
                                            shape_bbox.getMax()[1]);

    primal::Point<double, 2> query({x, y});

    int nevals_boxes = 0;
    int nevals_bisection = 0;
    int nevals_clipping = 0;

    double wn_boxes = 0.0;
    double wn_bisection = 0.0;
    double wn_clipping = 0.0;
    for(int i = 0; i < shape.size(); ++i)
    {
      //desmos_print(shape[i]);
      //desmos_print(query);
      wn_boxes += primal::winding_number(query, shape[i], nevals_boxes, tol, tol);
      wn_clipping +=
        primal::winding_number_clipping(query, shape[i], nevals_clipping, tol, tol);
      wn_bisection += primal::winding_number_bisection(query,
                                                       shape[i],
                                                       nevals_bisection,
                                                       tol,
                                                       tol);
    }

    if(!axom::utilities::isNearlyEqual(wn_bisection, wn_boxes, 1e-5))
    {
      std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                << std::endl;
      std::cout << "bisection machine broke" << std::endl;
    }

    if(!axom::utilities::isNearlyEqual(wn_boxes, wn_clipping, 1e-5))
    {
      std::cout << std::setprecision(16) << n << ": " << x << ", " << y
                << std::endl;
      std::cout << "clipping machine broke" << std::endl;
    }

    random_boxes_nevals_out << x << "," << y << "," << nevals_boxes << std::endl;
    random_clipping_nevals_out << x << "," << y << "," << nevals_clipping
                               << std::endl;
    random_bisection_nevals_out << x << "," << y << "," << nevals_bisection
                                << std::endl;
  }
}

TEST(primal_2d_paper_figure_data, orientation_bad)
{
  return;
  std::cout << "Running: \"orientation_bad\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = 
  "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\orientation_bad\\";

  int npts_x = 1000;
  int npts_y = 1000;

  std::string name = "butterfly";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  //BoundingBox bbox(shape[100].boundingBox());
  //double max_dim = axom::utilities::max(bbox.getMax()[0] - bbox.getMin()[0],
  //                                      bbox.getMax()[1] - bbox.getMin()[1]);
  //bbox.addPoint(primal::Point<double, 2> {bbox.getMin()[0] + max_dim,
  //                                        bbox.getMin()[1] + max_dim});
  //bbox.scale(2.5);

  BoundingBox bbox = curves_bbox(shape, 1.1, false);

  // Reverse one curve
  shape[100][1] = shape[100][0];  //.reverseOrientation();
  shape[100][2] = shape[100][0];  //.reverseOrientation();
  shape[100][3] = shape[100][0];  //.reverseOrientation();

  // Print them all to the file
  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  simple_grid_test(shape, bbox, npts_x, npts_y, wn_out);
}

TEST(primal_2d_paper_figure_data, quadrature_bad)
{
  return;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\quadrature_bad\\";


  int npts_x = 600;
  int npts_y = 600;

  std::string name = "curve";
  CurveArray curve;
  primal::convert_from_svg(data_dir + name + ".svg", curve);
  std::ofstream curve_out(data_dir + name + ".txt");
  for(int i = 0; i < curve.size(); ++i)
  {
    curve_out << curve[i] << std::endl;
  }
  BoundingBox curve_bbox = curves_bbox(curve, 1.5, true);
  std::ofstream curve_15_wn_out(data_dir + name + "_wn_15.csv");
  std::ofstream curve_30_wn_out(data_dir + name + "_wn_30.csv");
  std::ofstream curve_50_wn_out(data_dir + name + "_wn_50.csv");
  std::ofstream curve_true_wn_out(data_dir + name + "_wn_true.csv");

  for(double x = curve_bbox.getMin()[0]; x <= curve_bbox.getMax()[0];
      x += (curve_bbox.getMax()[0] - curve_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = curve_bbox.getMin()[1]; y <= curve_bbox.getMax()[1];
        y += (curve_bbox.getMax()[1] - curve_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn_15 = 0.0;
      double wn_30 = 0.0;
      double wn_50 = 0.0;
      double wn_true = 0.0;

      for(int i = 0; i < curve.size(); ++i)
      {
        wn_15 += primal::winding_number_quad(query, curve[i], 15);
        wn_30 += primal::winding_number_quad(query, curve[i], 30);
        wn_50 += primal::winding_number_quad(query, curve[i], 50);
        wn_true += primal::winding_number(query, curve[i], nevals);
      }

      curve_15_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_15
                      << std::endl;
      curve_30_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_30
                      << std::endl;
      curve_50_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_50
                      << std::endl;
      curve_true_wn_out << std::setprecision(16) << x << ',' << y << ","
                        << wn_true << std::endl;
    }
  }

  std::string spade_name = "twirl";
  CurveArray spade;
  primal::convert_from_svg(data_dir + spade_name + ".svg", spade);
  std::ofstream spade_out(data_dir + spade_name + ".txt");
  for(int i = 0; i < spade.size(); ++i)
  {
    spade_out << spade[i] << std::endl;
  }
  //BoundingBox spade_bbox = curves_bbox(spade, 1.0);
  primal::Point<double, 2> spade_center({83.73, 75});
  double spade_radius = 12.0;
  BoundingBox spade_bbox;
  spade_bbox.addPoint(primal::Point<double, 2> {spade_center[0] + spade_radius,
                                                spade_center[1] + spade_radius});
  spade_bbox.addPoint(primal::Point<double, 2> {spade_center[0] - spade_radius,
                                                spade_center[1] - spade_radius});

  std::ofstream spade_15_wn_out(data_dir + spade_name + "_wn_15.csv");
  std::ofstream spade_30_wn_out(data_dir + spade_name + "_wn_30.csv");
  std::ofstream spade_50_wn_out(data_dir + spade_name + "_wn_50.csv");
  std::ofstream spade_true_wn_out(data_dir + spade_name + "_wn_true.csv");

  for(double x = spade_bbox.getMin()[0]; x <= spade_bbox.getMax()[0];
      x += (spade_bbox.getMax()[0] - spade_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = spade_bbox.getMin()[1]; y <= spade_bbox.getMax()[1];
        y += (spade_bbox.getMax()[1] - spade_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn_15 = 0.0;
      double wn_30 = 0.0;
      double wn_50 = 0.0;
      double wn_true = 0.0;

      for(int i = 0; i < spade.size(); ++i)
      {
        wn_15 += primal::winding_number_quad(query, spade[i], 15);
        wn_30 += primal::winding_number_quad(query, spade[i], 30);
        wn_50 += primal::winding_number_quad(query, spade[i], 50);
        wn_true += primal::winding_number(query, spade[i], nevals);
      }

      spade_15_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_15
                      << std::endl;
      spade_30_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_30
                      << std::endl;
      spade_50_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_50
                      << std::endl;
      spade_true_wn_out << std::setprecision(16) << x << ',' << y << ","
                        << wn_true << std::endl;
    }
  }

  std::string zoomed_name = "zoomed";

  primal::Point<double, 2> zoomed_center = spade[2][3];
  double zoomed_radius = 1.0;
  BoundingBox zoomed_bbox;
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] + zoomed_radius,
                              zoomed_center[1] + zoomed_radius / 2.0});
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] - zoomed_radius,
                              zoomed_center[1] - zoomed_radius / 2.0});

  std::ofstream zoomed_15_wn_out(data_dir + zoomed_name + "_wn_15.csv");
  std::ofstream zoomed_30_wn_out(data_dir + zoomed_name + "_wn_30.csv");
  std::ofstream zoomed_50_wn_out(data_dir + zoomed_name + "_wn_50.csv");
  std::ofstream zoomed_true_wn_out(data_dir + zoomed_name + "_wn_true.csv");

  npts_x = 600;
  npts_y = 300;

  for(double x = zoomed_bbox.getMin()[0]; x <= zoomed_bbox.getMax()[0];
      x += (zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << std::endl;

    for(double y = zoomed_bbox.getMin()[1]; y <= zoomed_bbox.getMax()[1];
        y += (zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn_15 = 0.0;
      double wn_30 = 0.0;
      double wn_50 = 0.0;
      double wn_true = 0.0;

      for(int i = 0; i < spade.size(); ++i)
      {
        wn_15 += primal::winding_number_quad(query, spade[i], 15);
        wn_30 += primal::winding_number_quad(query, spade[i], 30);
        wn_50 += primal::winding_number_quad(query, spade[i], 50);
        wn_true += primal::winding_number(query, spade[i], nevals);
      }

      zoomed_15_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_15
                       << std::endl;
      zoomed_30_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_30
                       << std::endl;
      zoomed_50_wn_out << std::setprecision(16) << x << ',' << y << "," << wn_50
                       << std::endl;
      zoomed_true_wn_out << std::setprecision(16) << x << ',' << y << ","
                         << wn_true << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, graphical_abstract)
{
  return;
  std::cout << "Running: \"graphical_abstract\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\graphical_abstract\\";
//    "E:\\Code\\winding_number\\figures\\\\";

  int npts_x = 1000;
  int npts_y = 1000;

  std::string name = "butterfly";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  primal::convert_from_svg(data_dir + name + ".svg", shape);
  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
  }

  // Replace curve 100 with a point

  CurveArray exploded_shape;
  std::ofstream exploded_shape_out(data_dir + name + "_exploded.txt");
  std::ofstream exploded_shape_cn_out(data_dir + name + "_exploded_cn.csv");
  std::ofstream exploded_shape_wn_out(data_dir + name + "_exploded_wn.csv");
  primal::convert_from_svg(data_dir + name + "_exploded.svg", exploded_shape);

  for(int i = 0; i < exploded_shape.size(); ++i)
  {
    exploded_shape_out << exploded_shape[i] << std::endl;
  }

  BoundingBox shape_bbox = curves_bbox(shape, 1.1, false);

  // Get winding number and crossing number for the exploded shape
  BoundingBox exploded_shape_bbox = curves_bbox(exploded_shape, 1.1, false);

  for(double x = exploded_shape_bbox.getMin()[0];
      x <= exploded_shape_bbox.getMax()[0];
      x += (exploded_shape_bbox.getMax()[0] - exploded_shape_bbox.getMin()[0]) /
        npts_x)
  {
    std::cout << x << "/" << exploded_shape_bbox.getMax()[0] << std::endl;

    for(double y = exploded_shape_bbox.getMin()[1];
        y <= exploded_shape_bbox.getMax()[1];
        y += (exploded_shape_bbox.getMax()[1] - exploded_shape_bbox.getMin()[1]) /
          npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int crossing_number =
        primal::crossing_number(query, exploded_shape, exploded_shape_bbox);
      double winding_number = primal::winding_number(query, exploded_shape);

      exploded_shape_cn_out << std::setprecision(16) << x << ',' << y << ","
                            << crossing_number << std::endl;

      exploded_shape_wn_out << std::setprecision(16) << x << ',' << y << ","
                            << winding_number << std::endl;
    }
  }

  // Get a linearized version of the shape
  std::ofstream linearized_shape_out(data_dir + name + "_linearized.txt");

  // Construct the linear approximation with the given level of refinement
  CurveArray linearized_shape;
  for(int i = 0; i < exploded_shape.size(); ++i)
  {
    int n = 5;
    for(double j = 0; j < n; ++j)
    {
      primal::BezierCurve<double, 2> line(1);
      line[0] = exploded_shape[i].evaluate(j / n);
      line[1] = exploded_shape[i].evaluate((j + 1) / n);

      linearized_shape.push_back(line);
    }
  }

  std::ofstream linearized_shape_wn_out(data_dir + name + "_linearized_wn.csv");

  for(int i = 0; i < linearized_shape.size(); ++i)
  {
    linearized_shape_out << linearized_shape[i] << std::endl;
  }

  // Zoom in around curve 106
  std::string zoomed_name = "zoomed";

  primal::Point<double, 2> zoomed_center = exploded_shape[106].evaluate(0.85);
  double zoomed_radius = 400.0;
  BoundingBox zoomed_bbox;
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] + zoomed_radius,
                              zoomed_center[1] + zoomed_radius});
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] - zoomed_radius,
                              zoomed_center[1] - zoomed_radius});

  std::ofstream zoomed_true_out(data_dir + zoomed_name + "_true.csv");
  std::ofstream zoomed_linearized_out(data_dir + zoomed_name +
                                      "_linearized.csv");
  std::ofstream zoomed_quadrature_out(data_dir + zoomed_name +
                                      "_quadrature.csv");

  for(double x = zoomed_bbox.getMin()[0]; x <= zoomed_bbox.getMax()[0];
      x += (zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << "/" << exploded_shape_bbox.getMax()[0] << std::endl;

    for(double y = zoomed_bbox.getMin()[1]; y <= zoomed_bbox.getMax()[1];
        y += (zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double true_wn = 0.0;
      double linearized = 0.0;
      double quadrature = 0.0;

      for(int i = 0; i < exploded_shape.size(); ++i)
      {
        true_wn += primal::winding_number(query, exploded_shape[i], nevals);
        quadrature += primal::winding_number_quad(query, exploded_shape[i], 30);
      }

      for(int i = 0; i < linearized_shape.size(); ++i)
      {
        linearized += primal::winding_number(query, linearized_shape[i], nevals);
      }

      zoomed_true_out << std::setprecision(16) << x << ',' << y << ","
                      << true_wn << std::endl;
      zoomed_linearized_out << std::setprecision(16) << x << ',' << y << ","
                            << linearized << std::endl;
      zoomed_quadrature_out << std::setprecision(16) << x << ',' << y << ","
                            << quadrature << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, graphical_abstract_trumpet)
{
  return;

  std::cout << "Running: \"graphical_abstract_trumpet\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\graphical_abstract\\";
  //"E:\\Code\\winding_number\\figures\\graphical_abstract\\";

  int npts_x = 1000;
  int npts_y = 1000;

  CurveArray closed_trumpet;
  std::ofstream closed_trumpet_out(data_dir + "closed_trumpet.txt");
  primal::convert_from_svg(data_dir + "closed_trumpet.svg", closed_trumpet);
  for(int i = 0; i < closed_trumpet.size(); ++i)
  {
    closed_trumpet_out << closed_trumpet[i] << std::endl;
  }
  BoundingBox closed_bbox = curves_bbox(closed_trumpet, 1.1, false);

  CurveArray open_trumpet;
  std::ofstream open_trumpet_out(data_dir + "open_trumpet.txt");
  primal::convert_from_svg(data_dir + "open_trumpet.svg", open_trumpet);
  for(int i = 0; i < open_trumpet.size(); ++i)
  {
    open_trumpet_out << open_trumpet[i] << std::endl;
  }
  BoundingBox open_bbox = curves_bbox(open_trumpet, 1.1, false);

  std::ofstream closed_trumpet_wn_out(data_dir + "closed_trumpet_wn.csv");
  std::ofstream closed_trumpet_cn_out(data_dir + "closed_trumpet_cn.csv");
  std::ofstream open_trumpet_wn_out(data_dir + "open_trumpet_wn.csv");
  std::ofstream open_trumpet_cn_out(data_dir + "open_trumpet_cn.csv");

  // Compute winding numbers and crossing numbers for both trumpets
  for(double x = open_bbox.getMin()[0]; x <= open_bbox.getMax()[0];
      x += (open_bbox.getMax()[0] - open_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << "/" << open_bbox.getMax()[0] << std::endl;

    for(double y = open_bbox.getMin()[1]; y <= open_bbox.getMax()[1];
        y += (open_bbox.getMax()[1] - open_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int crossing_number =
        primal::crossing_number(query, open_trumpet, open_bbox);
      double winding_number = primal::winding_number(query, open_trumpet);

      open_trumpet_wn_out << std::setprecision(16) << x << ',' << y << ","
                          << winding_number << std::endl;
      open_trumpet_cn_out << std::setprecision(16) << x << ',' << y << ","
                          << crossing_number << std::endl;
    }
  }

  for(double x = closed_bbox.getMin()[0]; x <= closed_bbox.getMax()[0];
      x += (closed_bbox.getMax()[0] - closed_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << "/" << closed_bbox.getMax()[0] << std::endl;

    for(double y = closed_bbox.getMin()[1]; y <= closed_bbox.getMax()[1];
        y += (closed_bbox.getMax()[1] - closed_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int crossing_number =
        primal::crossing_number(query, closed_trumpet, closed_bbox);
      double winding_number = primal::winding_number(query, closed_trumpet);

      closed_trumpet_wn_out << std::setprecision(16) << x << ',' << y << ","
                            << winding_number << std::endl;
      closed_trumpet_cn_out << std::setprecision(16) << x << ',' << y << ","
                            << crossing_number << std::endl;
    }
  }

  // Get a linearized version of the shape
  std::ofstream linearized_shape_out(data_dir + "zoomed_linearized.txt");

  // Construct the linear approximation with the given level of refinement
  CurveArray linearized_shape;
  for(int i = 0; i < open_trumpet.size(); ++i)
  {
    int n = 5;
    for(double j = 0; j < n; ++j)
    {
      primal::BezierCurve<double, 2> line(1);
      line[0] = open_trumpet[i].evaluate(j / n);
      line[1] = open_trumpet[i].evaluate((j + 1) / n);

      linearized_shape.push_back(line);
    }
  }

  for(int i = 0; i < linearized_shape.size(); ++i)
  {
    linearized_shape_out << linearized_shape[i] << std::endl;
  }

  primal::Point<double, 2> zoomed_center =
    open_trumpet[open_trumpet.size() - 1].evaluate(0.21);
  double zoomed_radius = 22.0;
  BoundingBox zoomed_bbox;
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] + zoomed_radius,
                              zoomed_center[1] + zoomed_radius});
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_center[0] - zoomed_radius,
                              zoomed_center[1] - zoomed_radius});

  std::ofstream zoomed_true_out(data_dir + "zoomed_true.csv");
  std::ofstream zoomed_linearized_out(data_dir + "zoomed_linearized.csv");
  std::ofstream zoomed_quadrature_out(data_dir + "zoomed_quadrature.csv");

  npts_x = 1000;
  npts_y = 1000;

  for(double x = zoomed_bbox.getMin()[0]; x <= zoomed_bbox.getMax()[0];
      x += (zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0]) / npts_x)
  {
    std::cout << x << "/" << zoomed_bbox.getMax()[0] << std::endl;

    for(double y = zoomed_bbox.getMin()[1]; y <= zoomed_bbox.getMax()[1];
        y += (zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double true_wn = 0.0;
      double linearized = 0.0;
      double quadrature = 0.0;

      for(int i = 0; i < open_trumpet.size(); ++i)
      {
        true_wn += primal::winding_number(query, open_trumpet[i], nevals);
        quadrature += primal::winding_number_quad(query, open_trumpet[i], 30);
      }

      for(int i = 0; i < linearized_shape.size(); ++i)
      {
        linearized += primal::winding_number(query, linearized_shape[i], nevals);
      }

      zoomed_true_out << std::setprecision(16) << x << ',' << y << ","
                      << true_wn << std::endl;
      zoomed_linearized_out << std::setprecision(16) << x << ',' << y << ","
                            << linearized << std::endl;
      zoomed_quadrature_out << std::setprecision(16) << x << ',' << y << ","
                            << quadrature << std::endl;
    }
  }
}

TEST(primal_2d_paper_figure_data, animation_still)
{
  //return;
  std::cout << "Running: \"animation_still\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
  "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\animation\\";

  int npts_x = 200 * 11;
  int npts_y = 200 * 8.5;

  std::string name = "all_logos";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
  }

  BoundingBox bbox = curves_bbox(shape, 1.2, false);
  // Scale the box so that it has a ratio of 11x8.5
  primal::Point<double, 2> min_pt = bbox.getMin();
  primal::Point<double, 2> max_pt = bbox.getMax();
  double ratio = 11.0 / 8.5;
  //double ratio = 8.5 / 11.0;
  double width = max_pt[0] - min_pt[0];
  double height = max_pt[1] - min_pt[1];
  if(width / height > ratio)
  {
    double new_height = width / ratio;
    min_pt[1] -= (new_height - height) / 2.0;
    max_pt[1] += (new_height - height) / 2.0;
  }
  else
  {
    double new_width = height * ratio;
    min_pt[0] -= (new_width - width) / 2.0;
    max_pt[0] += (new_width - width) / 2.0;
  }
  simple_grid_test(shape, bbox, npts_x, npts_y, wn_out);
}

TEST(primal_2d_paper_figure_data, background)
{
  return;
  std::cout << "Running: \"background\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\animation\\";

  int npts_x = 400 * 16;
  int npts_y = 400 * 9;

  std::string name = "background";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
  }

  BoundingBox bbox = curves_bbox(shape, 1.0, false);

  simple_grid_test(shape, bbox, npts_x, npts_y, wn_out);
}

TEST(primal_2d_paper_figure_data, animation)
{
  return;
  std::cout << "Running: \"animation\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir = "E:\\Code\\winding_number\\figures\\animation\\";

  int permutation[] = {
    111, 32,  88, 120, 101, 66,  68, 7,   54,  113, 65,  117, 129, 67,  110,
    125, 92,  85, 95,  8,   23,  83, 122, 26,  64,  3,   104, 51,  22,  30,
    63,  123, 84, 0,   91,  99,  62, 40,  61,  42,  94,  102, 33,  20,  100,
    126, 49,  53, 69,  107, 56,  78, 58,  71,  127, 108, 13,  130, 41,  43,
    19,  5,   31, 15,  86,  80,  79, 128, 90,  11,  24,  9,   16,  121, 116,
    14,  89,  74, 4,   25,  28,  57, 97,  47,  10,  82,  87,  133, 46,  44,
    36,  55,  1,  81,  119, 114, 35, 21,  27,  77,  98,  124, 6,   38,  45,
    34,  106, 29, 76,  96,  52,  48, 39,  105, 17,  72,  73,  37,  75,  103,
    70,  109, 2,  18,  132, 131, 12, 59,  60,  93,  112, 50,  118, 115};

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "butterfly";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream shuffled_shape_out(data_dir + name + "shuffled_.txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  CurveArray shuffled_shape;
  primal::convert_from_svg(data_dir + name + "_simple_bonus_edger.svg",
                           shuffled_shape);

  BoundingBox bbox = curves_bbox(shape, 1.25, true);
  primal::Point<double, 2> min_pt = bbox.getMin();
  primal::Point<double, 2> max_pt = bbox.getMax();
  auto length = bbox.getMax()[0] - bbox.getMin()[0];
  min_pt[0] -= length * 7.0 / 9.0 * 0.9;
  max_pt[0] += length * 7.0 / 9.0 * 0.1;
  bbox.addPoint(min_pt);
  bbox.addPoint(max_pt);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape_out << shape[i] << std::endl;
    shuffled_shape_out << shuffled_shape[i] << std::endl;
  }

  // Curves one at a time
  if(false)
  {
    CurveArray shape_increment;

    // Print them all to the file
    for(int i = 0; i < shape.size(); ++i)
    {
      shape_increment.push_back(shape[permutation[i]]);
      std::ofstream increment_shape_out(data_dir + "animation1_frames\\" +
                                        name + "_" + std::to_string(i) + ".txt");
      std::ofstream increment_wn_out(data_dir + "animation1_frames\\" + name +
                                     "_" + std::to_string(i) + "_wn.csv");

      for(int j = 0; j < shape_increment.size(); ++j)
      {
        increment_shape_out << shape_increment[j] << std::endl;
      }

      simple_grid_test(shape_increment, bbox, npts_x, npts_y, increment_wn_out);
    }
  }

  // Explode n Back
  if(false)
  {
    auto interp_curve = [](const primal::BezierCurve<double, 2>& curve1,
                           const primal::BezierCurve<double, 2>& curve2,
                           double t) -> primal::BezierCurve<double, 2> {
      primal::BezierCurve<double, 2> result(curve1.getOrder());
      for(int i = 0; i <= curve1.getOrder(); ++i)
      {
        result[i][0] = curve1[i][0] * (1.0 - t) + curve2[i][0] * t;
        result[i][1] = curve1[i][1] * (1.0 - t) + curve2[i][1] * t;
      }

      return result;
    };

    double a = 10;
    auto speed_func = [&a](double t) -> double {
      double g0 = 1.0 / (1 + std::exp(0.5 * a));
      return (g0 - 1.0 / (1 + std::exp(-a * (t - 0.5)))) / (2 * g0 - 1);
    };

    for(int frame = 0; frame <= 100; ++frame)
    {
      double t;
      if(frame < 50)
      {
        t = speed_func(frame / 50.0);
      }
      else
      {
        t = speed_func((100 - frame) / 50.0);
      }

      CurveArray midshift_shape;
      std::ofstream midshift_shape_out(data_dir + "animation3_frames\\frame_" +
                                       std::to_string(frame) + "_shape.txt");
      std::ofstream midshift_wn_out(data_dir + "animation3_frames\\frame_" +
                                    std::to_string(frame) + "_wn.csv");

      for(int i = 0; i < shape.size(); ++i)
      {
        primal::BezierCurve<double, 2> interp =
          interp_curve(shape[i], shuffled_shape[i], t);
        midshift_shape_out << interp << std::endl;
        midshift_shape.push_back(interp);
      }

      simple_grid_test(midshift_shape, bbox, npts_x, npts_y, midshift_wn_out);
    }
  }

  // Slide a curve across the screen
  if(true)
  {
    auto interp_curve = [](const primal ::BezierCurve<double, 2>& curve1,
                           const primal::BezierCurve<double, 2>& curve2,
                           double t) -> primal::BezierCurve<double, 2> {
      primal::BezierCurve<double, 2> result(curve1.getOrder());
      for(int i = 0; i <= curve1.getOrder(); ++i)
      {
        result[i][0] = curve1[i][0] * (1.0 - t) + curve2[i][0] * t;
        result[i][1] = curve1[i][1] * (1.0 - t) + curve2[i][1] * t;
      }

      return result;
    };

    double a = 5;
    auto speed1 = [&a](double t) -> double {
      double g0 = 1.0 / (1 + std::exp(0.5 * a));
      return (g0 - 1.0 / (1 + std::exp(-a * (t - 0.5)))) / (2 * g0 - 1);
    };

    auto speed2 = [](double t) -> double {
      if(t < 0.5)
        return (std::pow(20, t) - 1) / (std::sqrt(20) - 1);
      else
        return (std::pow(20, 1 - t) - 1) / (std::sqrt(20) - 1);
    };

    shuffled_shape[103].reverseOrientation();
    primal::BezierCurve<double, 2> midtwist_curve(3);
    midtwist_curve[0] = primal::Point<double, 2> {1417.0, 67.0};
    midtwist_curve[1] = primal::Point<double, 2> {1920.88, -95.604};
    midtwist_curve[2] = primal::Point<double, 2> {1920.88, -95.604};
    midtwist_curve[3] = primal::Point<double, 2> {2473.0, -240.0};

    double nframes = 150;
    for(int frame = 0; frame <= nframes; ++frame)
    {
      double t = speed1(frame / nframes);

      CurveArray midshift_shape;
      std::ofstream midshift_shape_out(data_dir + "animation3_frames\\frame_" +
                                       std::to_string(frame) + "_shape.txt");
      std::ofstream midshift_wn_out(data_dir + "animation3_frames\\frame_" +
                                    std::to_string(frame) + "_wn.csv");

      primal::BezierCurve<double, 2> motion_curve1 =
        shuffled_shape[shuffled_shape.size() - 2];
      primal::BezierCurve<double, 2> motion_curve2 =
        shuffled_shape[shuffled_shape.size() - 1];

      primal::Point<double, 2> new_point, reference_point = shape[0][0];

      if(frame < (nframes / 2))
        new_point = motion_curve1.evaluate(speed2(frame / nframes));
      else
        new_point = motion_curve2.evaluate(1.0 - speed2(frame / nframes));

      primal::Point<double, 2> the_shift {new_point[0] - reference_point[0],
                                          new_point[1] - reference_point[1]};

      primal::BezierCurve<double, 2> new_curve_1(3), new_curve_2(3);

      for(int i = 0; i <= 3; ++i)
      {
        new_curve_1[i][0] = shape[0][i][0] + the_shift[0];
        new_curve_1[i][1] = shape[0][i][1] + the_shift[1];

        new_curve_2[i][0] = shape[1][i][0] + the_shift[0];
        new_curve_2[i][1] = shape[1][i][1] + the_shift[1];
      }

      midshift_shape_out << new_curve_1 << std::endl;
      midshift_shape_out << new_curve_2 << std::endl;

      midshift_shape.push_back(new_curve_1);
      midshift_shape.push_back(new_curve_2);

      //shuffled_shape[103][0] = shape[103][3];
      //shuffled_shape[103][1][0] += 2.5;
      //shuffled_shape[103][2][0] += 2.5;
      //shuffled_shape[103][3] = shape[103][0];

      for(int i = 2; i < shape.size(); ++i)
      {
        primal::BezierCurve<double, 2> interp(3);

        if(i == 103)
          interp =
            interp_curve(interp_curve(shape[i], midtwist_curve, t),
                         interp_curve(midtwist_curve, shuffled_shape[i], t),
                         t);
        else
          interp = interp_curve(shape[i], shuffled_shape[i], t);
        midshift_shape_out << interp << std::endl;
        midshift_shape.push_back(interp);
      }

      simple_grid_test(midshift_shape, bbox, npts_x, npts_y, midshift_wn_out);

      if(frame == nframes) shape = midshift_shape;
    }
  }

  // zoom in
  if(false)
  {
    int idx = 109;
    auto zoom_pt = shape[idx][shape[idx].getOrder()];
    zoom_pt[1] -= 0.5;
    double zoom_r = 1;
    BoundingBox zoomed_bbox;
    zoomed_bbox.addPoint(
      primal::Point<double, 2> {zoom_pt[0] + zoom_r, zoom_pt[1] + zoom_r});
    zoomed_bbox.addPoint(
      primal::Point<double, 2> {zoom_pt[0] - zoom_r, zoom_pt[1] - zoom_r});

    primal::Point<double, 2> zoomed_min_pt = zoomed_bbox.getMin();
    primal::Point<double, 2> zoomed_max_pt = zoomed_bbox.getMax();

    zoomed_max_pt[0] += 2 * zoom_r * 7.0 / 9.0 * 0.1;
    zoomed_min_pt[0] -= 2 * zoom_r * 7.0 / 9.0 * 0.9;
    zoomed_bbox.addPoint(zoomed_min_pt);
    zoomed_bbox.addPoint(zoomed_max_pt);

    double a = 15;
    auto speed_func = [&a](double t) -> double {
      double g0 = 1.0 / (1 + std::exp(0.5 * a));
      return (g0 - 1.0 / (1 + std::exp(-a * (t - 0.5)))) / (2 * g0 - 1);
    };

    auto interp_box = [](const primal::BoundingBox<double, 2>& box1,
                         const primal::BoundingBox<double, 2>& box2,
                         double t) -> primal::BoundingBox<double, 2> {
      primal::BoundingBox<double, 2> result;
      result.addPoint(primal::Point<double, 2> {
        box1.getMin()[0] * (1.0 - t) + box2.getMin()[0] * t,
        box1.getMin()[1] * (1.0 - t) + box2.getMin()[1] * t});

      result.addPoint(primal::Point<double, 2> {
        box1.getMax()[0] * (1.0 - t) + box2.getMax()[0] * t,
        box1.getMax()[1] * (1.0 - t) + box2.getMax()[1] * t});

      return result;
    };

    double nframes = 150;
    for(int frame = 0; frame <= nframes; ++frame)
    {
      double t = speed_func(frame / nframes);

      std::ofstream zoomed_wn_out(data_dir + "animation2_frames\\frame_" +
                                  std::to_string(frame) + "_wn.csv");

      auto semizoomed_bbox = interp_box(bbox, zoomed_bbox, t);
      simple_grid_test(shape, semizoomed_bbox, npts_x, npts_y, zoomed_wn_out);
    }
  }
}

TEST(primal_2d_paper_figure_data, fish_zoom)
{
  return;
  std::cout << "Running: \"fish zoom\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\animation\\";

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "fish";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  // Jiggle the curves, and print them all to the file
  srand(100);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();

    shape[i][0][0] += (rand() % 30) - 15;
    shape[i][0][1] += (rand() % 30) - 15;

    shape[i][shape[i].getOrder()][0] += (rand() % 30) - 15;
    shape[i][shape[i].getOrder()][0] += (rand() % 30) - 15;

    shape_out << shape[i] << std::endl;
  }

  BoundingBox bbox = curves_bbox(shape, 1.1, true);

  // zoom in
  BoundingBox zoomed_bbox(shape[240].boundingBox());
  double max_dim =
    axom::utilities::max(zoomed_bbox.getMax()[0] - zoomed_bbox.getMin()[0],
                         zoomed_bbox.getMax()[1] - zoomed_bbox.getMin()[1]);
  zoomed_bbox.addPoint(
    primal::Point<double, 2> {zoomed_bbox.getMin()[0] + max_dim,
                              zoomed_bbox.getMin()[1] + max_dim});
  zoomed_bbox.scale(12.0);

  double a = 15;
  auto speed_func = [&a](double t) -> double {
    double g0 = 1.0 / (1 + std::exp(0.5 * a));
    return (g0 - 1.0 / (1 + std::exp(-a * (t - 0.5)))) / (2 * g0 - 1);
  };

  auto interp_box = [](const primal::BoundingBox<double, 2>& box1,
                       const primal::BoundingBox<double, 2>& box2,
                       double t) -> primal::BoundingBox<double, 2> {
    primal::BoundingBox<double, 2> result;
    result.addPoint(primal::Point<double, 2> {
      box1.getMin()[0] * (1.0 - t) + box2.getMin()[0] * t,
      box1.getMin()[1] * (1.0 - t) + box2.getMin()[1] * t});

    result.addPoint(primal::Point<double, 2> {
      box1.getMax()[0] * (1.0 - t) + box2.getMax()[0] * t,
      box1.getMax()[1] * (1.0 - t) + box2.getMax()[1] * t});

    return result;
  };

  double nframes = 150;
  for(int frame = 0; frame <= nframes; ++frame)
  {
    double t = speed_func(frame / nframes);

    std::ofstream zoomed_wn_out(data_dir + "fish_animation\\frame_" +
                                std::to_string(frame) + "_wn.csv");

    auto semizoomed_bbox = interp_box(bbox, zoomed_bbox, t);
    simple_grid_test(shape, semizoomed_bbox, npts_x, npts_y, zoomed_wn_out);
  }
}

TEST(primal_2d_paper_figure_data, axom_trace)
{
  return;
  std::cout << "Running: \"axom trace\"" << std::endl;

  using CurveArray = axom::Array<primal::BezierCurve<double, 2>>;
  using BoundingBox = primal::BoundingBox<double, 2>;

  std::string data_dir =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\figures\\axom_animation\\";

  int npts_x = 500;
  int npts_y = 500;

  std::string name = "siggraph_logo";
  CurveArray shape;
  std::ofstream shape_out(data_dir + name + ".txt");
  std::ofstream wn_out(data_dir + name + "_wn.csv");
  primal::convert_from_svg(data_dir + name + ".svg", shape);

  for(int i = 0; i < shape.size(); ++i)
  {
    shape[i].reverseOrientation();
    shape_out << shape[i] << std::endl;
  }

  BoundingBox bbox = curves_bbox(shape, 1.2, false);

  double nframes = 200;
  for(int frame = 0; frame <= nframes; ++frame)
  {
    double t = frame / nframes;

    std::ofstream traced_wn_out(data_dir + "frames_2\\frame_" +
                                    std::to_string(frame) + "_wn.csv");
    std::ofstream traced_shape_out(data_dir + +"frames_2\\frame_" +
                                   std::to_string(frame) + "_curves.csv");

    CurveArray traced_shape;
    
    primal::BezierCurve<double, 2> real_curve, dummy;
    for(int i = 0; i < shape.size(); ++i)
    {
      shape[i].split(t, real_curve, dummy);
      traced_shape_out << real_curve << std::endl;

      traced_shape.push_back(real_curve);
    }

    simple_grid_test(traced_shape, bbox, npts_x, npts_y, traced_wn_out);
  }
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
