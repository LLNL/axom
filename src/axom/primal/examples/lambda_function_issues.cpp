// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

// C++ headers
#include <cmath>
#include <iostream>

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using Vector2D = primal::Vector<double, 2>;

using Bezier = primal::BezierCurve<double, 2>;
using CPolygon = primal::CurvedPolygon<double, 2>;

void parabola_test();
void square_test();
std::function<Vector2D(Point2D)> get_vector_field(Point2D p);
std::function<Vector2D(Point2D)> get_lame_winding_func(Point2D p);

  int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);
  
  axom::slic::SimpleLogger logger;
  parabola_test();

  return 0;
}

void parabola_test() 
{

  Point2D paranodes1[] = {Point2D { 1.0, 0.0},
                          Point2D { 0.0, 2.0},
                          Point2D {-1.0, 0.0}};
  Bezier para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {-1.0,  0.0},
                          Point2D { 0.0, -2.0},
                          Point2D { 1.0,  0.0}};
  Bezier para2(paranodes2, 2);

  Bezier pedges[2] = {para1, para2};
  CPolygon parabola_polygon(pedges, 2);

  Point2D p({0.0, 0.0});

  int npts = 15;

  auto func1 = get_vector_field(p);
  evaluate_line_integral(parabola_polygon, func1, npts);

  auto func2 = get_lame_winding_func(p);
  evaluate_line_integral(parabola_polygon, func2, npts);

  auto func3 = [p](Point2D x) -> Vector2D { return Vector2D({p[1], p[0]}); };
  evaluate_line_integral(parabola_polygon, func3, npts);
  
  auto func4 = [](Point2D x) -> Vector2D { return Vector2D({0.0, 0.0}); };
  evaluate_line_integral(parabola_polygon, func4, npts);
}

void square_test() { }

std::function<Vector2D(Point2D)> get_vector_field(Point2D p) 
{
  return [p](Point2D x) -> Vector2D { return Vector2D({p[1], p[0]}); };
}

std::function<Vector2D(Point2D)> get_lame_winding_func(Point2D p)
{
  return [](Point2D x) -> Vector2D { return Vector2D({0.0, 0.0}); };
}