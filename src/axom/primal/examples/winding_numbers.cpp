// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

// C++ headers
#include <cmath>
#include <iostream>
#include <fstream>

namespace primal = axom::primal;
using Point2D = primal::Point<double, 2>;
using Vector2D = primal::Vector<double, 2>;

using Segment = primal::Segment<double, 2>;
using Triangle = primal::Triangle<double, 2>;
using Bezier = primal::BezierCurve<double, 2>;
using CPolygon = primal::CurvedPolygon<double, 2>;

using vector_field = std::function<Vector2D(Point2D)>;


void winding_number_grid();

CPolygon get_self_intersecting_shape();
CPolygon get_large_shape();

int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  axom::slic::SimpleLogger logger;
  winding_number_grid();

  return 0;
}

CPolygon get_self_intersecting_shape() 
{ 
  // self intersecting cubic
  Point2D big_nodes[] = {Point2D {0.0, 0.0},
                         Point2D {0.0, 2.0},
                         Point2D {-3.0, -1.0},
                         Point2D {1.5, -1.0},
                         Point2D {1.5, 2.0},
                         Point2D {-3.0, 2.0},
                         Point2D {0.0, -1.0},
                         Point2D {0.0, 1.0}};
  Bezier super_intersecting(big_nodes, 7);

  Point2D small_nodes[] = {Point2D {0.0, 1.0},
                           Point2D {-2.0, 0.5},
                           Point2D {0.0, 0.0}};

  Bezier intersecting_closure(small_nodes, 2);
  Bezier super_intersecting_nodes[] = {super_intersecting, intersecting_closure};
  return CPolygon(super_intersecting_nodes, 2);
}

CPolygon get_large_shape() 
{
  // self intersecting cubic
  Point2D nodes1[] = {Point2D {0.0, 0.0},
                      Point2D {2.0, 1.0},
                      Point2D {-1.0, 1.0},
                      Point2D {1.0, 0.0}};
  Bezier self_intersecting_cubic(nodes1, 3);

  // kind-of circlish cubic
  Point2D nodes2[] = {Point2D {1.0, 0.0},
                      Point2D {2.0, 0.0},
                      Point2D {2.0, -1.0},
                      Point2D {1.0, -1.0}};
  Bezier semicircleish_cubic(nodes2, 3);

  // straight line segment
  Point2D nodes3[] = {Point2D {1.0, -1.0}, Point2D {-1.0, -2.0}};
  Bezier segment_linear(nodes3, 1);

  // parabola that intersects with the straight line
  Point2D nodes4[] = {Point2D {-1.0, -2.0},
                      Point2D {1.0, -3.0},
                      Point2D {-1.0, -1.0}};
  Bezier intersecting_parabola(nodes4, 2);

  // semi-circleish parabola that closes the curve
  Point2D nodes5[] = {Point2D {-1.0, -1.0},
                      Point2D {-1.0, 0.0},
                      Point2D {-0.1, 0.0}};
  Bezier semicircleish_parabola(nodes5, 2);

  Bezier pedges[] = {self_intersecting_cubic,
                     semicircleish_cubic,
                     segment_linear,
                     intersecting_parabola,
                     semicircleish_parabola};
  return CPolygon(pedges, 5);
}

void winding_number_grid()
{
  CPolygon cpoly = get_self_intersecting_shape();
  
  // Get big ol grid of query points
  Bezier::BoundingBoxType cpbb(cpoly.boundingBox().scale(1.1));
  const int num_pts = 400;
  double xpts[num_pts];
  double ypts[num_pts];

  axom::numerics::linspace(cpbb.getMin()[0], cpbb.getMax()[0], xpts, num_pts);
  axom::numerics::linspace(cpbb.getMin()[1], cpbb.getMax()[1], ypts, num_pts);
  //axom::numerics::linspace(-1.5, 1.0, xpts, num_pts);
  //axom::numerics::linspace(-0.5, 1.5, ypts, num_pts);
  //axom::numerics::linspace(-0.4, -0.2, xpts, num_pts);
  //axom::numerics::linspace(-1.8, -1.6, ypts, num_pts);
  //axom::numerics::linspace(0.0, 1.0, xpts, num_pts);
  //axom::numerics::linspace(0.0, 1.0, ypts, num_pts);

  // Get file storage syntax
  std::ofstream outfile(
    "C:/Users/spainhour1/Documents/bezier_plotting/datafile.csv");
  if(!outfile.good())
  {
    std::cout << "Could not write to the file" << std::endl;
    return;
  }

  // Store quadrature order
  int qnodes = 15;

  // Loop over each query point, store it and the computed winding number there
  for(double& x : xpts)
  {
    for(double& y : ypts)
    {
      Point2D qpoint({x, y});
      int max_depth = 0;
      double winding_num =
        winding_number(cpoly, qpoint, qnodes, max_depth);
      // clang-format off
      outfile << axom::fmt::format(
        "{0},{1},{2},{3}\n",
        qpoint[0],
        qpoint[1],
        winding_num, 
        max_depth
      );
      // clang-format on
    }
    std::cout << ".";
  }
  std::cout << std::endl;
}

