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
void split_parabola_test();
void square_test();
void big_square_test();
void segment_test();
std::function<Vector2D(Point2D)> get_winding_func(Point2D p);

int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  axom::slic::SimpleLogger logger;
  //square_test();
  parabola_test();
  //split_parabola_test();
  //segment_test();

  return 0;
}

void segment_test()
{
  Point2D segnodes[] = {Point2D {1.0, -1.0}, Point2D {1.0, 1.0}};
  Bezier segment(segnodes, 1);

  int npts = 4;
  for(int i = 0; i < 5; i++)
  {
    Point2D qpoint({1 - std::pow(10, -i), 0.0});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(segment, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
  
  Point2D badpoint({1.0, 0.0});

  auto bad_winding_func = get_winding_func(badpoint);
  double bad_winding_num =
    evaluate_line_integral(segment, bad_winding_func, npts);

  std::cout << badpoint << ": " << bad_winding_num << std::endl;

  for(int i = 4; i > 0; i--)
  {
    Point2D qpoint({1 + std::pow(10, -i), 0.0});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(segment, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
}



void parabola_test()
{
  SLIC_INFO("Parabola Test");
  Point2D paranodes1[] = {Point2D {1.0, 0.0},
                          Point2D {0.0, 2.0},
                          Point2D {-1.0, 0.0}};
  Bezier para1(paranodes1, 2);

  Point2D paranodes2[] = {Point2D {-1.0, 0.0},
                          Point2D {0.0, -2.0},
                          Point2D {1.0, 0.0}};
  Bezier para2(paranodes2, 2);

  Bezier pedges[2] = {para1, para2};
  CPolygon parabola_polygon(pedges, 2);
  int npts = 5;

  Bezier closer = para1.linear_closure();
  std::cout << para1 << std::endl;
  std::cout << closer << std::endl;

  for(int i = 0; i < 5; i++)
  {
    Point2D qpoint({0.0, 1 - std::pow(10, -i)});
    
    auto winding_func = get_winding_func(qpoint);
    double winding_num = evaluate_line_integral(parabola_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }

  for(int i = 4; i > 0; i--)
  {
    Point2D qpoint({0.0, 1 + std::pow(10, -i)});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(parabola_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
}

void split_parabola_test()
{
  SLIC_INFO("Split Parabola Test");
  Point2D paranodes1[] = {Point2D {1.0, 0.0},
                          Point2D {0.5, 1.0},
                          Point2D {0.0, 1.0}};
  Bezier para1(paranodes1, 2);

  Point2D paranodes15[] = {Point2D { 0.0, 1.0},
                           Point2D {-0.5, 1.0},
                           Point2D {-1.0, 0.0}};
  Bezier para15(paranodes15, 2);

  Point2D paranodes2[] = {Point2D {-1.0, 0.0},
                          Point2D {0.0, -2.0},
                          Point2D {1.0, 0.0}};
  Bezier para2(paranodes2, 2);

  Bezier pedges[3] = {para1, para15, para2};
  CPolygon parabola_polygon(pedges, 3);
  int npts = 30;
  
  for(int i = 0; i < 5; i++)
  {
    Point2D qpoint({0.0, 1 + std::pow(10, -i)});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(parabola_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }

  for(int i = 0; i < 5; i++)
  {
    Point2D qpoint({0.0, 1 - std::pow(10, -i)});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(parabola_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
}

void square_test()
{
  SLIC_INFO("Square Test");
  Point2D snodes1[] = {Point2D {-1.0, -1.0},
                       Point2D { 1.0, -1.0}};
  Bezier square1(snodes1, 1);

  Point2D snodes2[] = {Point2D {1.0, -1.0}, Point2D {1.0, 1.0}};
  Bezier square2(snodes2, 1);

  Point2D snodes3[] = {Point2D {1.0, 1.0}, Point2D {-1.0, 1.0}};
  Bezier square3(snodes3, 1);

  Point2D snodes4[] = {Point2D {-1.0, 1.0}, Point2D {-1.0, -1.0}};
  Bezier square4(snodes4, 1);

  Bezier pedges[4] = {square1, square2, square3, square4};
  CPolygon square_polygon(pedges, 4);
  int npts = 15;

  for(int i = 0; i < 5; i++)
  {
    double dist = std::pow(10, -i);
    Point2D qpoint({1 - dist, 0.0});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(square_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
}

void big_square_test()
{
  SLIC_INFO("Big Square Test");
  Point2D snodes1[] = {Point2D {-2.0, -2.0}, Point2D {2.0, -2.0}};
  Bezier square1(snodes1, 1);

  Point2D snodes2[] = {Point2D {2.0, -2.0}, Point2D {2.0, 2.0}};
  Bezier square2(snodes2, 1);

  Point2D snodes3[] = {Point2D {2.0, 2.0}, Point2D {-2.0, 2.0}};
  Bezier square3(snodes3, 1);

  Point2D snodes4[] = {Point2D {-2.0, 2.0}, Point2D {-2.0, -2.0}};
  Bezier square4(snodes4, 1);

  Bezier pedges[4] = {square1, square2, square3, square4};
  CPolygon square_polygon(pedges, 4);
  int npts = 15;

  for(int i = 0; i < 5; i++)
  {
    double dist = std::pow(10, -i);
    Point2D qpoint({2 - 2*dist, 0.0});

    auto winding_func = get_winding_func(qpoint);
    double winding_num =
      evaluate_line_integral(square_polygon, winding_func, npts);

    std::cout << qpoint << ": " << winding_num << std::endl;
  }
}

std::function<Vector2D(Point2D)> get_winding_func(Point2D p)
{
  return [p](Point2D x) -> Vector2D {
    double denom =
      2 * M_PI * ((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]));
    //std::cout << x << std::endl;
    //std::cout << denom << std::endl;
    return Vector2D({-(x[1] - p[1]) / denom, (x[0] - p[0]) / denom});
  };
}