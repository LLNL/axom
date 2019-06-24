// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/* /file bezier_test.cpp
 * /brief This file tests the BezierCurve.hpp and eval_bezier.hpp files
*/

// _prims_header_start
// Axom primitives
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
// _prims_header_end

// _eval_bez_start
#include "axom/primal/operators/eval_bezier.hpp"
// _eval_bez_end

// C++ headers
#include <iostream>
#include <fstream>

#include "fmt/fmt.hpp"

// _using_start
// "using" directives to simplify code
using namespace axom;
using namespace primal;

// almost all our examples are in 3D
constexpr int in3D = 3;

// primitives represented by doubles in 3D
typedef Point<double, in3D> PointType;
typedef BezierCurve<double, in3D> BezierCurveType;
// _using_end

BezierCurveType testBezier()
{
  //_ctrlpts_start
  BezierCurveType bCurve(1);
  std::cout << "----------------------Checking Bezier Functions-----------------------" << std::endl;
  std::cout << "Checking the order constructor:" << std::endl;
  std::cout << "Expected order: " << 1 << ". Order obtained from getOrder: " << bCurve.getOrder() << "." <<  std::endl;
  std::cout << "Adding points (.6, 1.2, 1.0) and (0.0 , 1.6, 1.8)" << std::endl;
  
  bCurve.addControlpoint( PointType::make_point( 0.6, 1.2, 1.0 ) );
  bCurve.addControlpoint( PointType::make_point( 0.0, 1.6, 1.8 ) );


  std::cout << bCurve << std::endl;

  const int nbr_points = 4;
  PointType data[nbr_points];
  data[0] = PointType::make_point(0.6, 1.2, 1.0);
  data[1] = PointType::make_point(1.3, 1.6, 1.8);
  data[2] = PointType::make_point(2.9, 2.4, 2.3);
  data[3] = PointType::make_point(3.2, 3.5, 3.0);
  BezierCurveType b2Curve(data, 4);
  std::cout << "Checking the control point constructor:" << std::endl;
  std::cout << b2Curve << std::endl;

  std::cout << "Checking indexing operator: " << std::endl;
  std::cout << "The final control points of the above two bezier curves are " << bCurve[1] << " and " << b2Curve[3] << "." << std::endl; 
 
  std::cout << "Checking the evaluation of bezier curves above: " << std::endl; 
  std::cout << "Curve 1 at t=0 is " << eval_bezier(bCurve,0.0) << " and Curve 2 at t=.5 is " << eval_bezier(b2Curve,.5) <<  std::endl;
  std::cout << "------------------End checking Bezier Functions---------------------" << std::endl;
  //_ctrlpts_end

  return bCurve;
}

int main(int argc, char** argv)
{

  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);

  testBezier();

  return 0;
}

