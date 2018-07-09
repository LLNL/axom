/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*! \file primal_introduction.cpp
 *  \brief This example code is a basic demonstration of how to use
 *  Primal to represent geometric primitives, perform geometric operations,
 *  and use a spatial index.
 */

/* This example code contains snippets used in the Primal Sphinx documentation.
 * They begin and end with comments
 *
 * prims_header_start
 * prims_header_end
 * clip_header_start
 * clip_header_end
 * closest_point_header_start
 * closest_point_header_end
 * bbox_header_start
 * bbox_header_end
 * using_start
 * using_end
 * clip_start
 * clip_end
 * closest_point_start
 * closest_point_end
 * bbox_start
 * bbox_end
 * naive_intersection_start
 * naive_intersection_end
 * ugrid_intersection_start
 * ugrid_intersection_end
 *
 * each prepended with an underscore.
 */

// _prims_header_start
// Axom primitives
#include "primal/BoundingBox.hpp"
#include "primal/Plane.hpp"
#include "primal/Point.hpp"
#include "primal/Polygon.hpp"
#include "primal/Ray.hpp"
#include "primal/Triangle.hpp"
// _prims_header_end

// Axom operations
// Each header is used in its own example, so each one is bracketed for separate inclusion.
// _clip_header_start
#include "primal/clip.hpp"
// _clip_header_end
// _closest_point_header_start
#include "primal/closest_point.hpp"
// _closest_point_header_end
// _bbox_header_start
#include "primal/compute_bounding_box.hpp"
// _bbox_header_end
// _intersect_header_start
#include "primal/intersect.hpp"
// _intersect_header_end
// _orient_header_start
#include "primal/orientation.hpp"
// _orient_header_end
// _sqdist_header_start
#include "primal/squared_distance.hpp"
// _sqdist_header_end

// Axom spatial index
#include "primal/UniformGrid.hpp"

// C++ headers
#include <cmath> // do we need this?

// _using_start
// "using" directives to simplify code
using namespace axom;
using namespace primal;

// all our examples are in 3D
const int NDIMS = 3;

// all our primitives are represented by doubles, in 3D
typedef Point<double, NDIMS> PointType;
typedef Triangle<double, NDIMS> TriangleType;
typedef BoundingBox<double, NDIMS> BoundingBoxType;
typedef Polygon<double, NDIMS> PolygonType;
// _using_end

PolygonType showClip()
{
  // _clip_start
  TriangleType tri(PointType::make_point(0.7,   0,   0),
                   PointType::make_point(  0, 1.8,   0),
                   PointType::make_point(  0,   0, 1.4));

  BoundingBoxType bbox(PointType::make_point(0, -0.5, 0),
                       PointType::make_point(1,    1, 1));

  PolygonType poly = clip(tri, bbox);
  // _clip_end

  std::cout << poly << std::endl;

  return poly;
}

void showClosestPoint()
{
  // _closest_point_start
  TriangleType tri(PointType::make_point(1, 0, 0),
                   PointType::make_point(0, 1, 0),
                   PointType::make_point(0, 0, 1));

  PointType pto = PointType::make_point( 0, 0, 0);
  PointType pta = PointType::make_point(-3, 2, 2);

  // Query point o lies at the origin.  Its closest point lies in the
  // interior of tri.
  PointType cpto = closest_point(pto, tri);

  // Query point a lies farther from the triangle.  Its closest point
  // is on tri's edge.
  int lcpta = 0;
  PointType cpta = closest_point(pta, tri, &lcpta);
  // _closest_point_end

  std::cout << pto << std::endl << cpto << std::endl << std::endl;
  std::cout << pta << std::endl << cpta << std::endl << lcpta << std::endl;
}

void showBoundingBoxes()
{
  // _bbox_start
  // An array of Points to include in the bounding boxes
  const int nbr_points = 6;
  PointType data[nbr_points];
  data[0] = PointType::make_point(0.6, 1.2, 1.0);
  data[1] = PointType::make_point(1.3, 1.6, 1.8);
  data[2] = PointType::make_point(2.9, 2.4, 2.3);
  data[3] = PointType::make_point(3.2, 3.5, 3.0);
  data[4] = PointType::make_point(3.6, 3.2, 4.0);
  data[5] = PointType::make_point(4.3, 4.3, 4.5);

  // A BoundingBox constructor takes an array of Point objects
  BoundingBoxType bbox(data, nbr_points);
  // Make an OrientedBoundingBox with the compute_oriented_bounding_box function
  OrientedBoundingBoxType obbox = compute_oriented_bounding_box(data, nbr_points);
  // _bbox_end

  for (int i = 0; i < nbr_points; ++i)
  {
    std::cout << data[i] << std::endl;
  }
  std::cout << bbox << std::endl;
  std::cout << obbox << std::endl;
}

int main(int argc, char** argv)
{

  // Deal with unused variables
  // _primitive_start
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);
  // _primitive_end

  int region[3375];

  DataStore* ds = create_datastore(region);
  access_datastore(ds);

  DataStore* tds = create_tiny_datastore();
  save_as_blueprint(tds);

  return 0;
}
