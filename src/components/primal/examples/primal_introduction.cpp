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
#include "primal/OrientedBoundingBox.hpp"
#include "primal/Plane.hpp"
#include "primal/Point.hpp"
#include "primal/Polygon.hpp"
#include "primal/Ray.hpp"
#include "primal/Segment.hpp"
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
#include <iostream>
#include <fstream>

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
typedef OrientedBoundingBox<double, NDIMS> OrientedBoundingBoxType;
typedef Polygon<double, NDIMS> PolygonType;
typedef Ray<double, NDIMS> RayType;
typedef Segment<double, NDIMS> SegmentType;
typedef Vector<double, NDIMS> VectorType;
// _using_end
typedef Plane<double, NDIMS> PlaneType;

PolygonType showClip()
{
  // _clip_start
  TriangleType tri(PointType::make_point(1.2,   0,   0),
                   PointType::make_point(  0, 1.8,   0),
                   PointType::make_point(  0,   0, 1.4));

  BoundingBoxType bbox(PointType::make_point(0, -0.5, 0),
                       PointType::make_point(1,    1, 1));

  PolygonType poly = clip(tri, bbox);
  // _clip_end

  std::cout << poly << std::endl;

  // Now write out an Asymptote file showing what we did.
  std::string fname = "showClip.asy";
  std::ofstream asy(fname);
  if (!asy.good()) {
    std::cout << "Could not write to " << fname << std::endl;
  } else {
    asy << "// To turn this Asymptote source file into an image for inclusion in\n"
           "// Axom's documentation,\n"
           "// 1. run Asymptote:\n"
           "//    asy -f png " << fname << std::endl <<
           "// 2. Optionally, use ImageMagick to convert the white background to transparent:\n"
           "//    convert " << fname << " -transparent white " << fname <<
      std::endl << std::endl;
    asy << "// preamble" << std::endl;
    asy << "settings.render = 6;" << std::endl;
    asy << "import three;" << std::endl;
    asy << "size(6cm, 0);" << std::endl << std::endl;

    asy << "// axes" << std::endl;
    asy << "draw(O -- 1.7X, arrow=Arrow3(DefaultHead2), " <<
      "L=Label(\"$x$\", position=EndPoint, align=W));" << std::endl;
    asy << "draw(O -- 2.4Y, arrow=Arrow3(), " <<
      "L=Label(\"$y$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2Z, arrow=Arrow3(), " <<
      "L=Label(\"$z$\", position=EndPoint));" << std::endl << std::endl;

    asy << "// polygon" << std::endl;
    asy << "path3 pgon = ";
    for (int i = 0; i < poly.numVertices(); ++i) {
      asy << poly[i] << "--";
    }
    asy << "cycle;" << std::endl << std::endl;

    asy << "// triangle" << std::endl;
    asy << "path3 tri = " << tri[0] << "--" << tri[1] << "--" <<
      tri[2] << "--cycle;" << std::endl << std::endl;

    asy << "// draw triangle then polygon" << std::endl;
    asy << "draw(surface(tri), surfacepen=blue+opacity(0.4));" << std::endl;
    asy << "draw(tri);" << std::endl << std::endl;
    asy << "draw(surface(pgon), surfacepen=yellow+opacity(0.4));" << std::endl;
    asy << "draw(pgon, yellow);" << std::endl << std::endl;

    asy << "// bounding box" << std::endl;
    asy << "draw(box(" << bbox.getMin() << ", " << bbox.getMax() <<
      "));" << std::endl;
  }

  return poly;
}

void showClosestPoint()
{
  // _closest_point_start
  TriangleType tri(PointType::make_point(1, 0, 0),
                   PointType::make_point(0, 1, 0),
                   PointType::make_point(0, 0, 1));

  PointType pto = PointType::make_point( 0, 0, 0);
  PointType pta = PointType::make_point(-1, 2, 1);

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

  // Now write out an Asymptote file showing what we did.
  // Helper: the XY plane.
  double normal[3] = {0, 0, 1};
  PlaneType XYPlane(normal, 0.);
  // Projected points
  double ppta[3];
  double pcpta[3];
  double pcpto[3];
  XYPlane.projectPoint(pta.array().data(),  ppta);
  XYPlane.projectPoint(cpta.array().data(), pcpta);
  XYPlane.projectPoint(cpto.array().data(), pcpto);
  std::string fname = "showClosestPoint.asy";
  std::ofstream asy(fname);
  if (!asy.good()) {
    std::cout << "Could not write to " << fname << std::endl;
  } else {
    asy << "// To turn this Asymptote source file into an image for inclusion in\n"
           "// Axom's documentation,\n"
           "// 1. run Asymptote:\n"
           "//    asy -f png " << fname << std::endl <<
           "// 2. Optionally, use ImageMagick to convert the white background to transparent:\n"
           "//    convert " << fname << " -transparent white " << fname <<
      std::endl << std::endl;
    asy << "// preamble" << std::endl;
    asy << "settings.render = 6;" << std::endl;
    asy << "import three;" << std::endl;
    asy << "size(6cm, 0);" << std::endl << std::endl;

    asy << "// axes" << std::endl;
    asy << "draw(-4.5X -- 1.7X, arrow=Arrow3(DefaultHead2), " <<
      "L=Label(\"$x$\", position=EndPoint, align=W));" << std::endl;
    asy << "draw(O -- 2.4Y, arrow=Arrow3(), " <<
      "L=Label(\"$y$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2Z, arrow=Arrow3(), " <<
      "L=Label(\"$z$\", position=EndPoint));" << std::endl << std::endl;

    asy << "// triangle" << std::endl;
    asy << "path3 tri = " << tri[0] << "--" << tri[1] << "--" <<
      tri[2] << "--cycle;" << std::endl << std::endl;

    asy << "// triangle" << std::endl;
    asy << "triple pto = " << pto << ";" << std::endl;
    asy << "triple pta = " << pta << ";" << std::endl;
    asy << "triple cpto = " << cpto << ";" << std::endl;
    asy << "triple cpta = " << cpta << ";" << std::endl;
    asy << "triple ppta = (" << ppta[0] << "," << ppta[1] << "," <<
      ppta[2] << ");" << std::endl;
    asy << "triple pcpto = (" << pcpto[0] << "," << pcpto[1] << "," <<
      pcpto[2] << ");" << std::endl;
    asy << "triple pcpta = (" << pcpta[0] << "," << pcpta[1] << "," <<
      pcpta[2] << ");" << std::endl << std::endl;

    asy << "// draw triangle then points and projections" << std::endl;
    asy << "draw(tri);" << std::endl;
    asy << "dot(pto, blue);" << std::endl;
    asy << "label(\"$o$\", pto, align=W);" << std::endl;
    asy << "dot(cpto, mediumblue);" << std::endl;
    asy << "label(\"$o'$\", cpto, align=N);" << std::endl;
    asy << "draw(cpto--pcpto, dotted);" << std::endl;
    asy << "dot(pta, lightolive);" << std::endl;
    asy << "label(\"$a$\", pta, align=W);" << std::endl;
    asy << "draw(pta--ppta, dotted);" << std::endl;
    asy << "dot(cpta, yellow);" << std::endl;
    asy << "label(\"$a'$\", cpta, align=NE);" << std::endl;
    asy << "draw(cpta--pcpta, dotted);" << std::endl;
  }
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

  // Now write out an Asymptote file showing what we did.
  std::string fname = "showBoundingBoxes.asy";
  std::ofstream asy(fname);
  if (!asy.good()) {
    std::cout << "Could not write to " << fname << std::endl;
  } else {
     asy << "// To turn this Asymptote source file into an image for inclusion in\n"
           "// Axom's documentation,\n"
           "// 1. run Asymptote:\n"
           "//    asy -f png " << fname << std::endl <<
           "// 2. Optionally, use ImageMagick to convert the white background to transparent:\n"
           "//    convert " << fname << " -transparent white " << fname <<
      std::endl << std::endl;
   asy << "// preamble" << std::endl;
    asy << "settings.render = 6;" << std::endl;
    asy << "import three;" << std::endl;
    asy << "size(6cm, 0);" << std::endl << std::endl;

    asy << "// projection" << std::endl;
    asy << "currentprojection = perspective((4, -1.8, 3), (0.07, 0.07, 1));" <<
      std::endl << std::endl;

    asy << "// axes" << std::endl;
    asy << "draw(O -- 4X, arrow=Arrow3(DefaultHead2), " <<
      "L=Label(\"$x$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 7Y, arrow=Arrow3(), " <<
      "L=Label(\"$y$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 5Z, arrow=Arrow3(), " <<
      "L=Label(\"$z$\", position=EndPoint));" << std::endl << std::endl;

    asy << "// points" << std::endl;
    asy << "triple[] points = new triple[6];" << std::endl;
    for (int i = 0; i < nbr_points; ++i) {
      asy << "points[" << i << "] = " << data[i] << ";" << std::endl;
    }
    asy << std::endl;

    asy << "// bbox" << std::endl;
    asy << "triple bboxmin = " << bbox.getMin() << ";" << std::endl;
    asy << "triple bboxmax = " << bbox.getMax() << ";" << std::endl;
    asy << std::endl;

    asy << "// oriented bounding box" << std::endl;
    asy << "triple[] obpts = new triple[8];" << std::endl;
    std::vector< PointType > obboxpts = obbox.vertices();
    for (int i = 0; i < 8; ++i) {
      asy << "obpts[" << i << "] = " << obboxpts[i] << ";" << std::endl;
    }
    asy << std::endl;

    asy << "// draw points" << std::endl;
    for (int i = 0; i < nbr_points; ++i) {
      asy << "dot(points[" << i << "], blue);" << std::endl;
    }
    asy << std::endl;

    asy << "// draw bbox" << std::endl;
    asy << "draw(box(bboxmin, bboxmax));" << std::endl << std::endl;

    asy << "// draw oriented bounding box" << std::endl;
    asy << "path3[] obboxpath = obpts[0]--obpts[1]--obpts[3]--obpts[2]--cycle"
      << std::endl << "     ^^ obpts[4]--obpts[5]--obpts[7]--obpts[6]--cycle"
      << std::endl << "    ";
    for (int i = 0; i < 4; ++i) {
      asy << " ^^ obpts[" << i << "]--obpts[" << (i+4) << "]";
    }
    asy << ";" << std::endl;
    asy << "draw(obboxpath, orange);" << std::endl << std::endl;
  }

  for (int i = 0; i < nbr_points; ++i)
  {
    std::cout << data[i] << std::endl;
  }
  std::cout << bbox << std::endl;
  std::cout << obbox << std::endl;
}

void showIntersect()
{
  // _intersect_start
  // Two triangles
  TriangleType tri1(PointType::make_point(1.2,   0,   0),
                    PointType::make_point(  0, 1.8,   0),
                    PointType::make_point(  0,   0, 1.4));

  TriangleType tri2(PointType::make_point(  0,   0, 0.5),
                    PointType::make_point(0.8, 0.1, 1.2),
                    PointType::make_point(0.8, 1.4, 1.2));

  // tri1 and tri2 should intersect
  if (intersect(tri1, tri2)) {
    std::cout << "Triangles intersect as expected." << std::endl;
  } else {
    std::cout << "There's an error somewhere..." << std::endl;
  }

  // A vertical ray constructed from origin and point
  RayType ray(SegmentType(PointType::make_point(0.4, 0.4, 0),
                          PointType::make_point(0.4, 0.4, 1)));

  // t will hold the intersection point between ray and tri1,
  // as parameterized along ray.
  double rt1t = 0;
  // rt1b will hold the intersection point barycentric coordinates,
  // and rt1p will hold the physical coordinates.
  PointType rt1b, rt1p;

  // The ray should intersect tri1 and tri2.
  if (intersect(tri1, ray, rt1t, rt1b) && intersect(tri2, ray)) {
    rt1p = tri1.baryToPhysical(rt1b);
    std::cout << "Ray intersects tri1 as expected.  Parameter t: " <<
      rt1t << std::endl << "  Intersect barycentric coordinates: " << rt1b <<
      std::endl << "  Intersect physical coordinates: " << rt1p << std::endl <<
      "Ray also intersects tri2 as expected." << std::endl;
  } else {
    std::cout << "There's an error somewhere..." << std::endl;
  }
  
  // A bounding box
  BoundingBoxType bbox(PointType::make_point(0.1, -0.23, 0.1),
                       PointType::make_point(0.8,  0.5,  0.4));

  // The bounding box should intersect tri1 and ray but not tr2.
  PointType bbtr1;
  if (intersect(ray, bbox, bbtr1) && intersect(tri1, bbox) &&
      !intersect(tri2, bbox)) {
    std::cout << "As hoped, bounding box intersects tri1 at " << bbtr1 <<
      " and ray, but not tri2." << std::endl;
  } else {
    std::cout << "There is at least one error somewhere..." << std::endl;
  }
  // _intersect_end

  // helper variables
  PolygonType poly = clip(tri1, bbox);
  // These are parametric coordinates ray-tri2-t, tri1-tri2 leg A-t, tri1-tri2 leg C-t
  // and corresponding physical points.
  double rt2t, t1t2at, t1t2ct;
  (void)intersect(tri2, ray, rt2t);
  PointType rt2p = ray.at(rt2t);
  SegmentType t2lega(tri2[0], tri2[1]);
  (void)intersect(tri1, t2lega, t1t2at);
  PointType t1t2ap = t2lega.at(t1t2at);
  SegmentType t2legc(tri2[2], tri2[0]);
  (void)intersect(tri1, t2legc, t1t2ct);
  PointType t1t2cp = t2legc.at(t1t2ct);
  // Project point C of tri2 onto the XY plane
  double normal[3] = {0, 0, 1};
  PlaneType XYPlane(normal, 0.);
  PointType tr2c = tri2[2];
  double pp[3];
  XYPlane.projectPoint(tr2c.array().data(), pp);

  // Now write out an Asymptote file showing what we did.
  std::string fname = "showIntersect.asy";
  std::ofstream asy(fname);
  if (!asy.good()) {
    std::cout << "Could not write to " << fname << std::endl;
  } else {
    asy << "// To turn this Asymptote source file into an image for inclusion in\n"
           "// Axom's documentation,\n"
           "// 1. run Asymptote:\n"
           "//    asy -f png " << fname << std::endl <<
           "// 2. Optionally, use ImageMagick to convert the white background to transparent:\n"
           "//    convert " << fname << " -transparent white " << fname <<
      std::endl << std::endl;
    asy << "// preamble" << std::endl;
    asy << "settings.render = 6;" << std::endl;
    asy << "import three;" << std::endl;
    asy << "size(6cm, 0);" << std::endl << std::endl;

    asy << "// axes" << std::endl;
    asy << "draw(O -- 1.7X, arrow=Arrow3(DefaultHead2), " <<
      "L=Label(\"$x$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2.4Y, arrow=Arrow3(), " <<
      "L=Label(\"$y$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2Z, arrow=Arrow3(), " <<
      "L=Label(\"$z$\", position=EndPoint, align=W));" << std::endl << std::endl;

    asy << "// triangle 1" << std::endl;
    asy << "path3 tri1 = " << tri1[0] << "--" << tri1[1] << "--" <<
      tri1[2] << "--cycle;" << std::endl << std::endl;

    asy << "// triangle 2" << std::endl;
    asy << "path3 tri2 = " << tri2[0] << "--" << tri2[1] << "--" <<
      tri2[2] << "--cycle;" << std::endl << std::endl;

    asy << "// ray" << std::endl;
    asy << "path3 ray = " << ray.origin() << "--" << ray.at(1.8) <<
      ";" << std::endl << std::endl;

    asy << "// polygon of intersection between bbox and triangle" << std::endl;
    asy << "path3 pgon = ";
    for (int i = 0; i < poly.numVertices(); ++i) {
      asy << poly[i] << "--";
    }
    asy << "cycle;" << std::endl << std::endl;

    asy << "// draw bounding box and other geometry" << std::endl;
    asy << "draw(box(" << bbox.getMin() << ", " << bbox.getMax() <<
      "), blue);" << std::endl;
    asy << "draw(pgon, deepblue);" << std::endl << std::endl;
    asy << "draw(ray, arrow=Arrow3(DefaultHead2), red);" << std::endl;
    asy << "dot(" << bbtr1 << ", red);" << std::endl << "dot(" << rt1p <<
      ", red);" << std::endl << "dot(" << rt2p << ", red);" << std::endl;
    asy << "draw(tri1);" << std::endl << "draw(tri2, blue);" << std::endl;
    asy << "draw(" << t1t2ap << "--" << t1t2cp << ", deepblue);" << std::endl;
    asy << "draw(" << tr2c << "--(" << pp[0] << "," << pp[1] << "," << 
      pp[2] <<"), dotted);" << std::endl;
  }
}

void showOrientation()
{
  // _orient_start
  // A triangle
  TriangleType tri(PointType::make_point(1.2,   0,   0),
                   PointType::make_point(  0, 1.8,   0),
                   PointType::make_point(  0,   0, 1.4));

  // Three points:
  //    one on the triangle's positive side,
  PointType pos = PointType::make_point(0.45, 1.5, 1);
  //    one coplanar to the triangle, the centroid,
  PointType cpl = PointType::lerp(PointType::lerp(tri[0], tri[1], 0.5),
                                   tri[2], 1./3.);
  //    and one on the negative side
  PointType neg = PointType::make_point(0, 0, 0.7);

  // Test orientation
  if (orientation(pos, tri)  == ON_POSITIVE_SIDE &&
      orientation(cpl, tri) == ON_BOUNDARY &&
      orientation(neg, tri)  == ON_NEGATIVE_SIDE) {
    std::cout << "As expected, point pos is on the positive side," <<
      std::endl << "    point cpl is on the boundary (on the triangle)," <<
      std::endl << "    and point neg is on the negative side." << std::endl;
  } else {
    std::cout << "Someone wrote this wrong." << std::endl;
  }
  // _orient_end

  // Helper variables
  // Project onto the XY plane
  PointType ppos = PointType::make_point(pos[0], pos[1], 0.);
  PointType pcpl = PointType::make_point(cpl[0], cpl[1], 0.);

  // Now write out an Asymptote file showing what we did.
  std::string fname = "showOrientation.asy";
  std::ofstream asy(fname);
  if (!asy.good()) {
    std::cout << "Could not write to " << fname << std::endl;
  } else {
    asy << "// To turn this Asymptote source file into an image for inclusion in\n"
           "// Axom's documentation,\n"
           "// 1. run Asymptote:\n"
           "//    asy -f png " << fname << std::endl <<
           "// 2. Optionally, use ImageMagick to convert the white background to transparent:\n"
           "//    convert " << fname << " -transparent white " << fname <<
      std::endl << std::endl;
    asy << "// preamble" << std::endl;
    asy << "settings.render = 6;" << std::endl;
    asy << "import three;" << std::endl;
    asy << "size(6cm, 0);" << std::endl << std::endl;

    asy << "// axes" << std::endl;
    asy << "draw(O -- 1.7X, arrow=Arrow3(DefaultHead2), " <<
      "L=Label(\"$x$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2.4Y, arrow=Arrow3(), " <<
      "L=Label(\"$y$\", position=EndPoint));" << std::endl;
    asy << "draw(O -- 2Z, arrow=Arrow3(), " <<
      "L=Label(\"$z$\", position=EndPoint, align=W));" << std::endl << std::endl;

    asy << "// triangle" << std::endl;
    asy << "path3 tri = " << tri[0] << "--" << tri[1] << "--" <<
      tri[2] << "--cycle;" << std::endl << std::endl;
    asy << "triple centroid = " << cpl << ";" << std::endl;

    asy << "draw(tri);" << std::endl;
    asy << "dot(" << neg << ", blue);" << std::endl;
    asy << "dot(" << cpl << ", blue);" << std::endl;
    asy << "draw(centroid--1.6centroid, arrow=Arrow3(DefaultHead2));" << std::endl;
    asy << "dot(" << pos << ", blue);" << std::endl;
    asy << "draw(" << pos << "--" << ppos << ", dotted);" << std::endl;
  }
}

int main(int argc, char** argv)
{

  // Deal with unused variables
  // _primitive_start
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);
  // _primitive_end

  showClip();
  showClosestPoint();
  showBoundingBoxes();
  showIntersect();
  showOrientation();

  return 0;
}
