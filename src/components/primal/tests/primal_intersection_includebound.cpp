/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Ray.hpp"
#include "primal/Segment.hpp"
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"

#define AXOM_TRI_INTERSECTION_INCLUDES_BOUNDARY
#include "primal/intersection.hpp"

using namespace axom;

namespace {

double randomDouble(double beg = 0., double end = 1.)
{
  double range = end-beg;

  if (range == 0) {
    range = 1.;
  }

  return beg + (rand() / ( RAND_MAX / range ) );
}

template < int NDIMS >
primal::Point< double,NDIMS > randomPt( double beg, double end )
{
  primal::Point< double,NDIMS > pt;
  for ( int i=0; i< NDIMS; ++i ) {
    pt[i] = randomDouble(beg,end);
  }

  return pt;
}

}

template < int DIM >
primal::Triangle< double, DIM > roll(const primal::Triangle< double, DIM > & t,
                                     const int i)
{
  return primal::Triangle< double, DIM > (t[i % 3], t[(i+1) % 3], t[(i+2) % 3]);
}

template < int DIM >
void permuteCornersTest(const primal::Triangle< double, DIM > & a,
                        const primal::Triangle< double, DIM > & b,
                        const std::string & whattest,
                        const bool testtrue)
{
  SCOPED_TRACE(whattest);

  bool allmatch = true;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      allmatch = allmatch &&
        (primal::intersect(roll(a, i), roll(b, j)) == testtrue);
    }
  }

  const primal::Triangle< double, DIM > ap(a[0], a[2], a[1]);
  const primal::Triangle< double, DIM > bp(b[0], b[2], b[1]);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      allmatch = allmatch &&
        (primal::intersect(roll(ap, i), roll(bp, j)) == testtrue);
    }
  }

  // Now swap a for b

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      allmatch = allmatch &&
        (primal::intersect(roll(b, i), roll(a, j)) == testtrue);
    }
  }

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      allmatch = allmatch &&
        (primal::intersect(roll(bp, i), roll(ap, j)) == testtrue);
    }
  }

  if (allmatch) {
    SUCCEED();
  } else {
    ADD_FAILURE() << (testtrue?
                      "Triangles should intersect but did not" :
                      "Triangles should not intersect but did");
  }
}

TEST( primal_intersection_includebound, 2D_triangle_triangle_intersection )
{

  typedef primal::Triangle< double,2 > Triangle2;
  typedef primal::Point< double,2 >    Point2;

  // Triangle 569
  Triangle2 triA( Point2::make_point(0.0,5.0),
                  Point2::make_point(5.0,5.0),
                  Point2::make_point( 0.0,0.0) );

  Triangle2 triB( Point2::make_point(0.0,5.0),
                  Point2::make_point(5.0,5.0),
                  Point2::make_point(0.0,0.0) );

  // axom::slic::setLoggingMsgLevel( axom::slic::message::Info);

  // Several intersection cases (and one non-intersection)

  permuteCornersTest(triA, triB, "identical 2D triangles", true);

  Triangle2 triC( Point2::make_point(-1.0,-1.0),
                  Point2::make_point(-5.0,-5.0),
                  Point2::make_point(-7.0,-8.0) );

  permuteCornersTest(triA, triC, "non-intersecting 2D triangles", false);

  triA = Triangle2( Point2::make_point(4.3,4.05),
                    Point2::make_point(-1.0,-0.06),
                    Point2::make_point( 7.3, -1.3) );

  triB = Triangle2( Point2::make_point(1.0, 0.0),
                    Point2::make_point(6.0,0.5),
                    Point2::make_point( 4.2,  2.1) );

  permuteCornersTest(triA, triB, "2D tri B completely contained in tri A",
                     true);

  triB = Triangle2( Point2::make_point(1.9,-2),
                    Point2::make_point(6.9,2.1),
                    Point2::make_point(0.8,5.1) );

  permuteCornersTest(triA, triB, "intersecting 2D triangles, no corner in",
                     true);

  triB = Triangle2( Point2::make_point(2.9,1.6),
                    Point2::make_point(-1.5,1.5),
                    Point2::make_point(0.8,5.1) );

  permuteCornersTest(triA, triB, "intersecting 2D triangles, one corner in",
                     true);

  triB = Triangle2( Point2::make_point(2.9,0),
                    Point2::make_point(2.1,0.1),
                    Point2::make_point(0.8,5.1) );

  permuteCornersTest(triA, triB, "intersecting 2D triangles, two corners in",
                     true);

  triB = Triangle2( Point2::make_point(2, -1),
                    Point2::make_point(-1.0,-0.06),
                    Point2::make_point(7.3,-1.3) );

  permuteCornersTest(triA, triB,
                     "2D t1 and t2 share a complete edge (and nothing else)",
                     true);

  Triangle2 triD( Point2::make_point(0, 0),
                  Point2::make_point(1,0),
                  Point2::make_point(1, 1) );

  Triangle2 triE( Point2::make_point(0, 0),
                  Point2::make_point(0.5,0),
                  Point2::make_point(-1, -1) );

  permuteCornersTest(triD, triE,
                     "2D t1 edge is a subset of t2's, and they share a corner (but nothing else)",
                     true);

  triE = Triangle2( Point2::make_point(0.5, 0),
                    Point2::make_point(1,0),
                    Point2::make_point(-1, -1) );

  permuteCornersTest(triD, triE,
                     "2D t1 edge is a subset of t2's, and they share the other corner (but nothing else)",
                     true);

  triE = Triangle2( Point2::make_point(0.5, 0),
                    Point2::make_point(1.5,0),
                    Point2::make_point(-1,-1) );

  permuteCornersTest(triD, triE,
                     "2D t1 edge overlaps t2 (no other intersection)",
                     true);

  triE = Triangle2( Point2::make_point(-0.5, 0),
                    Point2::make_point(0.5,0),
                    Point2::make_point(-1,-1) );

  permuteCornersTest(triD, triE,
                     "2D t1 edge overlaps t2 the other way (no other intersection)",
                     true);

  triE = Triangle2( Point2::make_point(-1, 0.5),
                    Point2::make_point(-1,-1),
                    Point2::make_point(2, -1) );

  permuteCornersTest(triD, triE,
                     "2D t1 point lands on t2 edge (no other intersection)",
                     true);

  triE = Triangle2( Point2::make_point(0, 0),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE,
                     "2D t1 point lands on t2 point (no other intersection)",
                     true);

  // Several non-intersection cases (and a few intersection)

  triE = Triangle2( Point2::make_point(0.2, -1e-3),
                    Point2::make_point(1,-1),
                    Point2::make_point(1.2, -1e-3) );

  permuteCornersTest(triD, triE, "2D disjunct, close parallel sides",
                     false);

  triE = Triangle2( Point2::make_point(0.2, -1e-3),
                    Point2::make_point(1,-1),
                    Point2::make_point(1, -1e-4) );

  permuteCornersTest(triD, triE, "2D disjunct, close converging sides",
                     false);

  triE = Triangle2( Point2::make_point(10, 1),
                    Point2::make_point(2,0),
                    Point2::make_point(11, -0.3) );

  permuteCornersTest(triD, triE, "2D disjunct, fairly far-separated",
                     false);

  triE = Triangle2( Point2::make_point(0, 0.1),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D disjunct, point comes close",
                     false);

  triE = Triangle2( Point2::make_point(-0.001, 0),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 2",
                     false);

  triE = Triangle2( Point2::make_point(-0.5, 0),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 3",
                     false);

  triE = Triangle2( Point2::make_point(-1.7, 0),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 4",
                     false);

  triE = Triangle2( Point2::make_point(-5.1, 0),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D disjunct, point comes close 5",
                     false);

  triE = Triangle2( Point2::make_point(0.5, 0.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point lands on side 2", true);

  triE = Triangle2( Point2::make_point(0.49999, 0.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side",
                     false);

  triE = Triangle2( Point2::make_point(0.49, 0.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side 2",
                     false);

  triE = Triangle2( Point2::make_point(0.4, 0.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side 3",
                     false);

  triE = Triangle2( Point2::make_point(-0.1, 0.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side 4",
                     false);

  triE = Triangle2( Point2::make_point(-2.6, 2.5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side 5",
                     false);

  triE = Triangle2( Point2::make_point(-6, 5),
                    Point2::make_point(-40,-0.7),
                    Point2::make_point(-23, 1.3) );

  permuteCornersTest(triD, triE, "2D point comes close to side 6",
                     false);
}

bool makeTwoRandomIntersecting3DTriangles(primal::Triangle< double, 3 > & l,
                                          primal::Triangle< double, 3 > & r)
{
  typedef primal::Triangle< double,3 > Triangle3;
  typedef primal::Point< double,3 >   Point3;
  typedef primal::Vector< double, 3 > Vector3;

  //Step 1: Construct a random triangle
  Point3 A= randomPt< 3 >(0.,1.);
  Point3 B= randomPt< 3 >(0.,1.);
  Point3 C= randomPt< 3 >(0.,1.);
  l = Triangle3(A,B,C);

  //Step 2: Construct two random points on the triangle.
  Point3 P;
  Point3 Q;

  double a1= randomDouble();
  double a2= randomDouble();
  double a3= randomDouble();

  double n1= (a1/(a1+a2+a3));
  double n2= (a2/(a1+a2+a3));
  double n3= (a3/(a1+a2+a3));

  double P_x= n1*A[0]+n2*B[0]+n3*C[0];
  double P_y= n1*A[1]+n2*B[1]+n3*C[1];
  double P_z= n1*A[2]+n2*B[2]+n3*C[2];
  P= Point3::make_point(P_x,P_y,P_z);

  a1= randomDouble();
  a2= randomDouble();
  a3= randomDouble();

  n1= (a1/(a1+a2+a3));
  n2= (a2/(a1+a2+a3));
  n3= (a3/(a1+a2+a3));

  double Q_x= n1*A[0]+n2*B[0]+n3*C[0];
  double Q_y= n1*A[1]+n2*B[1]+n3*C[1];
  double Q_z= n1*A[2]+n2*B[2]+n3*C[2];

  Q= Point3::make_point(Q_x,Q_y,Q_z);

  /*PQ is so random segment on the triangle.  We create a vertex called vertex1 and
     use
     it to create the triangle formed by P', Q' and vertex1. */

  //Step 3: choose some vertex away from the triangle
  Point3 vertex1 = randomPt< 3 >(0.,1.);

  //Step 4:
  //we scale the segments formed by both vertex 1 and P and by vertex 1 and Q so that
  // we now
  //have a triangle whose base is not necessarily on the plane formed by ABC
  Vector3 vertex2Direction = Vector3(Q, vertex1);
  Vector3 vertex3Direction = Vector3(P, vertex1);

  //construct the other two vertices of the triangle
  Point3 vertex2 =
    Point3::make_point(vertex1[0]-2*vertex2Direction[0],
                       vertex1[1]-2*vertex2Direction[1],
                       vertex1[2]-2*vertex2Direction[2]);
  Point3 vertex3 =
    Point3::make_point(vertex1[0]-2*vertex3Direction[0],
                       vertex1[1]-2*vertex3Direction[1],
                       vertex1[2]-2*vertex3Direction[2]);

  Triangle3 before = Triangle3(vertex1, Q, P);
  r = Triangle3(vertex1, vertex2, vertex3);

  return !l.degenerate() && !r.degenerate();
}

TEST( primal_intersection_includebound, 3D_triangle_triangle_intersection )
{

  typedef primal::Triangle< double,3 > Triangle3;
  typedef primal::Point< double,3 >   Point3;

  Triangle3 tri3d_1( Point3::make_point(-1.0,-1.0,-1.0),
                     Point3::make_point(-2.0,-5.0,-5.0),
                     Point3::make_point(-4.0,-8.0, -8.0) );

  Triangle3 tri3d_2( Point3::make_point(-1.0,-1.0,-1.0),
                     Point3::make_point(-2.0,-5.0,-5.0),
                     Point3::make_point(-4.0,-8.0, -8.0) );

  permuteCornersTest(tri3d_1, tri3d_2, "3D identical triangles", true);

  Triangle3 tri3d_3( Point3::make_point(1.0,1.0,1.0),
                     Point3::make_point(5.0,5.0,5.0),
                     Point3::make_point(8.0,7.0, 92.0) );

  permuteCornersTest(tri3d_1, tri3d_3, "3D disjunct triangles", false);

  Triangle3 tri3A( Point3::make_point(0, 0, 0),
                   Point3::make_point(1, 0,0),
                   Point3::make_point(0, 1.7, 2.3) );

  Triangle3 tri3B( Point3::make_point(0, 0, 0),
                   Point3::make_point(1, 0,0),
                   Point3::make_point(0, -2, 1.2) );

  permuteCornersTest(tri3A, tri3B, "3D tris sharing a segment",
                     true);

  tri3B = Triangle3( Point3::make_point(-0.2, 0, 0),
                     Point3::make_point(0.7, 0,0),
                     Point3::make_point(0, -2, 1.2) );

  permuteCornersTest(tri3A, tri3B, "3D tris sharing part of a segment",
                     true);

  tri3B = Triangle3( Point3::make_point(-1, 0, 0),
                     Point3::make_point(0, 4.3,6),
                     Point3::make_point(0, 1.7, 2.3) );

  permuteCornersTest(tri3A, tri3B, "3D tris sharing a vertex",
                     true);

  tri3B = Triangle3( Point3::make_point(0, -1, 0),
                     Point3::make_point(1, 1,0),
                     Point3::make_point(0, 1.7, -2.3) );

  permuteCornersTest(tri3A, tri3B, "3D tris; edges cross",
                     true);

  tri3B = Triangle3( Point3::make_point(0, -1, -1),
                     Point3::make_point(0.5, 0,0),
                     Point3::make_point(1, 1, -1) );

  permuteCornersTest(tri3A, tri3B, "3D tris; B vertex lands on A's edge",
                     true);

  tri3B = Triangle3( Point3::make_point(0.5, -1, 0.1),
                     Point3::make_point(0.5, 1,0.1),
                     Point3::make_point(1, 1, -1) );

  permuteCornersTest(tri3A, tri3B,
                     "3D tris intersect like two links in a chain", true);

  tri3B = Triangle3( Point3::make_point(-1, -1, 1),
                     Point3::make_point(0, 2,1),
                     Point3::make_point(5, 0, 1) );

  permuteCornersTest(tri3A, tri3B, "3D tri A pokes through B",
                     true);

  tri3B = Triangle3( Point3::make_point(1, -1, 1),
                     Point3::make_point(1, 2,1),
                     Point3::make_point(1, 0, -1) );

  permuteCornersTest(tri3A, tri3B, "3D tri A vertex tangent on B",
                     true);

  tri3B = Triangle3( Point3::make_point(1.00001, -1, 1),
                     Point3::make_point(1, 2,1),
                     Point3::make_point(1, 0, -1) );

  permuteCornersTest(tri3A, tri3B, "3D tri A vertex not quite tangent on B",
                     false);

  // 3D versions of 2D test cases (!)

  //future work: test triangle triangle with a bunch of random test cases

  srand(1);   //we want same random number sequence everytime to make sure our tests
              // don't differ on a case to case basis

  //Randomly generate a bunch of intersecting triangles (whose intersections form
  // segments) and test them
  // How many tests are we actually performing here?
  int rantests = 0;
  int skiptests = 0;
  for (int i=0; i<5000; i++) {
    Triangle3 randomTriangle, intersectingTriangle;

    if (makeTwoRandomIntersecting3DTriangles(randomTriangle,
                                             intersectingTriangle)) {
      permuteCornersTest(randomTriangle, intersectingTriangle, "random", true);
      rantests += 1;
    }
    else {
      skiptests += 1;
    }
  }

  SLIC_INFO( "Ran " << rantests << " and skipped " << skiptests <<
             " tests due to triangle degeneracy." );

  axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
