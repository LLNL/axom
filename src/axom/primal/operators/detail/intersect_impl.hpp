// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect_impl.hpp
 *
 * This file provides helper functions for testing whether
 * geometric primitives intersect
 */

#ifndef AXOM_PRIMAL_INTERSECT_IMPL_HPP_
#define AXOM_PRIMAL_INTERSECT_IMPL_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Triangle.hpp"

namespace axom
{
namespace primal
{
namespace detail
{
//--------------------TYPE ALIASES AND FUNCTION DECLARATIONS ------------------

using Vector3 = primal::Vector<double, 3>;
using Point3 = primal::Point<double, 3>;
using Triangle3 = primal::Triangle<double, 3>;
using Triangle2 = primal::Triangle<double, 2>;
using Point2 = primal::Point<double, 2>;

AXOM_HOST_DEVICE
bool isGt(double x, double y, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool isLt(double x, double y, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool isLeq(double x, double y, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool isLpeq(double x, double y, bool includeEqual = false, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool isGeq(double x, double y, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool isGpeq(double x, double y, bool includeEqual = false, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool nonzeroSignMatch(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool twoZeros(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
bool oneZeroOthersMatch(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
int countZeros(double x, double y, double z, double EPS = 1E-12);

AXOM_HOST_DEVICE
double twoDcross(const Point2& A, const Point2& B, const Point2& C);

/*!
 * This function finds where p1 lies in relation to the vertices of t2
 * and calls either checkEdge() or checkVertex().
 * \return status true iff triangle p1 q1 r1 intersects triangle p2 q2 r2.
 */
AXOM_HOST_DEVICE
bool intersectPermuted2DTriangles(const Point2& p1,
                                  const Point2& q1,
                                  const Point2& r1,
                                  const Point2& p2,
                                  const Point2& q2,
                                  const Point2& r2,
                                  bool includeBoundary,
                                  double EPS);

/*!
 * Triangle 1 vertices have been permuted to CCW: permute t2 to CCW
 * and call worker function to test for intersection.
 *
 * q1 and r1 both lie in the negative half-space defined by t2; p1 lies in
 * t2's plane or in its positive half-space.
 * The sign of dp2, dq2, and dr2 indicates whether the associated vertex
 * of t2 lies in the positive or negative half-space defined by t1.
 */
AXOM_HOST_DEVICE
bool intersectOnePermutedTriangle(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  double dp2,
                                  double dq2,
                                  double dr2,
                                  Vector3& normal,
                                  bool includeBoundary,
                                  double EPS);

AXOM_HOST_DEVICE
bool intersectTwoPermutedTriangles(const Point3& p1,
                                   const Point3& q1,
                                   const Point3& r1,
                                   const Point3& p2,
                                   const Point3& q2,
                                   const Point3& r2,
                                   bool includeBoundary,
                                   double EPS);

/*!
 * Project (nearly) coplanar triangles 1 and 2 on an axis; call 2D worker
 * function to test for intersection.
 */
AXOM_HOST_DEVICE
bool intersectCoplanar3DTriangles(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  Vector3 normal,
                                  bool includeBoundary,
                                  double EPS);

/*!
 * \brief Orients triangle vertices so both triangles are CCW; calls worker.
 * \return status true iff t1 and t2 intersect
 *
 * Determine triangle orientation, then call the worker function with
 * vertices from t1 and t2 permuted to ensure CCW orientation.
 */
AXOM_HOST_DEVICE
bool TriangleIntersection2D(const Triangle2& t1,
                            const Triangle2& t2,
                            bool includeBoundary,
                            double EPS);

//------------------------------ IMPLEMENTATIONS ------------------------------

/*! @{ @name 3D triangle-triangle intersection */

/*!
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * This algorithm is modeled after Devillers and Guigue (2002).  It computes
 * the line of intersection L of the triangles' planes and finds the segments
 * S1, S2 where t1 and t2 intersect L.  If those segments intersect, the
 * triangles must intersect.  Note that edge and vertex intersections are
 * reported as hits: an edge or vertex intersecting any part of another
 * triangle causes a return value of true.  This is consistent with the paper.
 *
 * Coplanar triangles are handled with a decision tree, testing the
 * relative position of triangle vertices.

 * Olivier Devillers and Phillipe Guigue, Faster Triangle-Triangle Intersection
 * Tests, RR-4488, INRIA (2002).  https://hal.inria.fr/inria-00072100/
 */
template <typename T>
AXOM_HOST_DEVICE bool intersect_tri3D_tri3D(const Triangle<T, 3>& t1,
                                            const Triangle<T, 3>& t2,
                                            bool includeBoundary,
                                            double EPS)
{
  typedef primal::Vector<T, 3> Vector3;

  SLIC_CHECK_MSG(!t1.degenerate(),
                 "\n\n WARNING \n\n Triangle " << t1 << " is degenerate");
  SLIC_CHECK_MSG(!t2.degenerate(),
                 "\n\n WARNING \n\n Triangle " << t2 << " is degenerate");

  // Step 1: Check if all the vertices of triangle 1 lie on the same side of
  // the plane created by triangle 2:

  // Vector3 t2Normal = Vector3::cross_product(Vector3(t2[2], t2[0]),
  //                                           Vector3(t2[2], t2[1]));
  Vector3 t2Normal = t2.normal().unitVector();
  double dp1 = (Vector3(t2[2], t1[0])).dot(t2Normal);
  double dq1 = (Vector3(t2[2], t1[1])).dot(t2Normal);
  double dr1 = (Vector3(t2[2], t1[2])).dot(t2Normal);

  if(nonzeroSignMatch(dp1, dq1, dr1, EPS))
  {
    return false;
  }

  if(!includeBoundary &&
     (twoZeros(dp1, dq1, dr1, EPS) || oneZeroOthersMatch(dp1, dq1, dr1, EPS)))
  {
    return false;
  }

  // Step 2: Check if all the vertices of triangle 2 lie on the same side of
  // the plane created by triangle 1:

  // Vector3 t1Normal = Vector3::cross_product(Vector3(t1[0], t1[1]),
  //                                           Vector3(t1[0], t1[2]));
  Vector3 t1Normal = t1.normal().unitVector();
  double dp2 = (Vector3(t1[2], t2[0])).dot(t1Normal);
  double dq2 = (Vector3(t1[2], t2[1])).dot(t1Normal);
  double dr2 = (Vector3(t1[2], t2[2])).dot(t1Normal);

  if(nonzeroSignMatch(dp2, dq2, dr2, EPS))
  {
    return false;
  }

  if(!includeBoundary &&
     (twoZeros(dp2, dq2, dr2, EPS) || oneZeroOthersMatch(dp2, dq2, dr2, EPS)))
  {
    return false;
  }

  /* Note: Because we know that all the vertices either triangle do not
     lay on the same side of the plane formed by the other triangle, we
     know that for each triangle, exactly 1 out of 3 points exists on one
     side of the plane formed by the other triangle.


     Step 3: We apply a circular permutation of triangle 1 such that its
     first point is the only point on the triangle that lies on one side of
     the plane formed by triangle 2 (with the other 2 on the other side),
     while handling the special case of one of the vertices lying on the
     plane formed by triangle 2.  We then perform a swap operation on the
     second and third points of triangle 2 to map the first point of
     triangle 1 to the positive halfspace formed by triangle 2's plane.
   */

  // compare the signs to create a convenient permutation of the vertices
  // of triangle 1

  if(isGt(dp1, 0.0, EPS))
  {
    if(isGt(dq1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[2],
                                          t1[0],
                                          t1[1],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else if(isGt(dr1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[1],
                                          t1[2],
                                          t1[0],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else
    {
      return intersectOnePermutedTriangle(t1[0],
                                          t1[1],
                                          t1[2],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
  }
  else if(isLt(dp1, 0.0, EPS))
  {
    if(isLt(dq1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[2],
                                          t1[0],
                                          t1[1],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else if(isLt(dr1, 0.0, EPS))
    {
      return intersectOnePermutedTriangle(t1[1],
                                          t1[2],
                                          t1[0],
                                          t2[0],
                                          t2[1],
                                          t2[2],
                                          dp2,
                                          dq2,
                                          dr2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
    else
    {
      return intersectOnePermutedTriangle(t1[0],
                                          t1[1],
                                          t1[2],
                                          t2[0],
                                          t2[2],
                                          t2[1],
                                          dp2,
                                          dr2,
                                          dq2,
                                          t1Normal,
                                          includeBoundary,
                                          EPS);
    }
  }
  else  //dp1 ~= 0
  {
    if(isLt(dq1, 0.0, EPS))
    {
      if(isGeq(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[1],
                                            t1[2],
                                            t1[0],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectOnePermutedTriangle(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
    else if(isGt(dq1, 0.0, EPS))
    {
      if(isGt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectOnePermutedTriangle(t1[1],
                                            t1[2],
                                            t1[0],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
    else
    {
      if(isGt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[2],
                                            t1[0],
                                            t1[1],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            dp2,
                                            dq2,
                                            dr2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else if(isLt(dr1, 0.0, EPS))
      {
        return intersectOnePermutedTriangle(t1[2],
                                            t1[0],
                                            t1[1],
                                            t2[0],
                                            t2[2],
                                            t2[1],
                                            dp2,
                                            dr2,
                                            dq2,
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
      else
      {
        return intersectCoplanar3DTriangles(t1[0],
                                            t1[1],
                                            t1[2],
                                            t2[0],
                                            t2[1],
                                            t2[2],
                                            t1Normal,
                                            includeBoundary,
                                            EPS);
      }
    }
  }
}

/*!
 * \brief Tests for general 3D triangle-triangle intersection.
 * \return status true iff triangles 1 and 2 intersect.
 *
 * Previous tests have ruled out cases where planes of triangles 1 and 2
 * are parallel or identical, as well as cases where plane 1 does not
 * intersect triangle 2 (and vice versa).  The vertices have been permuted
 * so p1 is across plane 2 from q1 and r1, and p2 is across plane 1
 * from q2 and r2.
 *
 * The core of the Devillers and Guigue's method examines the line l0
 * defined by the intersection of
 * planes 1 and 2.  Assume triangles 1 and 2 both intersect this line, which
 * they must if they intersect each other.  Then, the intersection of triangle
 * 1 with l0 is a segment s1, and the intersection of triangle 2 with l0 is
 * a segment s2.  If s1 and s2 overlap, triangles 1 and 2 intersect.
 *
 * This function implements Equation 1 from Devillers and Guigue (2002), p.8
 * using a hint from p.10 that greatly simplifies the computation.
 *
 * Helper function for T-T intersect.
 */
inline bool intersectTwoPermutedTriangles(const Point3& p1,
                                          const Point3& q1,
                                          const Point3& r1,
                                          const Point3& p2,
                                          const Point3& q2,
                                          const Point3& r2,
                                          bool includeBoundary,
                                          double EPS)
{
  /* Step 5: From step's 1 through 4, we now have two triangles that,
     if intersecting, have a line that intersects segments p1r1, p1q1,
     p2q2, and p2r2.  We check if these two intervals overlap:
   */
  const bool bdr = includeBoundary;

  return isLpeq(Vector3(q1, q2).dot(Triangle3(q1, p2, p1).normal()), 0.0, bdr, EPS) &&
    isLpeq(Vector3(p1, r2).dot(Triangle3(p1, p2, r1).normal()), 0.0, bdr, EPS);
}

/*! @} */

/*! @{ @name 2D triangle-triangle intersection */
/*
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * Note that edge and vertex intersections are reported as hits: an edge or
 * vertex intersecting any part of another triangle causes a return value of
 * true.  This is consistent with the paper.
 */
template <typename T>
bool intersect_tri2D_tri2D(const primal::Triangle<T, 2>& t1,
                           const primal::Triangle<T, 2>& t2,
                           bool includeBoundary,
                           double EPS)
{
  SLIC_CHECK_MSG(!t1.degenerate(),
                 "\n\n WARNING \n\n Triangle " << t1 << " is degenerate");
  SLIC_CHECK_MSG(!t2.degenerate(),
                 "\n\n WARNING \n\n Triangle " << t2 << " is degenerate");

  return TriangleIntersection2D(t1, t2, includeBoundary, EPS);
}

/*!
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2p2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 */
AXOM_HOST_DEVICE
inline bool checkEdge(const Point2& p1,
                      const Point2& q1,
                      const Point2& r1,
                      const Point2& p2,
                      const Point2& r2,
                      bool includeBoundary,
                      double EPS);

/*!
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 */
AXOM_HOST_DEVICE
inline bool checkVertex(const Point2& p1,
                        const Point2& q1,
                        const Point2& r1,
                        const Point2& p2,
                        const Point2& q2,
                        const Point2& r2,
                        bool includeBoundary,
                        double EPS);

/*!
 * \brief Compute cross product of two 2D vectors as if they were 3D.
 * \return Cross product of A C and B C.
 *
 * This function treats three Point2 values as corners of a 3D triangle with
 * zero Z-coordinate.  Thus we can calculate the cross product of A C with
 * B C using only the k-hat term, since the other terms go to zero.  A
 * positive value indicates CCW orientation.
 *
 * \note The result is equal to twice the signed area of a 2D triangle
 * with vertices (A,B,C) (in CCW order).
 */
inline double twoDcross(const Point2& A, const Point2& B, const Point2& C)
{
  return (((A[0] - C[0]) * (B[1] - C[1]) - (A[1] - C[1]) * (B[0] - C[0])));
}

/*!
 * \brief Checks if x > y, within a specified tolerance.
 */
inline bool isGt(double x, double y, double EPS)
{
  return ((x > y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/*!
 * \brief Checks if x < y, within a specified tolerance.
 */
inline bool isLt(double x, double y, double EPS)
{
  return ((x < y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/*!
 * \brief Checks if x <= y, within a specified tolerance.
 */
inline bool isLeq(double x, double y, double EPS) { return !(isGt(x, y, EPS)); }

/*!
 * \brief Checks if x < y, or possibly x == y, within a specified tolerance.
 *
 * The check for equality is controlled by parameter includeEqual.  This
 * lets users specify at compile time whether triangles intersecting only on
 * border points are reported as intersecting or not.
 *
 * Supports checkEdge and checkVertex
 */
inline bool isLpeq(double x, double y, bool includeEqual, double EPS)
{
  if(includeEqual && axom::utilities::isNearlyEqual(x, y, EPS))
  {
    return true;
  }

  return isLt(x, y, EPS);
}

/*!
 * \brief Checks if x >= y, within a specified tolerance.
 */
inline bool isGeq(double x, double y, double EPS) { return !(isLt(x, y, EPS)); }

/*!
 * \brief Checks if x > y, or possibly x == y, within a specified tolerance.
 *
 * The check for equality is controlled by parameter includeEqual.  This
 * lets users specify at compile time whether triangles intersecting only on
 * border points are reported as intersecting or not.
 *
 * Supports checkEdge and checkVertex
 */
inline bool isGpeq(double x, double y, bool includeEqual, double EPS)
{
  if(includeEqual && axom::utilities::isNearlyEqual(x, y, EPS))
  {
    return true;
  }

  return isGt(x, y, EPS);
}

/*!
 * \brief Check if x, y, and z all have the same sign.
 */
inline bool nonzeroSignMatch(double x, double y, double z, double EPS)
{
  return !(axom::utilities::isNearlyEqual(x, 0., EPS)) &&
    !(axom::utilities::isNearlyEqual(y, 0., EPS)) &&
    !(axom::utilities::isNearlyEqual(z, 0., EPS)) &&
    (0 < x) - (x < 0) == (0 < y) - (y < 0) &&
    (0 < x) - (x < 0) == (0 < z) - (z < 0);
}

/*!
 * \brief Check if two of x, y, and z are near zero.
 */
inline bool twoZeros(double x, double y, double z, double EPS)
{
  return countZeros(x, y, z, EPS) == 2;
}

/*!
 * \brief Check if one of x, y, and z is near zero and the others' signs match.
 */
inline bool oneZeroOthersMatch(double x, double y, double z, double EPS)
{
  namespace util = axom::utilities;
  return countZeros(x, y, z, EPS) == 1 &&
    ((util::isNearlyEqual(x, 0.0, EPS) && isGt(y * z, 0.0, EPS)) ||
     (util::isNearlyEqual(y, 0.0, EPS) && isGt(z * x, 0.0, EPS)) ||
     (util::isNearlyEqual(z, 0.0, EPS) && isGt(x * y, 0.0, EPS)));
}

/*!
 * \brief Count the number of arguments near zero.
 */
inline int countZeros(double x, double y, double z, double EPS)
{
  return (int)axom::utilities::isNearlyEqual(x, 0.0, EPS) +
    (int)axom::utilities::isNearlyEqual(y, 0.0, EPS) +
    (int)axom::utilities::isNearlyEqual(z, 0.0, EPS);
}

/*! @} */

/*!
 * \brief Helper function to find disjoint projections for the AABB-triangle
 * test
 * \param d0 The first value defining the test interval
 * \param d1 The second value defining the test interval
 * \param d2 The third value defining the test interval
 * \param r Radius of projection
 * \return True of the intervals are disjoint, false otherwise
 */
bool intervalsDisjoint(double d0, double d1, double d2, double r);
bool crossEdgesDisjoint(double d0, double d1, double r);

/*! @{ @name Triangle-bbox intersection */

/*!
 * \brief Determines if a triangle and a bounding box intersect
 *        (but does not find the intersections)
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 */
template <typename T>
bool intersect_tri_bbox(const primal::Triangle<T, 3>& tri,
                        const primal::BoundingBox<T, 3>& bb)
{
  // Note: Algorithm is derived from the one presented in chapter 5.2.9 of
  //   Real Time Collision Detection book by Christer Ericson
  // based on Akenine-Moller algorithm (Journal of Graphics Tools)
  //
  // It uses the Separating Axis Theorem to look for disjoint projections
  // along various axes associated with Faces and Edges of the AABB and
  // triangle.
  // There are 9 tests for the cross products of edges
  //           3 tests for the AABB face normals
  //           1 test for the triangle face normal
  // We use early termination if we find a separating axis between the shapes

  typedef typename BoundingBox<T, 3>::PointType PointType;
  typedef typename BoundingBox<T, 3>::VectorType VectorType;

  // Extent: vector center to max corner of BB
  VectorType e = 0.5 * bb.range();

  // Make the AABB center the origin by moving the triangle vertices
  PointType center(bb.getMin().array() + e.array());
  VectorType v[3] = {VectorType(center, tri[0]),
                     VectorType(center, tri[1]),
                     VectorType(center, tri[2])};

  // Create the edge vectors of the triangle
  VectorType f[3] = {v[1] - v[0], v[2] - v[1], v[0] - v[2]};

  /* clang-format off */

  // Test cross products of edges between triangle edge vectors f and cube normals (9 tests)
  // -- using separating axis theorem on the cross product of edges of triangle and face normals of AABB
  // Each test involves three cross products, two of which have the same value
  // The commented parameters highlights this symmetry.
#define XEDGE_R( _0, _1, _I )      e[ _0 ] * std::abs( f[ _I ][ _1 ]) + e[ _1 ] * std::abs(f[ _I ][ _0 ])
#define XEDGE_S( _0, _1, _V, _F ) -v[ _V ][ _0 ] * f[ _F ][ _1 ] + v[ _V ][ _1 ] * f[ _F ][ _0 ]

  if ( crossEdgesDisjoint(/*XEDGE_S(1,2,0,0),*/ XEDGE_S(1,2,1,0),   XEDGE_S(1,2,2,0),   XEDGE_R(1,2,0)) ||
       crossEdgesDisjoint(  XEDGE_S(1,2,0,1),/* XEDGE_S(1,2,1,1),*/ XEDGE_S(1,2,2,1),   XEDGE_R(1,2,1)) ||
       crossEdgesDisjoint(  XEDGE_S(1,2,0,2),   XEDGE_S(1,2,1,2),/* XEDGE_S(1,2,2,2),*/ XEDGE_R(1,2,2)) ||
       crossEdgesDisjoint(/*XEDGE_S(2,0,0,0),*/ XEDGE_S(2,0,1,0),   XEDGE_S(2,0,2,0),   XEDGE_R(0,2,0)) ||
       crossEdgesDisjoint(  XEDGE_S(2,0,0,1),/* XEDGE_S(2,0,1,1),*/ XEDGE_S(2,0,2,1),   XEDGE_R(0,2,1)) ||
       crossEdgesDisjoint(  XEDGE_S(2,0,0,2),   XEDGE_S(2,0,1,2),/* XEDGE_S(2,0,2,2),*/ XEDGE_R(0,2,2)) ||
       crossEdgesDisjoint(/*XEDGE_S(0,1,0,0),*/ XEDGE_S(0,1,1,0),   XEDGE_S(0,1,2,0),   XEDGE_R(0,1,0)) ||
       crossEdgesDisjoint(  XEDGE_S(0,1,0,1),/* XEDGE_S(0,1,1,1),*/ XEDGE_S(0,1,2,1),   XEDGE_R(0,1,1)) ||
       crossEdgesDisjoint(  XEDGE_S(0,1,0,2),   XEDGE_S(0,1,1,2),/* XEDGE_S(0,1,2,2),*/ XEDGE_R(0,1,2)) )
  {
    return false;
  }
  /* clang-format on */

#undef XEDGE_R
#undef XEDEG_S

  /// Test face normals of bounding box (3 tests)
  if(intervalsDisjoint(v[0][0], v[1][0], v[2][0], e[0]) ||
     intervalsDisjoint(v[0][1], v[1][1], v[2][1], e[1]) ||
     intervalsDisjoint(v[0][2], v[1][2], v[2][2], e[2]))
  {
    return false;
  }

  /// Final test -- face normal of triangle's plane
  VectorType planeNormal = VectorType::cross_product(f[0], f[1]);
  double planeDist = planeNormal.dot(tri[0]);

  double r = e[0] * std::abs(planeNormal[0]) + e[1] * std::abs(planeNormal[1]) +
    e[2] * std::abs(planeNormal[2]);
  double s = planeNormal.dot(center) - planeDist;

  return std::abs(s) <= r;
}

// ------------------------------------------------------------------------------

/*!
 * \brief Helper function for Triangle/BoundingBox intersection test
 */
inline bool crossEdgesDisjoint(double d0, double d1, double r)
{
  return axom::utilities::max(-axom::utilities::max(d0, d1),
                              axom::utilities::min(d0, d1)) > r;
}

/*! @} */

/*! @{ @name Triangle-ray intersection */

/*!
 * \brief Tests if 3D triangle tri intersects with 3D ray R.
 * \param [in] tri The input triangle
 * \param [in] R The input ray
 * \param [out] t Intersection point of tri and R, w.r.t. parametrization of R
 * \param [out] p Intersection point of tri and R, in **un-normalized**
 *   barycentric coordinates relative to tri.  To normalize, divide each
 *   component by the sum of the components.
 * \note If there is an intersection, the intersection point pt is:
 *                     pt = R.at(t)
 * \note If R is coplanar with tri, this routine will report a miss.  If R's
 *       origin lies within tri, this routine will report a miss.
 * \return status true iff tri intersects with R, otherwise, false.
 *
 * This algorithm is modeled after Woop, Benthin, and Wald (2013).  It
 * transforms the ray onto the unit z vector and applies the same transform
 * to the triangle's vertices.  This transform simplifies the calculation of
 * barycentric coordinates of the positive z-axis's intersection with the
 * triangle.  If any of these coordinates are less than zero, the ray misses.
 * Note that for efficiency, the barycentric coordinates of the intersection
 * are not normalized.  To normalize, divide each coordinate by the sum of
 * the coordinates.
 *
 * If any of the barycentric coordinates are equal to zero, more care is
 * needed to check if the ray hits the edge or misses.
 *
 * Sven Woop, Carsten Benthin, Ingo Wald, "Watertight Ray/Triangle
 * Intersection," Journal of Computer Graphics Techniques (JCGT), vol. 2,
 * no. 1, 65–82, 2013  http://jcgt.org/published/0002/01/05/
 */
template <typename T>
bool intersect_tri_ray(const Triangle<T, 3>& tri,
                       const Ray<T, 3>& R,
                       T& t,
                       Point<double, 3>& p)
{
  // Ray origins inside of the triangle are considered a miss.
  // This is a good thing, as pointed out by Matt Larsen in January 2017,
  // because of a common technique in ray chasing through a mesh:
  //  1. Cast the ray until it intersects a triangle (record the hit)
  //  2. Starting from the hit point, re-cast the ray to find the next hit
  // Origin-is-miss avoids finding the already-hit triangle as the next
  // hit.
  //
  // Coplanar rays, evidenced by det == 0, are also considered as a miss.
  // I (Arlie Capps, Jan. 2017) don't understand the motivation at this
  // point, but I'll accept this for now.

  typedef NumericArray<T, 3> NumArray;
  const T zero = T();

  //find out dimension where ray direction is maximal
  int kx, ky, kz;

  NumArray r = primal::abs(R.direction().array());

  //z-direction largest
  if((r[2] >= r[0]) && (r[2] >= r[1]))
  {
    kz = 2;
  }
  //y direction largest
  else if((r[1] >= r[0]) && (r[1] >= r[2]))
  {
    kz = 1;
  }
  //x direction largest
  else
  {
    kz = 0;
  }

  //assign other dimensions of the ray
  kx = (kz + 1) % 3;
  ky = (kz + 2) % 3;

  //if necessary swap  ky and kx to preserve triangle winding
  if(R.direction()[kz] < zero)
  {
    axom::utilities::swap(kx, ky);
  }

  //calculate shear constants
  NumericArray<T, 3> shear(1.0f / R.direction()[kz], 3);
  shear[0] *= R.direction()[kx];
  shear[1] *= R.direction()[ky];

  //A,B,C are the triangle vertices offset to Ray's origin
  NumArray A = tri[0].array() - R.origin().array();
  NumArray B = tri[1].array() - R.origin().array();
  NumArray C = tri[2].array() - R.origin().array();

  //shear and scale the vertices
  const T Ax = A[kx] - shear[0] * A[kz];
  const T Ay = A[ky] - shear[1] * A[kz];
  const T Bx = B[kx] - shear[0] * B[kz];
  const T By = B[ky] - shear[1] * B[kz];
  const T Cx = C[kx] - shear[0] * C[kz];
  const T Cy = C[ky] - shear[1] * C[kz];

  //scaled barycentric coordinates
  p[0] = Cx * By - Cy * Bx;
  p[1] = Ax * Cy - Ay * Cx;
  p[2] = Bx * Ay - By * Ax;
  const T& U = p[0];
  const T& V = p[1];
  const T& W = p[2];

  //edge testing
  if((U < zero || V < zero || W < zero) && (U > zero || V > zero || W > zero))
  {
    return false;
  }

  //calculate determinant
  const T det = U + V + W;

  if(det == zero)
  {
    return false;
  }

  //calculate scaled z-coordinates of the vertices and use them to calculate hit
  // distance
  const T Az = shear[2] * A[kz];
  const T Bz = shear[2] * B[kz];
  const T Cz = shear[2] * C[kz];
  t = (U * Az + V * Bz +
       W * Cz);  // save the parameter of the intersection w.r.t. ray
                 // R

  //make sure hit is in correct direction
  if(((t < zero) && !(det < zero)) || ((det < zero) && !(t < zero)))
  {
    return false;
  }

  t /= det;
  return true;
}

/*! @} */

/*!
 * \brief Tests if 3D triangle tri intersects with 3D ray S.
 * \param [in] tri The input triangle
 * \param [in] S The input segment
 * \param [out] t Intersection point of tri and S, w.r.t. parametrization of S
 * \param [out] p Intersection point of tri and S, in barycentric coordinates
 *   relative to tri
 * \return status true iff tri intersects with R, otherwise, false.
 *
 * This routine uses intersect_tri_ray(), which see.
 */
template <typename T>
bool intersect_tri_segment(const Triangle<T, 3>& tri,
                           const Segment<T, 3>& S,
                           T& t,
                           Point<double, 3>& p)
{
  typedef Vector<T, 3> Vector3;
  Ray<T, 3> r(S.source(), Vector3(S.source(), S.target()));

  //Ray-triangle intersection does not check endpoints, so we explicitly check
  // here
  if(tri.checkInTriangle(S.source()))
  {
    t = 0;
    p = tri.physToBarycentric(S.source());
    return true;
  }
  if(tri.checkInTriangle(S.target()))
  {
    t = 1;
    p = tri.physToBarycentric(S.target());
    return true;
  }

  // The triangle only intersects the segment if it intersects the ray defined
  // by one
  // of its endpoints and the direction defined by its two endpoints.
  // We can parametrize the line as:  r.origin() + t * r.direction()
  // Values of the parameter t between 0 and the length of the segment
  // correspond
  // to points on the segment.
  // Note: if intersect_tri_ray() is true, t must be greater than zero
  if(intersect_tri_ray(tri, r, t, p))
  {
    t = t / static_cast<T>(S.length());
    return t <= 1;
  }
  return false;
}

/*!
 * \brief Determines if a 1D OBB intersects a 1D OBB.
 * \param [in] b1 A 1D OrientedBoundingBox
 * \param [in] b2 A 1D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb1D_obb1D(const OrientedBoundingBox<T, 1>& b1,
                           const OrientedBoundingBox<T, 1>& b2)
{
  T c1 = b1.getCentroid()[0];
  T c2 = b2.getCentroid()[0];

  T e1 = b1.getExtents()[0];
  T e2 = b2.getExtents()[0];

  if(c1 + e1 > c2 - e2)
  {
    return true;
  }
  if(c2 + e2 > c1 - e1)
  {
    return true;
  }

  return false;
}

/*!
 * \brief Determines if a 2D OBB intersects a 2D OBB.
 * \param [in] b1 A 2D OrientedBoundingBox
 * \param [in] b2 A 2D OrientedBoundingBox
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb2D_obb2D(const OrientedBoundingBox<T, 2>& b1,
                           const OrientedBoundingBox<T, 2>& b2)
{
  Vector<T, 2> c1(b1.getCentroid());
  Vector<T, 2> c2(b2.getCentroid());

  Vector<T, 2> e1 = b1.getExtents();
  Vector<T, 2> e2 = b2.getExtents();

  const Vector<T, 2>* u1 = b1.getAxes();
  const Vector<T, 2>* u2 = b2.getAxes();

  Vector<T, 2> d = c2 - c1;

  for(int i = 0; i < 2; ++i)
  {
    if(utilities::abs<T>(d.dot(u1[i])) > e1[i] +
         utilities::abs<T>((e2[0] * u2[0]).dot(u1[i])) +
         utilities::abs<T>((e2[1] * u2[1]).dot(u1[i])))
    {
      return false;
    }
  }

  for(int i = 0; i < 2; ++i)
  {
    if(utilities::abs<T>(d.dot(u2[i])) > e2[i] +
         utilities::abs<T>((e1[0] * u1[0]).dot(u2[i])) +
         utilities::abs<T>((e1[1] * u1[1]).dot(u2[i])))
    {
      return false;
    }
  }

  return true;
}

/*!
 * \brief Determines if a 3D OBB intersects a 3D OBB.
 * \param [in] b1 A 3D OrientedBoundingBox
 * \param [in] b2 A 3D OrientedBoundingBox
 * \param [in] EPS error tolerance for intersection
 * \return true iff b1 intersects with b2, otherwise, false.
 */
template <typename T>
bool intersect_obb3D_obb3D(const OrientedBoundingBox<T, 3>& b1,
                           const OrientedBoundingBox<T, 3>& b2,
                           double EPS)
{
  Vector<T, 3> d =
    Vector<T, 3>(b1.getCentroid()) - Vector<T, 3>(b2.getCentroid());

  Vector<T, 3> e1 = b1.getExtents();
  Vector<T, 3> e2 = b2.getExtents();

  const Vector<T, 3>* u1 = b1.getAxes();
  const Vector<T, 3>* u2 = b2.getAxes();

  // compute r and r^T here:
  Vector<T, 3> r[3];
  Vector<T, 3> rt[3];
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      r[i][j] = utilities::abs<T>(u1[i].dot(u2[j]));
      rt[i][j] = utilities::abs<T>(u1[j].dot(u2[i]));
    }
  }

  // check for separating planes parallel to faces
  for(int i = 0; i < 3; ++i)
  {
    if(utilities::abs<T>(d.dot(u1[i])) > e1[i] + e2.dot(r[i]) + EPS)
    {
      return false;
    }
  }

  for(int i = 0; i < 3; ++i)
  {
    if(utilities::abs<T>(d.dot(u2[i])) > e2[i] + e1.dot(rt[i]) + EPS)
    {
      return false;
    }
  }

  // check for separating planes with normals parallel to cross product of edges
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      T left = utilities::abs<T>(d.dot(u1[(i + 2) % 3]) * r[(i + 1) % 3][j] -
                                 (d.dot(u1[(i + 1) % 3]) * r[(i + 2) % 3][j]));
      T right = (e1[(i + 1) % 3] * r[(i + 2) % 3][j] +
                 e1[(i + 2) % 3] * r[(i + 1) % 3][j]);
      right += (e2[(i + 1) % 3] * r[i][(j + 2) % 3] +
                e2[(i + 2) % 3] * r[i][(j + 1) % 3]);
      if(left > right + EPS)
      {
        return false;
      }
    }
  }

  // didn't find a separating anything
  return true;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif  // AXOM_PRIMAL_INTERSECT_IMPL_HPP_
