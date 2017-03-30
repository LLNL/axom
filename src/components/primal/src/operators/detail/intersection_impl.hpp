/*!
 *******************************************************************************
 * \file Intersection.hpp
 *
 * This file provides several functions to test whether geometric primitives
 * intersect.
 *******************************************************************************
 */

#ifndef INTERSECTION_IMPL_HPP_
#define INTERSECTION_IMPL_HPP_

#include "primal/BoundingBox.hpp"
#include "primal/Determinants.hpp"
#include "primal/Point.hpp"
#include "primal/Ray.hpp"
#include "primal/Segment.hpp"
#include "primal/Triangle.hpp"

#include "axom_common/Utilities.hpp"

namespace axom {
namespace primal {
namespace detail {

// ================================== FORWARD DECLARATIONS =================

typedef primal::Vector< double, 3 > Vector3;
typedef primal::Point< double, 3 > Point3;
typedef primal::Triangle< double, 3 > Triangle3;
typedef primal::Triangle< double, 2 > Triangle2;
typedef primal::Point< double, 2 > Point2;

bool isGt(double x, double y, double EPS=1.0e-12);
bool isLt(double x, double y, double EPS=1.0e-12);
bool isLeq(double x, double y, double EPS=1.0e-12);
bool isGeq(double x, double y, double EPS=1.0e-12);
bool signMatch(double x, double y, double z, double EPS=1.0e-12);
double twoDcross(const Point2& A, const Point2& B, const Point2& C);

bool checkEdge(const Point2 p1,
               const Point2 q1,
               const Point2 r1,
               const Point2 p2,
               const Point2 r2);
bool checkVertex(const Point2 p1,
                 const Point2 q1,
                 const Point2 r1,
                 const Point2 p2,
                 const Point2 q2,
                 const Point2 r2);

bool intersectPermuted2DTriangles(const Point2& p1,
                                  const Point2& q1,
                                  const Point2& r1,
                                  const Point2& p2,
                                  const Point2& q2,
                                  const Point2& r2);

bool intersectOnePermutedTriangle(
  const Point3 &p1, const Point3 &q1, const Point3 &r1,
  const Point3 &p2, const Point3 &q2, const Point3 &r2,
  double dp2, double dq2, double dr2,  Vector3 &normal);

bool intersectTwoPermutedTriangles(const Point3 p1,
                                   const Point3 q1,
                                   const Point3 r1,
                                   const Point3 p2,
                                   const Point3 q2,
                                   const Point3 r2);

bool intersectCoplanar3DTriangles(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  Vector3 normal);

bool TriangleIntersection2D(const Triangle2& t1,
                            const Triangle2& t2);

/** @{ @name 3D triangle-triangle intersection */

/*!
 *******************************************************************************
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
 *******************************************************************************
 */
template < typename T >
bool intersect_tri3D_tri3D( const Triangle< T, 3 >& t1, const Triangle< T,
                                                                        3 >& t2)
{
  typedef primal::Vector< T, 3 > Vector3;

  SLIC_CHECK_MSG(
    !t1.degenerate(),
    "\n\n WARNING \n\n Triangle " << t1 <<" is degenerate");
  SLIC_CHECK_MSG(
    !t2.degenerate(),
    "\n\n WARNING \n\n Triangle " << t2 <<" is degenerate");

  // Step 1: Check if all the vertices of triangle 1 lie on the same side of
  // the plane created by triangle 2:

  // Vector3 t2Normal = Vector3::cross_product(Vector3(t2[2], t2[0]),
  //                                           Vector3(t2[2], t2[1]));
  Vector3 t2Normal = t2.normal();
  double dp1 = (Vector3(t2[2],t1[0])).dot(t2Normal);
  double dq1 = (Vector3(t2[2],t1[1])).dot(t2Normal);
  double dr1 = (Vector3(t2[2],t1[2])).dot(t2Normal);
  if (signMatch(dp1, dq1, dr1)) {
    return false;
  }

  // Step 2: Check if all the vertices of triangle 2 lie on the same side of
  // the plane created by triangle 1:

  // Vector3 t1Normal = Vector3::cross_product(Vector3(t1[0], t1[1]),
  //                                           Vector3(t1[0], t1[2]));
  Vector3 t1Normal = t1.normal();
  double dp2 = (Vector3(t1[2],t2[0])).dot(t1Normal);
  double dq2 = (Vector3(t1[2],t2[1])).dot(t1Normal);
  double dr2 = (Vector3(t1[2],t2[2])).dot(t1Normal);
  if (signMatch(dp2, dq2, dr2)) {
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

  if (isGt(dp1, 0.0)) {
    if (isGt(dq1, 0.0)) {
      return intersectOnePermutedTriangle(t1[2], t1[0], t1[1],
                                          t2[0], t2[2], t2[1],
                                          dp2, dr2, dq2, t1Normal);
    }
    else if (isGt(dr1, 0.0)) {
      return intersectOnePermutedTriangle(t1[1], t1[2], t1[0],
                                          t2[0], t2[2], t2[1],
                                          dp2, dr2, dq2, t1Normal);
    }
    else{return intersectOnePermutedTriangle(t1[0], t1[1], t1[2],
                                             t2[0], t2[1], t2[2],
                                             dp2, dq2, dr2, t1Normal);
    }
  }
  else if (isLt(dp1, 0.0)) {
    if (isLt(dq1, 0.0)) {
      return intersectOnePermutedTriangle(t1[2], t1[0], t1[1],
                                          t2[0], t2[1], t2[2],
                                          dp2, dq2, dr2, t1Normal);
    }
    else if (isLt(dr1, 0.0f)) {
      return intersectOnePermutedTriangle(t1[1], t1[2], t1[0],
                                          t2[0], t2[1], t2[2],
                                          dp2, dq2, dr2, t1Normal);
    }
    else{return intersectOnePermutedTriangle(t1[0], t1[1], t1[2],
                                             t2[0], t2[2], t2[1],
                                             dp2, dr2, dq2, t1Normal);
    }
  }
  else { //dp1 ~= 0
    if (isLt(dq1, 0.0)) {
      if (isGeq(dr1, 0.0)) {
        return intersectOnePermutedTriangle(t1[1], t1[2], t1[0],
                                            t2[0], t2[2], t2[1],
                                            dp2, dr2, dq2, t1Normal);
      }
      else {
        return intersectOnePermutedTriangle(t1[0], t1[1], t1[2],
                                            t2[0], t2[1], t2[2],
                                            dp2, dq2, dr2, t1Normal);
      }
    }
    else if (isGt(dq1, 0.0)) {
      if (isGt(dr1, 0.0)) {
        return intersectOnePermutedTriangle(t1[0], t1[1], t1[2],
                                            t2[0], t2[2], t2[1],
                                            dp2, dr2, dq2, t1Normal);
      }
      else {
        return intersectOnePermutedTriangle(t1[1], t1[2], t1[0],
                                            t2[0], t2[1], t2[2],
                                            dp2, dq2, dr2, t1Normal);
      }
    }
    else {
      if (isGt(dr1, 0.0)) {
        return intersectOnePermutedTriangle(t1[2], t1[0], t1[1],
                                            t2[0], t2[1], t2[2],
                                            dp2, dq2, dr2, t1Normal);
      }
      else if (isLt(dr1, 0.0)) {
        return intersectOnePermutedTriangle(t1[2], t1[0], t1[1],
                                            t2[0], t2[2], t2[1],
                                            dp2, dr2, dq2, t1Normal);
      }
      else{return intersectCoplanar3DTriangles(t1[0], t1[1], t1[2],
                                               t2[0], t2[1], t2[2], t1Normal);
      }
    }
  }
}

/*!
 *****************************************************************************
 * Triangle 1 vertices have been permuted to CCW: permute t2 to CCW
 * and call worker function to test for intersection.
 *
 * q1 and r1 both lie in the negative half-space defined by t2; p1 lies in
 * t2's plane or in its positive half-space.
 * The sign of dp2, dq2, and dr2 indicates whether the associated vertex
 * of t2 lies in the positive or negative half-space defined by t1.
 *****************************************************************************
 */
inline bool intersectOnePermutedTriangle(
  const Point3 &p1, const Point3 &q1, const Point3 &r1,
  const Point3 &p2, const Point3 &q2, const Point3 &r2,
  double dp2, double dq2, double dr2,  Vector3 &normal)
{
  /* Step 4: repeat Step 3, except doing it for triangle 2
     instead of triangle 1 */
  if (isGt(dp2, 0.0)) {
    if (isGt(dq2, 0.0)) {
      return intersectTwoPermutedTriangles(p1,r1,q1,r2,p2,q2);
    }
    else if (isGt(dr2, 0.0)) {
      return intersectTwoPermutedTriangles(p1,r1,q1,q2,r2,p2);
    }
    else{
      return intersectTwoPermutedTriangles(p1,q1,r1,p2,q2,r2);
    }
  }
  else if (isLt(dp2,  0.0)) {
    if (isLt(dq2, 0.0)) {
      return intersectTwoPermutedTriangles(p1,q1,r1,r2,p2,q2);
    }
    else if (isLt(dr2, 0.0)) {
      return intersectTwoPermutedTriangles(p1,q1,r1,q2,r2,p2);
    }
    else{
      return intersectTwoPermutedTriangles(p1,r1,q1,p2,q2,r2);
    }
  }
  else {
    if (isLt(dq2, 0.0)) {
      if (isGeq(dr2, 0.0)) {
        return intersectTwoPermutedTriangles(p1,r1,q1,q2,r2,p2);
      }
      else{
        return intersectTwoPermutedTriangles(p1,q1,r1,p2,q2,r2);
      }
    }
    else if (isGt(dq2, 0.0)) {
      if (isGt(dr2, 0.0)) {
        return intersectTwoPermutedTriangles(p1,r1,q1,p2,q2,r2);
      }
      else{
        return intersectTwoPermutedTriangles(p1,q1,r1,q2,r2,p2);
      }
    }
    else {
      if (isGt(dr2, 0.0)) {
        return intersectTwoPermutedTriangles(p1,q1,r1,r2,p2,q2);
      }
      else if (isLt(dr2, 0.0)) {
        return intersectTwoPermutedTriangles(p1,r1,q1,r2,p2,q2);
      }
      else{
        return intersectCoplanar3DTriangles(p1,q1,r1,p2,q2,r2,normal);
      }
    }
  }
}

/*!
 *****************************************************************************
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
 *****************************************************************************
 */
inline bool intersectTwoPermutedTriangles(const Point3 p1,
                                          const Point3 q1,
                                          const Point3 r1,
                                          const Point3 p2,
                                          const Point3 q2,
                                          const Point3 r2)
{
  /* Step 5: From step's 1 through 4, we now have two triangles that,
     if intersecting, have a line that intersects segments p1r1, p1q1,
     p2q2, and p2r2.  We check if these two intervals overlap:
   */

  if (isGt(Vector3(q1, q2).dot(Triangle3(q1, p2, p1).normal()), 0.0)) {
    return false;
  }
  if (isGt((Vector3(p1, r2).dot(Triangle3(p1, p2, r1).normal())), 0.0)) {
    return false;
  }

  return true;
}

/*!
 *****************************************************************************
 * Project (nearly) coplanar triangles 1 and 2 on an axis; call 2D worker
 * function to test for intersection.
 *****************************************************************************
 */
inline bool intersectCoplanar3DTriangles(const Point3& p1,
                                         const Point3& q1,
                                         const Point3& r1,
                                         const Point3& p2,
                                         const Point3& q2,
                                         const Point3& r2,
                                         Vector3 normal)
{
  /* Co-planar triangles are projected onto the axis that maximizes their
     area and the 2d intersection used to check if they intersect.
   */

  //find triangle with maximum area:
  for (int i=0; i<3; i++) {
    normal[i] = std::abs(normal[i]);
  }

  if ((isGt(normal[0], normal[2])) && (isGeq(normal[0], normal[1]))) {
    //if x projection area greatest, project on YZ and return 2D checker

    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[2],q1[1]),
                                       Point2::make_point( p1[2],p1[1]),
                                       Point2::make_point( r1[2],r1[1]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[2],q2[1]),
                                       Point2::make_point( p2[2],p2[1]),
                                       Point2::make_point( r2[2],r2[1]));

    return TriangleIntersection2D(t1_2da, t2_2da);
  }
  else if (isGt(normal[1],normal[2]) && isGeq(normal[1],normal[0])) {
    //if y projection area greatest, project on XZ and return 2D checker
    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[0],q1[2]),
                                       Point2::make_point( p1[0],p1[2]),
                                       Point2::make_point( r1[0],r1[2]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[0],q2[2]),
                                       Point2::make_point( p2[0],p2[2]),
                                       Point2::make_point( r2[0],r2[2]));

    return TriangleIntersection2D(t1_2da, t2_2da);
  }
  else{
    //if z projection area greatest, project on XY and return 2D checker
    const Triangle2 t1_2da = Triangle2(Point2::make_point(p1[0],p1[1]),
                                       Point2::make_point( q1[0],q1[1]),
                                       Point2::make_point( r1[0],r1[1]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(p2[0],p2[1]),
                                       Point2::make_point( q2[0],q2[1]),
                                       Point2::make_point( r2[0],r2[1]));

    return TriangleIntersection2D(t1_2da, t2_2da);
  }
  return false;
}

/** @} */

/** @{ @name 2D triangle-triangle intersection */
/*
 *******************************************************************************
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *
 * Note that edge and vertex intersections are reported as hits: an edge or
 * vertex intersecting any part of another triangle causes a return value of
 * true.  This is consistent with the paper.
 ********************************************************************************
 */
template < typename T >
bool intersect_tri2D_tri2D( const primal::Triangle< T, 2 >& t1,
                            const primal::Triangle< T, 2 >& t2)
{
  SLIC_CHECK_MSG(
    !t1.degenerate(),
    "\n\n WARNING \n\n Triangle " << t1 <<" is degenerate");
  SLIC_CHECK_MSG(
    !t2.degenerate(),
    "\n\n WARNING \n\n Triangle " << t2 <<" is degenerate");

  return TriangleIntersection2D(t1, t2);
}

/*!
 *****************************************************************************
 * \brief Orients triangle vertices so both triangles are CCW; calls worker.
 * \return status true iff t1 and t2 intersect
 *
 * Determine triangle orientation, then call the worker function with
 * vertices from t1 and t2 permuted to ensure CCW orientation.
 *****************************************************************************
 */
inline bool TriangleIntersection2D(const Triangle2& t1,
                                   const Triangle2& t2)
{
  if (isLt(twoDcross(t1[0],t1[1],t1[2]),0.0)) {
    if ((isLt(twoDcross(t2[0], t2[1], t2[2]),0.0))) {
      return intersectPermuted2DTriangles(t1[0], t1[2], t1[1],
                                          t2[0], t2[2], t2[1]);
    }
    else{return intersectPermuted2DTriangles(t1[0], t1[2], t1[1],
                                             t2[0], t2[1], t2[2]);
    }
  }
  else {
    if (isLt(twoDcross(t2[0], t2[1], t2[2]),0.0)) {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2],
                                          t2[0], t2[2], t2[1]);
    }
    else {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2],
                                          t2[0], t2[1], t2[2]);
    }
  }
}

/*!
 *****************************************************************************
 * This function finds where p1 lies in relation to the vertices of t2
 * and calls either checkEdge() or checkVertex().
 * \return status true iff triangle p1 q1 r1 intersects triangle p2 q2 r2.
 *****************************************************************************
 */
inline bool intersectPermuted2DTriangles(const Point2& p1,
                                         const Point2& q1,
                                         const Point2& r1,
                                         const Point2& p2,
                                         const Point2& q2,
                                         const Point2& r2)
{
  // Step 2: Orient triangle 2 to be counter clockwise and break the problem
  // into two generic cases (where we test the vertex for intersection or the
  // edges).
  //
  // See paper at https://hal.inria.fr/inria-00072100/document for more details

  if (isGeq(twoDcross(p2,q2,p1), 0.0 )) {
    if (isGeq(twoDcross(q2,r2,p1), 0.0 )) {
      if (isGeq(twoDcross(r2,p2,p1), 0.0)) {
        return true;
      }
      else{return checkEdge(p1,q1,r1,p2,r2); //T1 clockwise
      }
    }
    else {
      if (isGeq(twoDcross(r2,p2,p1), 0.0)) {
        //5 region decomposition with p1 in the +-- region
        return checkEdge(p1,q1,r1,r2,q2);
      }
      else{return checkVertex(p1,q1,r1,p2,q2,r2);
      }
    }
  }
  else {
    if (isGeq(twoDcross(q2,r2,p1), 0.0)) {
      if (isGeq(twoDcross(r2,p2,p1), 0.0)) {
        //four region decomposition.  ++- region
        return checkEdge(p1,q1,r1,q2,p2);
      }
      else{return checkVertex(p1,q1,r1,q2,r2,p2);
      }
    }
    else{return checkVertex(p1,q1,r1,r2,p2,q2);
    }
  }
}

/*!
 *****************************************************************************
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2p2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 *****************************************************************************
 */
inline bool checkEdge(const Point2 p1,
                      const Point2 q1,
                      const Point2 r1,
                      const Point2 p2,
                      const Point2 r2)
{
  if (isGeq(twoDcross(r2,p2,q1),0.0)) {
    if (isGeq(twoDcross(p1,p2,q1), 0.0)) {
      if (isGeq(twoDcross(p1,q1,r2), 0.0)) {
        return true;
      }
      else{
        return false;
      }
    }
    else {
      if (isGeq(twoDcross(q1,r1,p2), 0.0)) {
        if (isGeq(twoDcross(r1,p1,p2), 0.0)) {
          return true;
        }
        else{
          return false;
        }
      }
      else{
        return false;
      }
    }
  }
  else {
    if (isGeq(twoDcross(r2,p2,r1), 0.0)) {
      if (isGeq(twoDcross(p1,p2,r1), 0.0)) {
        if (isGeq(twoDcross(q1,r1,r2), 0.0)) {
          return true;
        }
        else{
          return false;
        }
      }
      else{
        return false;
      }
    }
    else{
      return false;
    }
  }
}

/*!
 *****************************************************************************
 * \brief Check for 2D triangle-edge intersection, given p1 close to r2.
 * \return status true iff coplanar CCW triangles 1 and 2 intersect.
 *****************************************************************************
 */
inline bool checkVertex(const Point2 p1,
                        const Point2 q1,
                        const Point2 r1,
                        const Point2 p2,
                        const Point2 q2,
                        const Point2 r2)
{
  if (isGeq(twoDcross(r2,p2,q1),0.0)) {
    if (isGeq(twoDcross(q2,r2,q1),0.0)) {
      if (isGeq(twoDcross(p1,p2,q1),0.0)) {
        if (isLeq(twoDcross(p1,q2,q1),0.0)) {
          return true;
        }
        else{
          return false;
        }
      }
      else {
        if (isGeq(twoDcross(p1,p2,r1),0.0)) {
          if (isGeq(twoDcross(r2,p2,r1),0.0)) {
            return true;
          }
          else{
            return false;
          }
        }
        else{
          return false;
        }
      }
    }
    else {
      if (isLeq(twoDcross(p1,q2,q1),0.0)) {
        if (isGeq(twoDcross(q2,r2,r1),0.0)) {
          if (isGeq(twoDcross(q1,r1,q2),0.0)) {
            return true;
          }
          else{
            return false;
          }
        }
        else{
          return false;
        }
      }
      else{
        return false;
      }
    }
  }
  else {
    if (isGeq(twoDcross(r2,p2,r1),0.0)) {
      if (isGeq(twoDcross(q1,r1,r2),0.0)) {
        if (isGeq(twoDcross(r1,p1,p2),0.0)) {
          return true;
        }
        else{
          return false;
        }
      }
      else {
        if (isGeq(twoDcross(q1,r1,q2),0.0)) {
          if (isGeq(twoDcross(q2,r2,r1),0.0)) {
            return true;
          }
          else{
            return false;
          }
        }
        else{
          return false;
        }
      }
    }
    else{
      return false;
    }
  }
}

/*!
 *****************************************************************************
 * \brief Compute cross product of two 2D vectors as if they were 3D.
 * \return Cross product of A C and B C.
 *
 * This function treats three Point2 values as corners of a 3D triangle with
 * zero Z-coordinate.  Thus we can calculate the cross product of A C with
 * B C using only the k-hat term, since the other terms go to zero.  A
 * positive value indicates CCW orientation.
 *****************************************************************************
 */
inline double twoDcross(const Point2& A, const Point2& B, const Point2& C)
{
  return  (((A[0]-C[0])*(B[1]-C[1])-(A[1]-C[1])*(B[0]-C[0])));
}

/*!
 *****************************************************************************
 * \brief Checks if x > y, within a specified tolerance.
 *****************************************************************************
 */
inline bool isGt(double x, double y, double EPS)
{
  return ((x > y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/*!
 *****************************************************************************
 * \brief Checks if x < y, within a specified tolerance.
 *****************************************************************************
 */
inline bool isLt(double x, double y, double EPS)
{
  return ((x < y) && !(axom::utilities::isNearlyEqual(x, y, EPS)));
}

/*!
 *****************************************************************************
 * \brief Checks if x <= y, within a specified tolerance.
 *****************************************************************************
 */
inline bool isLeq(double x, double y, double EPS)
{
  return !(isGt(x,y,EPS));
}

/*!
 *****************************************************************************
 * \brief Checks if x >= y, within a specified tolerance.
 *****************************************************************************
 */
inline bool isGeq(double x, double y, double EPS)
{
  return !(isLt(x,y,EPS));
}

/*!
 *****************************************************************************
 * \brief Check if x, y, and z all have the same sign.
 ****************************************************************************
 */
inline bool signMatch(double x, double y, double z, double EPS)
{
  return ((isGt(x*y, 0.0, EPS)) && (isGt(x*z, 0.0, EPS)));
}

/** @} */

/*!
 *******************************************************************************
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *      ip returns the intersection point on S.
 * \return status true iff R intersects with S, otherwise, false.
 *******************************************************************************
 */
template < typename T >
bool intersect_ray_seg( const primal::Ray< T,2 >& R,
                        const primal::Segment< T,2 >& S,
                        primal::Point< T,2 >& ip )
{
  // STEP 0: Construct a ray from the segment, i.e., represent the
  // segment in parametric form S(t1)=A+td, t \in [0,1]
  Ray< T,2 > R2( S );

  // Step 1: Equating R(t0)=S(t1) yields a system of two equations and
  // two unknowns, namely, t0 and t1. We can solve this system directly
  // using Cramer's Rule.
  const double denom = primal::determinant(
    R.direction()[0], (-1.0)*R2.direction()[0],
    R.direction()[1], (-1.0)*R2.direction()[1]     );

  // STEP 2: if denom is zero, the system is singular, which implies that the
  // ray and the segment are parallel
  const double parepsilon = 1.0e-9;
  if ( axom::utilities::isNearlyEqual( denom, 0.0, parepsilon ) ) {

    // ray and segment are parallel
    return false;

  }

  // STEP 3: Solve for t0 and t1 directly using cramer's rule
  const double alpha = S.source()[0] - R.origin()[0];
  const double beta  = S.source()[1] - R.origin()[1];

  const double t0 = primal::determinant( alpha, (-1.0)*R2.direction()[0],
                                         beta,
                                         (-1.0)*R2.direction()[1] )/denom;

  const double t1 = primal::determinant( R.direction()[0], alpha,
                                         R.direction()[1], beta   )/denom;

  // STEP 4: Define lower/upper threshold
  const double tlow  = 0.0-1.0e-9;
  const double thigh = 1.0+1.0e-9;

  // STEP 5: Necessary and sufficient criteria for an intersection between
  // ray, R(t0),  and a finite segment S(t1) are:
  // 1. t0 >= tlow w.r.t. the ray R(t0).
  // 2. tlow >= t1 >= thigh w.r.t. the segment S(t1).
  if ( (t0 >= tlow) && (t1 >= tlow) && (t1 <= thigh) ) {
    ip = R2.at( t1 );
    return true;
  }

  // STEP 6: Ray does not intersect the segment
  return false;
}

/*!
 *******************************************************************************
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *      ip the point of intersection on R.
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * Computes Ray Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 *******************************************************************************
 */
template < typename T, int DIM >
bool intersect_ray_bbox(const primal::Ray< T,DIM > & R,
                        const primal::BoundingBox< T,DIM > & bb,
                        primal::Point< T,DIM > & ip)
{
  T tmin = std::numeric_limits< T >::min();
  SLIC_ASSERT(tmin>=0.0);
  T tmax = std::numeric_limits< T >::max();

  for (int i=0; i<DIM; i++) {
    if (axom::utilities::isNearlyEqual(R.direction()[i],
                                       std::numeric_limits< T >::min(),
                                       1.0e-9 )) {
      T pointDim =  R.origin()[i];
      if ((pointDim<bb.getMin()[i]) || (pointDim>bb.getMax()[i])) {
        return false;
      }
    }
    else{
      T ood = (static_cast< T >(1.0)) / (R.direction()[i]);
      T t1 = ((bb.getMin()[i]- R.origin()[i])*ood);
      T t2 = ((bb.getMax()[i]- R.origin()[i])*ood);

      if (t1>t2) {
        std::swap(t1,t2);
      }

      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      if (tmin > tmax) {
        return false;
      }
    }
  }

  for (int i = 0; i < DIM; i++) {
    ip.data()[i] = R.origin()[i] + R.direction()[i] * tmin;
  }

  return true;
}

/*!
 *******************************************************************************
 * \brief Computes the intersection of the given segment, S, with the Box, bb.
 *     ip the point of intersection on S.
 * \return status true iff bb intersects with S, otherwise, false.
 *
 * Computes Segment Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 * WIP: More test cases for this
 *******************************************************************************
 */
template < typename T, int DIM >
bool intersect_seg_bbox( const primal::Segment< T,DIM > & S,
                         const primal::BoundingBox< T,DIM > & bb,
                         primal::Point< T,DIM > & ip)
{
  T tmin = std::numeric_limits< T >::min();
  primal::Vector< T,DIM > direction(S.source(), S.target());
  T tmax = direction.norm();
  primal::Ray< T,DIM > R(S.source(), direction);

  // These operations constrain the parameter specifying ray-slab intersection
  // points to exclude points not within the segment.
  tmin = static_cast< T >(0);
  tmax = static_cast< T >(1);

  for (int i=0; i<DIM; i++) {
    if (axom::utilities::isNearlyEqual(R.direction()[i],
                                       std::numeric_limits< T >::min(),
                                       1.0e-9 )) {
      T pointDim =  R.origin()[i];
      if ((pointDim<bb.getMin()[i]) || (pointDim>bb.getMax()[i])) {
        return false;
      }
    }
    else{
      T ood = (static_cast< T >(1.0)) / (R.direction()[i]);
      T t1 = ((bb.getMin()[i]- R.origin()[i])*ood);
      T t2 = ((bb.getMax()[i]- R.origin()[i])*ood);

      if (t1>t2) {
        std::swap(t1,t2);
      }

      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      if (tmin > tmax) {
        return false;
      }
    }
  }

  for (int i=0; i< DIM; i++) {
    ip.data()[i] = R.origin()[i] + R.direction()[i]*tmin;
  }

  return true;
}

typedef primal::Vector< double, 3 > Vector3;
bool intervalsDisjoint(double d0, double d1, double d2, double r);
bool crossEdgesDisjoint(double d0, double d1, double r);

/** @{ @name Triangle-bbox intersection */

/*!
 *******************************************************************************
 * \brief Determines if a triangle and a bounding box intersect
 *        (but does not find the intersections)
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 *******************************************************************************
 */
template < typename T >
bool intersect_tri_bbox( const primal::Triangle< T, 3 >& tri,
                         const primal::BoundingBox< T, 3 >& bb)
{
  // Note: Algorithm is derived from the one presented in chapter 5.2.9 of
  //   Real Time Collision Detection book by Christer Ericson
  // based on Akenine-Moller algorithm (Journal of Graphics Tools)
  //
  // It uses the Separating Axis Theorem to look for disjoint projections
  // along various axes associated with Faces and Edges of the AABB and triangle.
  // There are 9 tests for the cross products of edges
  //           3 tests for the AABB face normals
  //           1 test for the triangle face normal
  // We use early termination if we find a separating axis between the shapes

  typedef typename BoundingBox< T,3 >::PointType PointType;
  typedef typename BoundingBox< T,3 >::VectorType VectorType;

  // Extent: vector center to max corner of BB
  VectorType e = 0.5 * bb.range();

  // Make the AABB center the origin by moving the triangle vertices
  PointType center(bb.getMin().array() + e.array());
  VectorType v[3] = { VectorType(center, tri[0]),
                      VectorType(center, tri[1]),
                      VectorType(center, tri[2]) };

  // Create the edge vectors of the triangle
  VectorType f[3] = { v[1] - v[0], v[2] - v[1],  v[0] - v[2] };

/* *INDENT-OFF* */

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
/* *INDENT-ON* */

#undef XEDGE_R
#undef XEDEG_S

  /// Test face normals of bounding box (3 tests)
  if ( intervalsDisjoint(v[0][0], v[1][0], v[2][0], e[0]) ||
       intervalsDisjoint(v[0][1], v[1][1], v[2][1], e[1]) ||
       intervalsDisjoint(v[0][2], v[1][2], v[2][2], e[2]) ) {
    return false;
  }

  /// Final test -- face normal of triangle's plane
  VectorType planeNormal  = VectorType::cross_product(f[0],f[1]);
  double planeDist    = planeNormal.dot(tri[0]);

  double r = e[0]* std::abs( planeNormal[0]) +
             e[1]* std::abs( planeNormal[1]) +
             e[2]* std::abs( planeNormal[2]);
  double s = planeNormal.dot(center) - planeDist;

  return std::abs(s) <= r;
}

// ------------------------------------------------------------------------------

/**
 * \brief Helper function to find disjoint projections for the AABB-triangle test
 * \param d0 The first value defining the test interval
 * \param d1 The second value defining the test interval
 * \param d2 The third value defining the test interval
 * \param r Radius of projection
 * \return True of the intervals are disjoint, false otherwise
 */
bool intervalsDisjoint(double d0, double d1, double d2, double r)
{
  if (d1 < d0) {
    std::swap(d1,d0);  // d0 < d1
  }
  if (d2 > d1) {
    std::swap(d2,d1);  // d1 is max(d0,d1,d2)
  }
  else if (d2 < d0) {
    std::swap(d2,d0);  // d0 is min(d0,d1,d2)

  }
  SLIC_ASSERT(  d0 <= d1 && d0 <= d2);
  SLIC_ASSERT(  d1 >= d0 && d1 >= d2);

  return d1 < -r || d0 > r;
}

/**
 * \brief Helper function for Triangle/BoundingBox intersection test
 */
bool crossEdgesDisjoint(double d0, double d1, double r)
{
  return std::max( -std::max(d0,d1), std::min(d0,d1) ) > r;
}

/** @} */

/** @{ @name Triangle-ray intersection */

/*!
 *******************************************************************************
 * \brief Tests if 3D triangle tri intersects with 3D ray R.
 * \return status true iff tri intersects with R, otherwise, false.
 *
 * This algorithm is modeled after Woop, Benthin, and Wald (2013).  It
 * transforms the ray onto the unit z vector and applies the same transform
 * to the triangle's vertices.  This transform simplifies the calculation of
 * barycentric coordinates of the positive z-axis's intersection with the
 * triangle.  If any of these coordinates are less than zero, the ray misses.
 *
 * If any are equal to zero, more care is needed to check if the ray hits
 * the edge or misses.  Additionally, the code falls back to double precision
 * for cases that call for it.
 *
 * Sven Woop, Carsten Benthin, Ingo Wald, "Watertight Ray/Triangle
 * Intersection," Journal of Computer Graphics Techniques (JCGT), vol. 2,
 * no. 1, 65–82, 2013  http://jcgt.org/published/0002/01/05/
 *******************************************************************************
 */
template < typename T >
bool intersect_tri_ray(const Triangle< T, 3 >& tri, const Ray< T,3 >& R)
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

  //find out dimension where ray direction is maximal
  int kx,ky,kz;
  float rX,rY,rZ;

  //assign initial directions
  rX = R.direction()[0];
  rY = R.direction()[1];
  rZ = R.direction()[2];

  //make all dimensions positive
  if (rX < 0.0f) {
    rX *= -1.0f;
  }
  if (rY < 0.0f) {
    rY *= -1.0f;
  }
  if (rZ < 0.0f) {
    rZ *= -1.0f;
  }

  //z-direction largest
  if ((rZ>=rX) && (rZ >= rY)) {
    kz=2;
  }
  //y direction largest
  else if ((rY >= rX) && (rY >= rZ)) {
    kz=1;
  }
  //x direction largest
  else{
    kz=0;
  }

  //assign other dimensions of the ray
  kx = (kz+1) % 3;
  ky = (kz+2) % 3;

  //if necessary swap  ky and kx to preserve triangle winding
  if (R.direction()[kz] < 0.0f) {
    std::swap(kx, ky);
  }

  //calculate shear constants
  float Sx = R.direction()[kx]/R.direction()[kz];
  float Sy = R.direction()[ky]/R.direction()[kz];
  float Sz = 1.0f/R.direction()[kz];

  //A,B,C are the triangle vertices
  float A[3],B[3],C[3];

  //calculate vertices relative to the ray origin
  A[kx] = tri[0].array()[kx] - R.origin()[kx];
  A[ky] = tri[0].array()[ky] - R.origin()[ky];
  A[kz] = tri[0].array()[kz] - R.origin()[kz];

  B[kx] = tri[1].array()[kx] - R.origin()[kx];
  B[ky] = tri[1].array()[ky] - R.origin()[ky];
  B[kz] = tri[1].array()[kz] - R.origin()[kz];

  C[kx] = tri[2].array()[kx] - R.origin()[kx];
  C[ky] = tri[2].array()[ky] - R.origin()[ky];
  C[kz] = tri[2].array()[kz] - R.origin()[kz];

  //shear and scale the vertices
  const float Ax = A[kx] - Sx*A[kz];
  const float Ay = A[ky] - Sy*A[kz];
  const float Bx = B[kx] - Sx*B[kz];
  const float By = B[ky] - Sy*B[kz];
  const float Cx = C[kx] - Sx*C[kz];
  const float Cy = C[ky] - Sy*C[kz];

  //scaled barycentric coordinates
  float U = Cx*By - Cy*Bx;
  float V = Ax*Cy - Ay*Cx;
  float W = Bx*Ay - By*Ax;

  //fallback to test against edges using double precision
  if (U == 0.0f || V == 0.0f || W == 0.0f) {
    double CxBy = (double)Cx*(double)By;
    double CyBx = (double)Cy*(double)Bx;
    U = (float)(CxBy - CyBx);

    double AxCy = (double)Ax*(double)Cy;
    double AyCx = (double)Ay*(double)Cx;
    V = (float)(AxCy - AyCx);

    double BxAy = (double)Bx*(double)Ay;
    double ByAx = (double)By*(double)Ax;
    W = (float)(BxAy - ByAx);
  }

  //edge testing
  if ((U<0.0f || V<0.0f || W<0.0f) && (U>0.0f || V>0.0f || W>0.0f)) {
    return false;
  }

  //calculate determinant
  float det = U + V + W;

  if (det == 0.0f) {
    return false;
  }

  //calculate scaled z-coordinates of the vertices and use them to calculate hit
  // distance
  const float Az = Sz*A[kz];
  const float Bz = Sz*B[kz];
  const float Cz = Sz*C[kz];
  const float Q = U*Az + V*Bz +W*Cz;

  //make sure hit is in correct direction
  bool dir = ((Q<0.0f) && !(det<0.0f)) || ((det < 0.0f) && !(Q<0.0f));
  if (dir) {
    return false;
  }
  return true;
}

/** @} */

template < typename T >
bool intersect_tri_segment(const Triangle< T, 3 >& tri, const Segment< T,3 >& S)
{
  typedef Vector< T,3 > Vector3;
  Ray< T,3 > r1(S.source(), Vector3(S.source(), S.target()));
  Ray< T,3 > r2(S.target(), Vector3(S.target(), S.source()));

  //Ray-triangle intersection does not check endpoints, so we explicitly check here
  if (tri.checkInTriangle(r1.origin())) {
    return true;
  }
  if (tri.checkInTriangle(r2.origin())) {
    return true;
  }

  //if the two rays formed by the endpoints of the segment intersect the triangle,
  //then the triangle must intersect the triangle
  if ((intersect_tri_ray(tri, r1)) && (intersect_tri_ray(tri, r2))) {
    return true;
  }
  return false;
}

} /* end namespace detail */
} /* end namespace primal */
} /* end namespace axom */

#endif /* INTERSECTION_IMPL_HPP_ */
