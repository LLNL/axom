/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*!
 *******************************************************************************
 * \file Intersection.hpp
 *
 * \date Jan 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef INTERSECTION_HPP_
#define INTERSECTION_HPP_

#include "quest/BoundingBox.hpp"
#include "quest/Determinants.hpp"
#include "quest/Point.hpp"
#include "quest/Ray.hpp"
#include "quest/Segment.hpp"
#include "quest/Triangle.hpp"

#include "common/Utilities.hpp"

namespace quest {

/*!
 *******************************************************************************
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 * \param [in] R user-supplied ray R.
 * \param [in] S user-supplied segment S.
 * \param [in,out] ip the point of intersection.
 * \return status true iff R intersects with S, otherwise, false.
 *******************************************************************************
 */
template < typename T >
bool intersect( const Ray<T,2>& R, const Segment<T,2>& S, Point<T,2>& ip )
{
   // STEP 0: Construct a ray from the segment, i.e., represent the
   // segment in parametric form S(t1)=A+td, t \in [0,1]
   Ray<T,2> R2( S );

   // Step 1: Equating R(t0)=S(t1) yields a system of two equations and
   // two unknowns, namely, t0 and t1. We can solve this system directly
   // using Cramer's Rule.
   const double denom = math::determinant(
                               R.direction()[0], (-1.0)*R2.direction()[0],
                               R.direction()[1], (-1.0)*R2.direction()[1]     );


   // STEP 2: if denom is zero, the system is singular, which implies that the
   // ray and the segment are parallel
   if ( asctoolkit::utilities::isNearlyEqual( denom, 0.0, 1.0e-9 ) ) {

       // ray and segment are parallel
       return false;

   }

   // STEP 3: Solve for t0 and t1 directly using cramer's rule
   const double alpha = S.source()[0] - R.origin()[0];
   const double beta  = S.source()[1] - R.origin()[1];

   const double t0 = math::determinant(alpha, (-1.0)*R2.direction()[0],
                                       beta,  (-1.0)*R2.direction()[1] )/denom;

   const double t1 = math::determinant( R.direction()[0], alpha,
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
 * \param [in] R user-supplied ray R.
 * \param [in] bb user-supplied box bb.
 * \param [out] ip the point of intersection.
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * Computes Ray Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 *******************************************************************************
 */
template < typename T, int DIM>
bool intersect (const Ray<T,DIM> & R,
		const BoundingBox<T,DIM> & bb,
		Point<T,DIM> & ip)
{
  T tmin = std::numeric_limits<T>::min();
  SLIC_ASSERT(tmin>=0.0);
  T tmax = std::numeric_limits<T>::max();

  for (int i=0; i<DIM; i++)
  {
    if (asctoolkit::utilities::isNearlyEqual(R.direction()[i],
					     std::numeric_limits<T>::min(),
					     1.0e-9 ))
    {
      T pointDim =  R.origin()[i];
      if ((pointDim<bb.getMin()[i]) || (pointDim>bb.getMax()[i]))
      {
	return false;
      }
    }
    else
    {
      T ood = (static_cast<T>(1.0)) / (R.direction()[i]);
      T t1 = ((bb.getMin()[i]- R.origin()[i])*ood);
      T t2 = ((bb.getMax()[i]- R.origin()[i])*ood);

      if (t1>t2)
      {
	std::swap(t1,t2);
      }

      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      if (tmin > tmax)
      {
	return false;
      }
    }
  }

  for (int i = 0; i < DIM; i++)
  {
    ip.data()[i] = R.origin()[i] + R.direction()[i] * tmin;
  }

  return true;
}


/*!
 *******************************************************************************
 * \brief Computes the intersection of the given segment, S, with the Box, bb.
 * \param [in] S user-supplied segment S.
 * \param [in] bb user-supplied box bb.
 * \param [out] ip the point of intersection.
 * \return status true iff bb intersects with S, otherwise, false.
 *
 * Computes Segment Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 * WIP: More test cases for this
 *******************************************************************************
 */
template < typename T, int DIM>
bool intersect (const Segment<T,DIM> & S,
		const BoundingBox<T,DIM> & bb,
		Point<T,DIM> & ip)
{
  T tmin = std::numeric_limits<T>::min();
  Vector<T,DIM> direction(S.source(), S.target());
  T tmax = direction.norm();
  Ray<T,DIM> R(S.source(), direction);

  // These operations constrain the parameter specifying ray-slab intersection
  // points to exclude points not within the segment.
  tmin = static_cast<T>(0);
  tmax = static_cast<T>(1);

  for (int i=0; i<DIM; i++)
  {
    if (asctoolkit::utilities::isNearlyEqual(R.direction()[i],
					     std::numeric_limits<T>::min(),
					     1.0e-9 ))
    {
      T pointDim =  R.origin()[i];
      if ((pointDim<bb.getMin()[i]) || (pointDim>bb.getMax()[i]))
      {
	return false;
      }
    }
    else
    {
      T ood = (static_cast<T>(1.0)) / (R.direction()[i]);
      T t1 = ((bb.getMin()[i]- R.origin()[i])*ood);
      T t2 = ((bb.getMax()[i]- R.origin()[i])*ood);

      if (t1>t2)
      {
	std::swap(t1,t2);
      }

      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      if (tmin > tmax)
      {
	return false;
      }
    }
  }

  for (int i=0; i< DIM; i++)
  {
    ip.data()[i] = R.origin()[i] + R.direction()[i]*tmin;
  }

  return true;
}


namespace {

  typedef quest::Vector<double, 3> Vector3;
  typedef quest::Point<double, 3> Point3;
  typedef quest::Triangle<double, 3> Triangle3;
  typedef quest::Triangle<double, 2> Triangle2;
  typedef quest::Point<double, 2> Point2;

  /*!
   *****************************************************************************
   * \brief Checks if x > y
   * \param [in] x first value
   * \param [in] y value to compare to x
   * \return status true iff x - 1.0e-12 > y
   *
   * Helper function for T-T intersect.  > paired with is_nearly_equal
   *****************************************************************************
   */
  inline bool isGt(double x, double y, double EPS=1.0e-12)
  {
    return ((x > y) && !(asctoolkit::utilities::isNearlyEqual(x, y, EPS)));
  }

  /*!
   *****************************************************************************
   * \brief Checks if x < y
   * \param [in] x first value
   * \param [in] y value to compare to x
   * \return status true iff x + 1.0e-12 < y
   *
   * Helper function for T-T intersect.  < paired with is_nearly_equal
   *****************************************************************************
   */
  inline bool isLt(double x, double y, double EPS=1.0e-12)
  {
    return ((x < y) && !(asctoolkit::utilities::isNearlyEqual(x, y, EPS)));
  }

  /*!
   *****************************************************************************
   * \brief Checks if x <= y
   * \param [in] x first value
   * \param [in] y value to compare to x
   * \return status true iff x - 1.0e-12 <= y
   *
   * Helper function for T-T intersect.
   *****************************************************************************
   */
  inline bool isLeq(double x, double y)
  {
    return !(isGt(x,y));
  }

  /*!
   *****************************************************************************
   * \brief Checks if x >= y
   * \param [in] x first value
   * \param [in] y value to compare to x
   * \return status true iff x + 1.0e-12 >= y
   *
   * Helper function for T-T intersect.  < paired with is_nearly_equal
   *****************************************************************************
   */
  inline bool isGeq(double x, double y)
  {
    return !(isLt(x,y));
  }

  /*!
   *****************************************************************************
    * \brief Check sign match
    * \param [in] x first test value
    * \param [in] y second test value
    * \param [in] z third test value
    * \return status true iff x and y have same sign, x and z have same sign
    *
    * Helper function for T-T intersect.
    ****************************************************************************
    */
  inline bool signMatch(double x, double y, double z)
  {
    return ((isGt(x*y, 0.0)) &&  (isGt(x*z, 0.0)));
  }

  /*!
   *****************************************************************************
   * \brief Check 2D triangle orientation.
   * \param [in] A first corner of triangle
   * \param [in] B second corner of triangle
   * \param [in] C third corner of triangle
   * \return Cross product of A C and B C.
   *
   * This function treats three Point2 values as corners of a 3D triangle with
   * zero Z-coordinate.  Thus we can calculate the cross product of A C with
   * B C using only the k-hat term, since the other terms go to zero.  A
   * positive value indicates CCW orientation.
   *
   * Helper function for T-T intersect.
   *****************************************************************************
   */
  inline double checkCCW(const Point2& A, const Point2& B, const Point2& C)
  {
    return  (((A[0]-C[0])*(B[1]-C[1])-(A[1]-C[1])*(B[0]-C[0])));
  }


  /*!
   *****************************************************************************
   * \brief Check for 2D triangle-edge intersection.
   * \param [in] Ai first corner of triangle
   * \param [in] Bi second corner of triangle
   * \param [in] Ci third corner of triangle
   * \param [in] Aj first point of segment
   * \param [in] Bj second corner of segment
   * \return status true iff the edge  Aj Cj intersects the triangle Ai Bi Ci.
   *
   * Helper function for T-T intersect.
   *****************************************************************************
   */
  inline bool checkEdge(const Point2 Ai,
			const Point2 Bi,
			const Point2 Ci,
			const Point2 Aj,
			const Point2 Cj)
  {
    if (isGeq(checkCCW(Cj,Aj,Bi),0.0)) {
      if (isGeq(checkCCW(Ai,Aj,Bi), 0.0)) {
	if (isGeq(checkCCW(Ai,Bi,Cj), 0.0))
	  return true;
	else
	  return false;
      }
      else {
	if (isGeq(checkCCW(Bi,Ci,Aj), 0.0)) {
	  if (isGeq(checkCCW(Ci,Ai,Aj), 0.0))
	    return true;
	  else
	    return false;
	}
	else
	  return false;
      }
    }
    else {
      if (isGeq(checkCCW(Cj,Aj,Ci), 0.0)) {
	if (isGeq(checkCCW(Ai,Aj,Ci), 0.0)) {
	  if (isGeq(checkCCW(Ai,Ci,Cj), 0.0))
	    return true;
	  else {
	    if (isGeq(checkCCW(Bi,Ci,Cj), 0.0))
	      return true;
	    else
	      return false;
	  }
	}
	else
	  return false;
      }
      else
	return false;
    }
  }

  /* Note that Evan's original code had checkVertex() with exactly the same
     text as checkEdge().  This is almost certainly not what he intended, since
     the two functions are called from CheckReorientedPoints2D() in different
     code paths.

     2016-09-16 Until I (Arlie Capps) understand this code better, I'm including
     checkVertex() as a stub that calls checkEdge() in order to let the code
     compile.
  */
  inline bool checkVertex(const Point2 Ai,
			const Point2 Bi,
			const Point2 Ci,
			const Point2 Aj,
			const Point2 Cj)
  {
    // TODO:  FIXME:  Properly implement the body of this function.

    if (isGeq(checkCCW(Cj,Aj,Bi),0.0)) {
      if (isGeq(checkCCW(Ai,Aj,Bi), 0.0)) {
	if (isGeq(checkCCW(Ai,Bi,Cj), 0.0))
	  return true;
	else
	  return false;
      }
      else {
	if (isGeq(checkCCW(Bi,Ci,Aj), 0.0)) {
	  if (isGeq(checkCCW(Ci,Ai,Aj), 0.0))
	    return true;
	  else
	    return false;
	}
	else
	  return false;
      }
    }
    else {
      if (isGeq(checkCCW(Cj,Aj,Ci), 0.0)) {
	if (isGeq(checkCCW(Ai,Aj,Ci), 0.0)) {
	  if (isGeq(checkCCW(Ai,Ci,Cj), 0.0))
	    return true;
	  else {
	    if (isGeq(checkCCW(Bi,Ci,Cj), 0.0))
	      return true;
	    else
	      return false;
	  }
	}
	else
	  return false;
      }
      else
	return false;
    }
  }

  /*!
   *****************************************************************************
   * \brief Runs tri-tri intersection tests, allowing for corner permutation.
   * \param [in] Ai first corner of triangle 1
   * \param [in] Bi second corner of triangle 1
   * \param [in] Ci third corner of triangle 1
   * \param [in] Aj first corner of triangle 2
   * \param [in] Bj second corner of triangle 2
   * \param [in] Cj third corner of triangle 2
   * \return status true iff triangle Ai Bi Ci intersects triangle Aj Bj Cj.
   *
   * Helper function for T-T intersect.
   *****************************************************************************
   */
  inline bool CheckReorientedPoints2D(const Point2& Ai,
				      const Point2& Bi,
				      const Point2& Ci,
				      const Point2& Aj,
				      const Point2& Bj,
				      const Point2& Cj)
  {
    // Step 2: Orient triangle 2 to be counter clockwise and break the problem
    // into two generic cases (where we test the vertex for intersection or the
    // edges).
    //
    // See paper at https://hal.inria.fr/inria-00072100/document for more details

    if (isGeq(checkCCW(Aj,Bj,Ai), 0.0 )) {
      if (isGeq(checkCCW(Bj,Cj,Ai), 0.0 )) {
	if (isGeq(checkCCW(Cj,Aj,Ai), 0.0)) {
	  return true;
	}
	else return checkEdge(Ai,Bi,Ci,Aj,Cj); //T1 clockwise
      }
      else {
	if (isGeq(checkCCW(Cj,Aj,Ai), 0.0)){
	  //5 region decomposistion with Ai in the +-- region
	  return checkEdge(Ai,Bi,Ci,Cj,Bj);
	}
	else return checkVertex(Ai,Bi,Ci,Aj,Cj);
      }
    }
    else {
      if (isGeq(checkCCW(Bj,Cj,Ai), 0.0)) {
	if (isGeq(checkCCW(Cj,Aj,Ai), 0.0)) {
	  //four region decomposistion.  ++- region
	  return checkEdge(Ai,Bi,Ci,Bj,Aj);
	}
	else return checkVertex(Ai,Bi,Ci,Bj,Aj);
      }
      else return checkVertex(Ai,Bi,Ci,Cj,Bj);
    }
  }

  /*!
   *****************************************************************************
   * \brief Detect triangle intersection in 2D
   * \param [in] t1 First Triangle<T, 2>
   * \param [in] t2 Second Triangle<T, 2>
   * \return status true iff t1 and t2 intersect
   *
   * Determine triangle orientation, then call the worker function with
   * vertices from t1 and t2 permuted to ensure CCW orientation.
   *
   * Helper function that does T-T intersect for 2d
   *****************************************************************************
   */
  inline bool TriangleIntersection2D(const Triangle2& t1,
				     const Triangle2& t2)
  {
    if (isLt(checkCCW(t1.A(),t1.B(),t1.C()),0.0)) {
      if ((isLt(checkCCW(t2.A(), t2.B(), t2.C()),0.0))) {
	return CheckReorientedPoints2D(t1.A(), t1.C(), t1.B(),
				       t2.A(), t2.C(), t2.B());
      }
      else return CheckReorientedPoints2D(t1.A(), t1.C(), t1.B(),
					  t2.A(), t2.B(), t2.C());
    }
    else {
      if (isLt(checkCCW(t2.A(), t2.B(), t2.C()),0.0)) {
	return CheckReorientedPoints2D(t1.A(), t1.B(), t1.C(),
				       t2.A(), t2.C(), t2.B());
      }
      else {
	return CheckReorientedPoints2D(t1.A(), t1.B(), t1.C(),
				       t2.A(), t2.B(), t2.C());
      }
    }
  }

  /*!
   *****************************************************************************
   * \brief Returns whether the two intervals overlap
   * \param [in] Ai first corner of triangle 1
   * \param [in] Bi second corner of triangle 1
   * \param [in] Ci third corner of triangle 1
   * \param [in] Aj first corner of triangle 2
   * \param [in] Bj second corner of triangle 2
   * \param [in] Cj third corner of triangle 2
   * \return status true iff triangles 1 and 2 intersect.
   *
   * This function implements Equation 1 from Devillers and Guigue (2002), p.8
   * using a hint from p.10 that greatly simplifies the computation.  The
   * triangles have been carefully rotated so Ai is across plane j from Bi and
   * Ci, and Aj is across plane i from Bj and Cj.  Previous tests have ruled
   * out cases where planes i and j are parallel or identical, as well as cases
   * where triangle 1 lies entirely off to one side of triangle 2 (and vice
   * versa).
   *
   * The core of the method examines the line l0 defined by the intersection of
   * planes i and j.  Assume triangles 1 and 2 both intersect this line, which
   * they must if they intersect each other.  Then the intersection of triangle
   * 1 with l0 is a segment s1, and the intersection of triangle 2 with l0 is
   * a segment s2.  If s1 and s2 overlap, triangles 1 and 2 intersect.  Hence
   * the name "intervalCheck".
   *
   * Helper function for T-T intersect.
   *****************************************************************************
   */
  inline bool intervalCheck(const Point3 Ai, const Point3 Bi, const Point3 Ci,
			    const Point3 Aj, const Point3 Bj, const Point3 Cj)
  {
    /* Step 5: From step's 1 through 4, we now have two triangles that,
       if intersecting, have a line that intersects segments AiCi, AiBi,
       AjBj, and AjCj.  We check if these two intervals overlap:
    */

    if (isGt(Vector3(Bi, Bj).dot(Triangle3(Bi, Aj, Ai).normal()), 0.0))
      return false;
    if (isGt((Vector3(Ai, Cj).dot(Triangle3(Ai, Aj, Ci).normal())), 0.0))
      return false;

    return true;
  }
}

  /*!
   *****************************************************************************
   * \brief Returns whether the two triangles are coplanar
   * \param [in] Ai first corner of triangle 1
   * \param [in] Bi second corner of triangle 1
   * \param [in] Ci third corner of triangle 1
   * \param [in] Aj first corner of triangle 2
   * \param [in] Bj second corner of triangle 2
   * \param [in] Cj third corner of triangle 2
   * \return status true iff triangle Ai Bi Ci is coplanar to triangle Aj Bj Cj
   *****************************************************************************
   */
  inline bool coplanarCheck(const Point3& Ai, const Point3& Bi, const Point3& Ci,
			    const Point3& Aj, const Point3& Bj, const Point3& Cj,
			    Vector3 normal)
  {
    /* Co-planar triangles are projected onto the axis that maximizes their
       area and the 2d intersection used to check if they intersect.
    */

    //find triangle with maximum area:
    for (int i=0; i<3; i++)
    {
      normal[i] = std::abs(normal[i]);
    }

    if ((isGt(normal[0], normal[2])) && (isGeq(normal[0], normal[1])))
    {
      //if x projection area greatest, project on YZ and return 2D checker

      const Triangle2 t1_2da = Triangle2(Point2::make_point(Bi[2],Bi[1]),
					 Point2::make_point(Ai[2],Ai[1]),
					 Point2::make_point(Ci[2],Ci[1]));

      const Triangle2 t2_2da = Triangle2(Point2::make_point(Bj[2],Bj[1]),
					 Point2::make_point(Aj[2],Aj[1]),
					 Point2::make_point(Cj[2],Cj[1]));

      return TriangleIntersection2D(t1_2da, t2_2da);
    }
    else if (isGt(normal[1],normal[2]) && isGeq(normal[1],normal[0]))
    {
      //if y projection area greatest, project on XZ and return 2D checker
      const Triangle2 t1_2da = Triangle2(Point2::make_point(Bi[0],Bi[2]),
					 Point2::make_point(Ai[0],Ai[2]),
					 Point2::make_point(Ci[0],Ci[2]));

      const Triangle2 t2_2da = Triangle2(Point2::make_point(Bj[0],Bj[2]),
					 Point2::make_point(Aj[0],Aj[2]),
					 Point2::make_point(Cj[0],Cj[2]));

      return TriangleIntersection2D(t1_2da, t2_2da);
    }
    else
    {
      //if z projection area greatest, project on XY and return 2D checker
      const Triangle2 t1_2da = Triangle2(Point2::make_point(Ai[0],Ai[1]),
					 Point2::make_point(Bi[0],Bi[1]),
					 Point2::make_point(Ci[0],Ci[1]));

      const Triangle2 t2_2da = Triangle2(Point2::make_point(Aj[0],Aj[1]),
					 Point2::make_point(Bj[0],Bj[1]),
					 Point2::make_point(Cj[0],Cj[1]));

      return TriangleIntersection2D(t1_2da, t2_2da);
    }
    return false;
  }

  /*!
   *****************************************************************************
   * \brief Worker function testing for 3D triangle intersection.
   * \param [in] Ai first corner of triangle 1
   * \param [in] Bi second corner of triangle 1
   * \param [in] Ci third corner of triangle 1
   * \param [in] Aj first corner of triangle 2
   * \param [in] Bj second corner of triangle 2
   * \param [in] Cj third corner of triangle 2
   * \param [in] dAj Dot product of normal with segment from t1 to Aj
   * \param [in] dBj Dot product of normal with segment from t1 to Bj
   * \param [in] dCj Dot product of normal with segment from t1 to Cj
   * \param [in] normal Normal vector of triangle 1
   * \return status true iff there is an interval overlap, false otherwise
   *
   * Bi and Ci both lie in the negative half-space defined by t2; Ai lies in
   * t2's plane or in its positive half-space.
   * The sign of dAj, dBj, and dCj indicates whether the associated vertex
   * of t2 lies in the positive or negative half-space defined by t1.
   *
   * Helper function for TT-intersect
   *****************************************************************************
   */
  inline bool checkTriangleIntersect(const Point3 &Ai, const Point3 &Bi, const Point3 &Ci,
				     const Point3 &Aj, const Point3 &Bj, const Point3 &Cj,
				     double dAj, double dBj, double dCj,  Vector3 &normal)
  {
    /*Step 4: repeat Step 3, except doing it for triangle 2 instead of triangle 1 */
    if (isGt(dAj, 0.0))
    {
      if (isGt(dBj, 0.0))
	return intervalCheck(Ai,Ci,Bi,Cj,Aj,Bj);
      else if (isGt(dCj, 0.0))
	return intervalCheck(Ai,Ci,Bi,Bj,Cj,Aj);
      else
	return intervalCheck(Ai,Bi,Ci,Aj,Bj,Cj);
    }
    else if (isLt(dAj,  0.0))
    {
      if (isLt(dBj, 0.0))
	return intervalCheck(Ai,Bi,Ci,Cj,Aj,Bj);
      else if (isLt(dCj, 0.0))
	return intervalCheck(Ai,Bi,Ci,Bj,Cj,Aj);
      else
	return intervalCheck(Ai,Ci,Bi,Aj,Bj,Cj);
    } else {
      if (isLt(dBj, 0.0)) {
	if (isGeq(dCj, 0.0))
	  return intervalCheck(Ai,Ci,Bi,Bj,Cj,Aj);
	else
	  return intervalCheck(Ai,Bi,Ci,Aj,Bj,Cj);
      }
      else if (isGt(dBj, 0.0)) {
	if (isGt(dCj, 0.0))
	  return intervalCheck(Ai,Ci,Bi,Aj,Bj,Cj);
	else
	  return intervalCheck(Ai,Bi,Ci,Bj,Cj,Aj);
      }
      else {
	if (isGt(dCj, 0.0))
	  return intervalCheck(Ai,Bi,Ci,Cj,Aj,Bj);
	else if (isLt(dCj, 0.0))
	  return intervalCheck(Ai,Ci,Bi,Cj,Aj,Bj);
	else
	  return coplanarCheck(Ai,Bi,Ci,Aj,Bj,Cj,normal);
      }
    }
  }


/*
 *******************************************************************************
 * \brief Tests if 2D Triangles t1 and t2 intersect.
 * \param [in] t1 First 2D triangle to test.
 * \param [in] t2 Second 2D triangle to test.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *******************************************************************************
 */
template < typename T>
bool intersect( const Triangle<T, 2>& t1, const Triangle<T, 2>& t2)
{
  if (t1.degenerate() || t2.degenerate()) {
    if (t1.degenerate())
      SLIC_INFO("\n\n WARNING \n\n Triangle " << t1 <<" is degenerate");
    if (t2.degenerate())
      SLIC_INFO("\n\n WARNING \n\n Triangle " << t2 <<" is degenerate");
  }
  return TriangleIntersection2D(t1, t2);
}


/*!
 *******************************************************************************
 * \brief Tests if 3D Triangles t1 and t2 intersect.
 * \param [in] t1 First 3D triangle to test.
 * \param [in] t2 Second 3D triangle to test.
 * \return status true iff t1 intersects with t2, otherwise, false.
 *******************************************************************************
 */
template < typename T>
bool intersect( const Triangle<T, 3>& t1, const Triangle<T, 3>& t2)
{
  typedef quest::Vector<T, 3> Vector3;

  /*............................................................................
    Compute the line of intersection between the two planes with a decision
    tree aproach.  Note that these nasty if statments are used for robustness
    to reduce numerical errors -- I found a nice paper which discussed this,
    and I modeled the following code after it.  By using a decision tree
    approach, they were able to reduce the number of arithmatic operations and
    thus reduce numerical imprecision.  What remains to be seen is whether the
    amount of overhead associated with our abstractions make this approach
    impractical for our application.

    See paper at https://hal.inria.fr/inria-00072100/document for more details.
    ............................................................................
  */
  if (t1.degenerate() || t2.degenerate()) {
    if (t1.degenerate())
      SLIC_INFO("\n\n WARNING \n\n Triangle " << t1 <<" is degenerate");
    if (t2.degenerate())
      SLIC_INFO("\n\n WARNING \n\n Triangle " << t2 <<" is degenerate");
  }

  // Step 1: Check if all the vertices of triangle 1 lay on the same side of
  // the plane created by triangle 2:

  Vector3 t2Normal = Vector3::cross_product(Vector3(t2.C(), t2.A()),
					    Vector3(t2.C(), t2.B()));
  double dAi = (Vector3(t2.C(), t1.A())).dot(t2Normal);
  double dBi = (Vector3(t2.C(),t1.B())).dot(t2Normal);
  double dCi = (Vector3(t2.C(),t1.C())).dot(t2Normal);
  if (signMatch(dAi, dBi, dCi)) {
    return false;
  }

  // Step 2: Check if all the vertices of triangle 2 lay on the same side of
  // the plane created by triangle 1:

  Vector3 t1Normal = Vector3::cross_product(Vector3(t1.A(), t1.B()),
					    Vector3(t1.A(), t1.C()));
  double dAj = (Vector3(t1.C(),t2.A())).dot(t1Normal);
  double dBj = (Vector3(t1.C(),t2.B())).dot(t1Normal);
  double dCj = (Vector3(t1.C(),t2.C())).dot(t1Normal);
  if (signMatch(dAj, dBj, dCj)) {
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

  if (isGt(dAi, 0.0)) {
    if (isGt(dBi, 0.0)) {
      return checkTriangleIntersect(t1.C(), t1.A(), t1.B(),
				    t2.A(), t2.C(), t2.B(),
				    dAj, dCj, dBj, t1Normal);
    }
    else if (isGt(dCi, 0.0)) {
      return checkTriangleIntersect(t1.B(), t1.C(), t1.A(),
				    t2.A(), t2.C(), t2.B(),
				    dAj, dCj, dBj, t1Normal);
    }
    else return checkTriangleIntersect(t1.A(), t1.B(), t1.C(),
				       t2.A(), t2.B(), t2.C(),
				       dAj, dBj, dCj, t1Normal);
  }
  else if (isLt(dAi, 0.0)) {
    if (isLt(dBi, 0.0)) {
      return checkTriangleIntersect(t1.C(), t1.A(), t1.B(),
				    t2.A(), t2.B(), t2.C(),
				    dAj, dBj, dCj, t1Normal);
    }
    else if (isLt(dCi, 0.0f)) {
      return checkTriangleIntersect(t1.B(), t1.C(), t1.A(),
				    t2.A(), t2.B(), t2.C(),
				    dAj, dBj, dCj, t1Normal);
    }
    else return checkTriangleIntersect(t1.A(), t1.B(), t1.C(),
				       t2.A(), t2.C(), t2.B(),
				       dAj, dCj, dBj, t1Normal);
  }
  else { //dAi ~= 0
    if (isLt(dBi, 0.0)) {
      if (isGeq(dCi, 0.0)) {
	return checkTriangleIntersect(t1.B(), t1.C(), t1.A(),
				      t2.A(), t2.C(), t2.B(),
				      dAj, dCj, dBj, t1Normal);
      }
      else {
	return checkTriangleIntersect(t1.A(), t1.B(), t1.C(),
				      t2.A(), t2.B(), t2.C(),
				      dAj, dBj, dCj, t1Normal);
      }
    }
    else if (isGt(dBi, 0.0)) {
      if (isGt(dCi, 0.0)) {
	return checkTriangleIntersect(t1.A(), t1.B(), t1.C(),
				      t2.A(), t2.C(), t2.B(),
				      dAj, dCj, dBj, t1Normal);
      }
      else {
	return checkTriangleIntersect(t1.B(), t1.C(), t1.A(),
				      t2.A(), t2.B(), t2.C(),
				      dAj, dBj, dCj, t1Normal);
      }
    }
    else  {
      if (isGt(dCi, 0.0)) {
	return checkTriangleIntersect(t1.C(), t1.A(), t1.B(),
				      t2.A(), t2.B(), t2.C(),
				      dAj, dBj, dCj, t1Normal);
      }
      else if (isLt(dCi, 0.0)) {
	return checkTriangleIntersect(t1.C(), t1.A(), t1.B(),
				      t2.A(), t2.C(), t2.B(),
				      dAj, dCj, dBj, t1Normal);
      }
      else return coplanarCheck(t1.A(), t1.B(), t1.C(),
				t2.A(), t2.B(), t2.C(), t1Normal);
    }
  }
}



namespace {

  typedef quest::Vector<double, 3> Vector3;

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
      if(d1 < d0)
          std::swap(d1,d0);  // d0 < d1
      if(d2 > d1)
          std::swap(d2,d1);  // d1 is max(d0,d1,d2)
      else if(d2 < d0)
          std::swap(d2,d0);  // d0 is min(d0,d1,d2)

      SLIC_ASSERT( d0 <= d1 && d0 <= d2);
      SLIC_ASSERT( d1 >= d0 && d1 >= d2);

      return d1 < -r || d0 > r;
  }

  /**
   * \brief Helper function for Triangle/BoundingBox intersection test
   */
  bool crossEdgesDisjoint(double d0, double d1, double r)
  {
      return std::max( -std::max(d0,d1), std::min(d0,d1) ) > r;
  }

}



/*!
 *******************************************************************************
 * \brief Determines if a triangle and a bounding box intersect
 *        (but does not find the point of intersection)
 * \param [in] tri user-supplied triangle (with three vertices).
 * \param [in] bb user-supplied axis aligned bounding box.
 * \return true iff tri intersects with bb, otherwise, false.
 *******************************************************************************
 */
template < typename T>
bool intersect( const Triangle<T, 3>& tri, const BoundingBox<T, 3>& bb)
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

    typedef typename BoundingBox<T,3>::PointType PointType;
    typedef typename BoundingBox<T,3>::VectorType VectorType;

    // Extent: vector center to max corner of BB
    VectorType e = 0.5 * bb.range();

    // Make the AABB center the origin by moving the triangle vertices
    PointType center(bb.getMin().array() + e.array());
    VectorType v[3] = { VectorType(center, tri.A())
                      , VectorType(center, tri.B())
                      , VectorType(center, tri.C()) };

    // Create the edge vectors of the triangle
    VectorType f[3] = { v[1] - v[0], v[2] - v[1],  v[0] - v[2] };


    // Test cross products of edges between triangle edge vectors f and cube normals (9 tests)
    // -- using separating axis theorem on the cross product of edges of triangle and face normals of AABB
    // Each test involves three cross products, two of which have the same value
    // The commented parameters highlights this symmetry.
    #define XEDGE_R( _E0, _E1, _F0, _F1, _IND )   e[ _E0 ] * std::abs(f[ _IND ][ _F0 ])            \
                                                + e[ _E1 ] * std::abs(f[ _IND ][ _F1 ])

    #define XEDGE_S( _V0, _V1, _F0, _F1, _VIND, _FIND) -v[ _VIND ][ _V0 ] * f[ _FIND ][ _F0 ]       \
                                                       +v[ _VIND ][ _V1 ] * f[ _FIND ][ _F1 ]

    if( crossEdgesDisjoint(/*XEDGE_S(1,2,2,1,0,0),*/ XEDGE_S(1,2,2,1,1,0),   XEDGE_S(1,2,2,1,2,0),   XEDGE_R(1,2,2,1,0))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(1,2,2,1,0,1),/* XEDGE_S(1,2,2,1,1,1),*/ XEDGE_S(1,2,2,1,2,1),   XEDGE_R(1,2,2,1,1))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(1,2,2,1,0,2),   XEDGE_S(1,2,2,1,1,2),/* XEDGE_S(1,2,2,1,2,2),*/ XEDGE_R(1,2,2,1,2))) return false;

    if( crossEdgesDisjoint(/*XEDGE_S(2,0,0,2,0,0),*/ XEDGE_S(2,0,0,2,1,0),   XEDGE_S(2,0,0,2,2,0),   XEDGE_R(0,2,2,0,0))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(2,0,0,2,0,1),/* XEDGE_S(2,0,0,2,1,1),*/ XEDGE_S(2,0,0,2,2,1),   XEDGE_R(0,2,2,0,1))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(2,0,0,2,0,2),   XEDGE_S(2,0,0,2,1,2),/* XEDGE_S(2,0,0,2,2,2),*/ XEDGE_R(0,2,2,0,2))) return false;

    if( crossEdgesDisjoint(/*XEDGE_S(0,1,1,0,0,0),*/ XEDGE_S(0,1,1,0,1,0),   XEDGE_S(0,1,1,0,2,0),   XEDGE_R(0,1,1,0,0))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(0,1,1,0,0,1),/* XEDGE_S(0,1,1,0,1,1),*/ XEDGE_S(0,1,1,0,2,1),   XEDGE_R(0,1,1,0,1))) return false;
    if( crossEdgesDisjoint(  XEDGE_S(0,1,1,0,0,2),   XEDGE_S(0,1,1,0,1,2),/* XEDGE_S(0,1,1,0,2,2),*/ XEDGE_R(0,1,1,0,2))) return false;

    #undef XEDGE_R
    #undef XEDEG_S


    /// Test face normals of bounding box (3 tests)
    if(intervalsDisjoint(v[0][0], v[1][0], v[2][0], e[0])) return false;
    if(intervalsDisjoint(v[0][1], v[1][1], v[2][1], e[1])) return false;
    if(intervalsDisjoint(v[0][2], v[1][2], v[2][2], e[2])) return false;


    /// Final test -- face normal of triangle's plane
    VectorType planeNormal  = VectorType::cross_product(f[0],f[1]);
    double planeDist    = planeNormal.dot( tri.A());

    double r = e[0]* std::abs( planeNormal[0]) + e[1]* std::abs( planeNormal[1]) + e[2]* std::abs( planeNormal[2]);
    double s = planeNormal.dot(center) - planeDist;

    return std::abs(s) <= r;
}


//------------------------------------------------------------------------------

/*!
 *******************************************************************************
 * \brief Determines if two axis aligned bounding boxes intersect
 * \param [in] bb1 user-supplied axis aligned bounding box.
 * \param [in] bb2 user-supplied axis aligned bounding box.
 * \return true iff bb1 intersects with bb2, otherwise, false.
 *******************************************************************************
 */
template < typename T, int DIM>
bool intersect( const BoundingBox<T, DIM>& bb1, const BoundingBox<T, DIM>& bb2)
{
    return bb1.intersects(bb2);
}


} /* end namespace quest */

#endif /* INTERSECTION_HPP_ */
