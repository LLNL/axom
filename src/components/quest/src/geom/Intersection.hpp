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
#include "quest/fuzzy_compare.hpp"
#include "quest/Point.hpp"
#include "quest/Ray.hpp"
#include "quest/Segment.hpp"
#include "quest/Triangle.hpp"

namespace quest {

/*!
 *******************************************************************************
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 * \param [in] R user-supplied ray R.
 * \param [in] S user-supplied segment S.
 * \param [in/out] ip the point of intersection.
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
   if ( math::fuzzy_compare( denom, 0.0, 1.0e-9 ) ) {

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


namespace {

  typedef quest::Vector<double, 3> Vector3;

  /**
   * \brief Helper function to find disjoint projections for the AABB-triangle test
   * \param {d0,d1,d2}  Values defining the test interval
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
    VectorType e = bb.range() / 2.;

    // Make the AABB center the origin by moving the triangle vertices
    PointType center = bb.centroid();
    VectorType v[3] = { VectorType(tri.A().array() - center.array())
                      , VectorType(tri.B().array() - center.array())
                      , VectorType(tri.C().array() - center.array()) };

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
