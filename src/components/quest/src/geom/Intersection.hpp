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

#include "quest/Determinants.hpp"
#include "quest/fuzzy_compare.hpp"
#include "quest/Point.hpp"
#include "quest/Ray.hpp"
#include "quest/Segment.hpp"

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

//------------------------------------------------------------------------------

} /* end namespace quest */

#endif /* INTERSECTION_HPP_ */
