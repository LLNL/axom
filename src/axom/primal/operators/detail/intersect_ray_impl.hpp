// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_INTERSECT_RAY_HPP_
#define PRIMAL_INTERSECT_RAY_HPP_

// primal includes
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

namespace axom
{
namespace primal
{
namespace detail
{

/*!
 * \brief Computes the intersection of the given ray, R, with the segment, S.
 *      ip returns the intersection point on S.
 * \return status true iff R intersects with S, otherwise, false.
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
  const double denom = numerics::determinant(
    R.direction()[0], (-1.0)*R2.direction()[0],
    R.direction()[1], (-1.0)*R2.direction()[1]     );

  // STEP 2: if denom is zero, the system is singular, which implies that the
  // ray and the segment are parallel
  const double parepsilon = 1.0e-9;
  if ( axom::utilities::isNearlyEqual( denom, 0.0, parepsilon ) )
  {

    // ray and segment are parallel
    return false;

  }

  // STEP 3: Solve for t0 and t1 directly using cramer's rule
  const double alpha = S.source()[0] - R.origin()[0];
  const double beta  = S.source()[1] - R.origin()[1];

  const double t0 = numerics::determinant( alpha, (-1.0)*R2.direction()[0],
                                           beta,
                                           (-1.0)*R2.direction()[1] )/denom;

  const double t1 = numerics::determinant( R.direction()[0], alpha,
                                           R.direction()[1], beta   )/denom;

  // STEP 4: Define lower/upper threshold
  const double tlow  = 0.0-1.0e-9;
  const double thigh = 1.0+1.0e-9;

  // STEP 5: Necessary and sufficient criteria for an intersection between
  // ray, R(t0),  and a finite segment S(t1) are:
  // 1. t0 >= tlow w.r.t. the ray R(t0).
  // 2. tlow >= t1 >= thigh w.r.t. the segment S(t1).
  if ( (t0 >= tlow) && (t1 >= tlow) && (t1 <= thigh) )
  {
    ip = R2.at( t1 );
    return true;
  }

  // STEP 6: Ray does not intersect the segment
  return false;
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *      ip the point of intersection on R.
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * Computes Ray Box intersection using the slab method from pg 180 of
 * Real Time Collision Detection by Christer Ericson.
 */
template < typename T, int DIM >
bool intersect_ray_bbox(const primal::Ray< T,DIM > & R,
                        const primal::BoundingBox< T,DIM > & bb,
                        primal::Point< T,DIM > & ip)
{
  T tmin = std::numeric_limits< T >::min();
  SLIC_ASSERT(tmin>=0.0);
  T tmax = std::numeric_limits< T >::max();

  for (int i=0 ; i<DIM ; i++)
  {
    if (axom::utilities::isNearlyEqual(R.direction()[i],
                                       std::numeric_limits< T >::min(),
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
      T ood = (static_cast< T >(1.0)) / (R.direction()[i]);
      T t1 = ((bb.getMin()[i]- R.origin()[i])*ood);
      T t2 = ((bb.getMax()[i]- R.origin()[i])*ood);

      if (t1>t2)
      {
        std::swap(t1,t2);
      }

      tmin = axom::utilities::max(tmin, t1);
      tmax = axom::utilities::min(tmax, t2);

      if (tmin > tmax)
      {
        return false;
      }
    }
  }

  for (int i = 0 ; i < DIM ; i++)
  {
    ip.data()[i] = R.origin()[i] + R.direction()[i] * tmin;
  }

  return true;
}

} /* namespace detail */
} /* namespace primal */
} /* namespace axom */

#endif
