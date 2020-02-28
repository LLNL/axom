// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_INTERSECT_RAY_HPP_
#define PRIMAL_INTERSECT_RAY_HPP_

// numerics includes
#include "axom/core/numerics/floating_point_limits.hpp"

// primal includes
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

#include <type_traits> // for std::is_floating_point< T >()
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
inline bool intersect_ray( const primal::Ray< T,2 >& R,
                           const primal::Segment< T,2 >& S,
                           primal::Point< T,2 >& ip )
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

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
 * \brief Helper routine for ray / AABB intersection test
 *
 * \param [in] x0 coordinate component of the ray origin.
 * \param [in] n normal component of the ray direction.
 * \param [in] min the AABB min coordinate along a direction.
 * \param [in] max the AABB max coordinate along a drection.
 *
 * \param [in,out] tmin
 * \param [in,out] tmax
 *
 * \param [in] TOL
 *
 * \return
 */
template < typename T >
AXOM_HOST_DEVICE
inline bool intersect_ray_bbox_test( const T& x0,
                                     const T& n,
                                     const T& min,
                                     const T& max,
                                     T& tmin,
                                     T& tmax,
                                     T TOL )
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

  constexpr T ZERO  = 0.0f;

  bool status = true;

  if ( axom::utilities::isNearlyEqual( n, ZERO, TOL ) )
  {
    status = ( ( (x0 < min) || (x0 > max) ) ? false : true );
  }
  else
  {
    const T invn = static_cast< T >( 1.0 ) / n;
    T t1 = ( min - x0 ) * invn;
    T t2 = ( max - x0 ) * invn;

    if ( t1 > t2 )
    {
      axom::utilities::swap( t1, t2 );
    }

    tmin = axom::utilities::max( tmin, t1 );
    tmax = axom::utilities::min( tmax, t2 );

    status = ( ( tmin > tmax ) ? false : true );
  }

  return status;

}

/*!
 * \brief
 *
 * \param [in] x0
 * \param [in] n
 * \param [in] xmin
 * \param [in] xmax
 * \param [in] ymin
 * \param [in] ymax
 *
 * \param [out] t
 *
 * \param [in] TOL optional tolerance. Defaults to 1.e-9 if not specified.
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template < typename T >
AXOM_HOST_DEVICE
inline bool intersect_ray( const T* x0,
                           const T* n,
                           const T& xmin,
                           const T& xmax,
                           const T& ymin,
                           const T& ymax,
                           T& t,
                           T TOL=1.e-9 )
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

  SLIC_ASSERT( x0 != nullptr );
  SLIC_ASSERT( n != nullptr );

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits< T >::max();

  bool status = true;
  status = status && intersect_ray_bbox_test( x0[0],n[0],xmin,xmax,t,tmax,TOL );
  status = status && intersect_ray_bbox_test( x0[1],n[1],ymin,ymax,t,tmax,TOL );

  return status;
}

/*!
 * \brief
 *
 * \param [in] x0
 * \param [in] n
 * \param [in] xmin
 * \param [in] xmax
 * \param [in] ymin
 * \param [in] ymax
 * \param [in] zmin
 * \param [in] zmax
 *
 * \param [out] t
 *
 * \param [in] TOL optional tolerance. Defaults to 1.e-9 if not specified.
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template < typename T >
AXOM_HOST_DEVICE
inline bool intersect_ray( const T* x0,
                           const T* n,
                           const T& xmin,
                           const T& xmax,
                           const T& ymin,
                           const T& ymax,
                           const T& zmin,
                           const T& zmax,
                           T& t,
                           T TOL=1.e-9)
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

  SLIC_ASSERT( x0 != nullptr );
  SLIC_ASSERT( n != nullptr );

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits< T >::max();

  bool status = true;
  status = status && intersect_ray_bbox_test( x0[0],n[0],xmin,xmax,t,tmax,TOL );
  status = status && intersect_ray_bbox_test( x0[1],n[1],ymin,ymax,t,tmax,TOL );
  status = status && intersect_ray_bbox_test( x0[2],n[2],zmin,zmax,t,tmax,TOL );

  return status;
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
inline bool intersect_ray( const primal::Ray< T,DIM > & R,
                           const primal::BoundingBox< T,DIM > & bb,
                           primal::Point< T,DIM > & ip,
                           T TOL=1.e-9 )
{
  AXOM_STATIC_ASSERT( std::is_floating_point< T >::value );

  T tmin = std::numeric_limits< T >::min();
  T tmax = std::numeric_limits< T >::max();

  bool intersects = true;
  for (int i=0 ; ( intersects && (i < DIM) ) ; i++)
  {
    intersects = intersect_ray_bbox_test( R.origin()[ i ],
                                          R.direction()[ i ],
                                          bb.getMin()[ i ],
                                          bb.getMax()[ i ],
                                          tmin,
                                          tmax,
                                          TOL );

  }

  if ( intersects )
  {
    for (int i = 0 ; i < DIM ; i++)
    {
      ip.data()[i] = R.origin()[i] + R.direction()[i] * tmin;
    }
  }

  return intersects;
}

} /* namespace detail */
} /* namespace primal */
} /* namespace axom */

#endif
