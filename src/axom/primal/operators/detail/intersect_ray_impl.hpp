// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

#include <type_traits>  // for std::is_floating_point< T >()

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Computes the intersection of the given ray \a R with the segment \a S.
 *
 * When there is a valid intersection (within toleranace \a EPS), \a ray_param 
 * returns the parametric coordinate of the intersection point along \a R
 *
 * \return status true iff R intersects with S, otherwise, false.
 */
template <typename T>
inline bool intersect_ray(const primal::Ray<T, 2>& R,
                          const primal::Segment<T, 2>& S,
                          T& ray_param,
                          const double EPS)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  // STEP 0: Find the parameterized direction of the segment
  const auto seg_dir = primal::Vector<T, 2>(S.source(), S.target());
  const auto& ray_dir = R.direction();

  // Step 1: Equating R(ray_param)=S(seg_param) yields a system of two equations and
  // two unknowns, namely, t0 and t1. We can solve this system directly using Cramer's Rule.
  const double denom =
    numerics::determinant(ray_dir[0], -seg_dir[0], ray_dir[1], -seg_dir[1]);

  // STEP 2: if denom is (nearly) zero (within tolerance EPS), the system is singular
  if(axom::utilities::isNearlyEqual(denom, 0.0, EPS))
  {
    // ray and segment are parallel
    return false;
  }

  // STEP 3: Solve for the ray_param and seg_param directly using cramer's rule
  const auto sol = S.source().array() - R.origin().array();

  // Note: ray_param is an OUT parameter of this function
  ray_param =
    numerics::determinant(sol[0], -seg_dir[0], sol[1], -seg_dir[1]) / denom;

  const T seg_param =
    numerics::determinant(ray_dir[0], sol[0], ray_dir[1], sol[1]) / denom;

  // STEP 4: Define lower/upper threshold
  const double tlow = 0.0 - EPS;
  const double thigh = 1.0 + EPS;

  // STEP 5: Necessary and sufficient criteria for an intersection between
  // ray, R(t0),  and a finite segment S(t1) are:
  // 1. ray_param >= tlow w.r.t. the ray R(ray_param).
  // 2. tlow >= seg_param >= thigh w.r.t. the segment S(seg_param).
  return ((ray_param >= tlow) && (seg_param >= tlow) && (seg_param <= thigh));
}

/*!
 * \brief Helper routine for ray / AABB intersection test.
 *
 * \param [in] x0 coordinate component of the ray origin.
 * \param [in] n normal component of the ray direction.
 * \param [in] min the AABB min coordinate along a direction.
 * \param [in] max the AABB max coordinate along a drection.
 *
 * \param [in,out] tmin length of closest intersection point (computed)
 * \param [in,out] tmax length of farthest intersection point (computed)
 *
 * \param [in] TOL user-supplied tolerance.
 *
 * \return status true if the ray intersects, otherwise false.
 *
 * \note This routine is called by the intersect ray/AABB methods for each
 *  spatial dimension.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray_bbox_test(const T& x0,
                                                     const T& n,
                                                     const T& min,
                                                     const T& max,
                                                     T& tmin,
                                                     T& tmax,
                                                     T TOL)
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  constexpr T ZERO = 0.0f;

  bool status = true;

  if(axom::utilities::isNearlyEqual(n, ZERO, TOL))
  {
    status = (((x0 < min) || (x0 > max)) ? false : true);
  }
  else
  {
    const T invn = static_cast<T>(1.0) / n;
    T t1 = (min - x0) * invn;
    T t2 = (max - x0) * invn;

    if(t1 > t2)
    {
      axom::utilities::swap(t1, t2);
    }

    tmin = axom::utilities::max(tmin, t1);
    tmax = axom::utilities::min(tmax, t2);

    status = ((tmin > tmax) ? false : true);
  }

  return status;
}

/*!
 * \brief Checks if a specified ray intersects with the supplied bounding box.
 *
 * \param [in] x0 x-coordinate of the ray source point
 * \param [in] nx x-component of the ray normal
 * \param [in] y0 y-coordinate of the ray source point
 * \param [in] ny y-component of the ray normal
 * \param [in] xmin x-coordinate of the lower bounding box corner
 * \param [in] xmax x-coordinate of the upper bounding box corner
 * \param [in] ymin y-coordinate of the lower bounding box corner
 * \param [in] ymax y-coordinate of the upper bounding box corner
 * \param [out] t length of intersection point from the ray source point.
 * \param [in] TOL optional tolerance. Defaults to 1.e-9 if not specified.
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray(
  const T& x0,
  const T& nx,
  const T& y0,
  const T& ny,
  const T& xmin,
  const T& xmax,
  const T& ymin,
  const T& ymax,
  T& t,
  T TOL = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool status = true;
  status = status && intersect_ray_bbox_test(x0, nx, xmin, xmax, t, tmax, TOL);
  status = status && intersect_ray_bbox_test(y0, ny, ymin, ymax, t, tmax, TOL);

  return status;
}

/*!
 * \brief Checks if a specified ray intersects with the supplied bounding box.
 *
 * \param [in] x0 x-coordinate of the ray source point
 * \param [in] nx x-component of the ray normal
 * \param [in] y0 y-coordinate of the ray source point
 * \param [in] ny y-component of the ray normal
 * \param [in] z0 z-coordinate of the ray source point
 * \param [in] nz z-component of the ray normal
 * \param [in] xmin x-coordinate of the lower bounding box corner
 * \param [in] xmax x-coordinate of the upper bounding box corner
 * \param [in] ymin y-coordinate of the lower bounding box corner
 * \param [in] ymax y-coordinate of the upper bounding box corner
 * \param [in] zmin z-coordinate of the lower bounding box corner
 * \param [in] zmax z-coordinate of the upper bounding box corner
 * \param [out] t length of intersection point from the ray source point.
 * \param [in] TOL optional tolerance. Defaults to 1.e-9 if not specified.
 *
 * \return status true if the ray intersects the bounding box, otherwise, false.
 */
template <typename T>
AXOM_HOST_DEVICE inline bool intersect_ray(
  const T& x0,
  const T& nx,
  const T& y0,
  const T& ny,
  const T& z0,
  const T& nz,
  const T& xmin,
  const T& xmax,
  const T& ymin,
  const T& ymax,
  const T& zmin,
  const T& zmax,
  T& t,
  T TOL = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  t = 0.0f;
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool status = true;
  status = status && intersect_ray_bbox_test(x0, nx, xmin, xmax, t, tmax, TOL);
  status = status && intersect_ray_bbox_test(y0, ny, ymin, ymax, t, tmax, TOL);
  status = status && intersect_ray_bbox_test(z0, nz, zmin, zmax, t, tmax, TOL);

  return status;
}

/*!
 * \brief Computes the intersection of the given ray, R, with the Box, bb.
 *
 * \param [in] R the specified ray
 * \param [in] bb the user-supplied axis-aligned bounding box
 *
 * \param [out] ip the intersection point where R intersects bb.
 *
 * \return status true iff bb intersects with R, otherwise, false.
 *
 * \see primal::Ray
 * \see primal::Segment
 * \see primal::BoundingBox
 *
 * \note Computes Ray Box intersection using the slab method from pg 180 of
 *  Real Time Collision Detection by Christer Ericson.
 */
template <typename T, int DIM>
AXOM_HOST_DEVICE inline bool intersect_ray(
  const primal::Ray<T, DIM>& R,
  const primal::BoundingBox<T, DIM>& bb,
  primal::Point<T, DIM>& ip,
  T TOL = numerics::floating_point_limits<T>::epsilon())
{
  AXOM_STATIC_ASSERT(std::is_floating_point<T>::value);

  T tmin = axom::numerics::floating_point_limits<T>::min();
  T tmax = axom::numerics::floating_point_limits<T>::max();

  bool intersects = true;
  for(int i = 0; (intersects && (i < DIM)); i++)
  {
    intersects = intersect_ray_bbox_test(R.origin()[i],
                                         R.direction()[i],
                                         bb.getMin()[i],
                                         bb.getMax()[i],
                                         tmin,
                                         tmax,
                                         TOL);
  }

  if(intersects)
  {
    for(int i = 0; i < DIM; i++)
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
