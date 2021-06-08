// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_QUERYACCESSOR_HPP_
#define AXOM_SPIN_BVH_QUERYACCESSOR_HPP_

#include "axom/core/Macros.hpp"  // for axom macros
#include "axom/core/Types.hpp"   // for axom::IndexType

// C/C++ includes
#include <type_traits>  // for std::is_floating_point()

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{
/*!
 * \brief QueryAccessor is a singleton class consisting of methods to extract
 *  the query input information for different input types, e.g., points, rays,
 *  etc.
 *
 *  \tparam NDIMS the spatial dimension
 *  \tparam FloatType the floating point precision, e.g., double or float
 *
 *  \pre NDIMS==2 || NDIMS==3
 *  \pre std::is_floating_point< FloatType >::value
 */
template <int NDIMS, typename FloatType>
class QueryAccessor
{
public:
  AXOM_STATIC_ASSERT_MSG(((NDIMS == 2) || (NDIMS == 3)),
                         "TraversalPredicates are defined for 2D and 3D only");
  AXOM_STATIC_ASSERT_MSG(
    std::is_floating_point<FloatType>::value,
    "A valid FloatingType must be used, e.g., double or float");

  /// \name Point Query Access Methods
  /// @{

  /*!
   * \brief Gets the coordinates of a query point from the supplied coordinates.
   *
   * \param [out] point the point where to pack the coordinates
   * \param [in] idx the index of the query point w.r.t. the supplied arrays.
   * \param [in] x user-supplied array of x-coordinates
   * \param [in] y user-supplied array of y-coordinates
   * \param [in] z user-supplied array of z-coordinates (3D only)
   *
   * \tparam PointType the data-structure to hold a point.
   *
   * \note The PointType must overload the `[]` operator for accessing the
   *  point coordinates.
   */
  template <typename PointType>
  AXOM_HOST_DEVICE static inline void getPoint(PointType& point,
                                               IndexType idx,
                                               const FloatType* x,
                                               const FloatType* y,
                                               const FloatType* z);

  /// @}

  /// \name Ray Query Access Methods
  /// @{

  /*!
   * \brief Gets the ray defined by a source point and normal direction.
   *
   * \param [out] ray the ray where to pack the source point and normal.
   * \param [in] idx the index of the ray w.r.t. the supplied arrays.
   * \param [in] x0 array of the x-coordinate source point locations.
   * \param [in] nx array of x-components of the ray normal directions.
   * \param [in] y0 array of the y-coordinate source point locations.
   * \param [in] ny array of y-components of the ray normal directions.
   * \param [in] z0 array of of the z-coorindate source point locations.
   * \param [in] nz array of z-components of the ray normal directions.
   *
   * \tparam RayType the data-structure to hold a ray.
   *
   * \note The RayType must overload the `[]` operator for accessing the
   *  ray source point coordinates and normal direction information, such that,
   *  the first NDIMS locations store the source points coordinates followed
   *  by the normal direction.
   */
  template <typename RayType>
  AXOM_HOST_DEVICE static inline void getRay(RayType& ray,
                                             IndexType idx,
                                             const FloatType* x0,
                                             const FloatType* nx,
                                             const FloatType* y0,
                                             const FloatType* ny,
                                             const FloatType* z0,
                                             const FloatType* nz);
  /// @}

  /// \name Bounding Box Query Access Methods
  /// @{

  /*!
   * \brief Gets the bounding box defined by lower and upper bound corners.
   *
   * \param [out] box the bounding box where to pack the min and max point
   * \param [in] idx the index of the bounding box w.r.t. the supplied arrays.
   * \param [in] xmin array of x-coordinate of the lower bounding box corner
   * \param [in] xmax array of x-coordinate of the upper bounding box corner
   * \param [in] ymin array of y-coordinate of the lower bounding box corner
   * \param [in] ymax array of y-coordinate of the upper bounding box corner
   * \param [in] zmin array of z-coordinate of the lower bounding box corner
   * \param [in] zmax array of z-coordinate of the upper bounding box corner
   *
   * \tparam BoundingBoxType the data-structure to hold a Bounding Box.
   *
   * \note The BoundingBoxType must overload the `[]` operator for accessing
   *  the lower and upper bound corner coordinates, such that,
   *  the first NDIMS locations store the lower bound coordinate min followed
   *  by the upper bound coordinate max.
   */
  template <typename BoundingBoxType>
  AXOM_HOST_DEVICE static inline void getBoundingBox(BoundingBoxType& box,
                                                     IndexType idx,
                                                     const FloatType* xmin,
                                                     const FloatType* xmax,
                                                     const FloatType* ymin,
                                                     const FloatType* ymax,
                                                     const FloatType* zmin,
                                                     const FloatType* zmax);
  /// @}
};

//------------------------------------------------------------------------------
// 2D Specialization
//------------------------------------------------------------------------------
template <typename FloatType>
class QueryAccessor<2, FloatType>
{
public:
  template <typename PointType>
  AXOM_HOST_DEVICE static inline void getPoint(PointType& point,
                                               IndexType idx,
                                               const FloatType* x,
                                               const FloatType* y,
                                               const FloatType* AXOM_NOT_USED(z))
  {
    point[0] = x[idx];
    point[1] = y[idx];
  }

  template <typename RayType>
  AXOM_HOST_DEVICE static inline void getRay(RayType& ray,
                                             IndexType idx,
                                             const FloatType* x0,
                                             const FloatType* nx,
                                             const FloatType* y0,
                                             const FloatType* ny,
                                             const FloatType* AXOM_NOT_USED(z0),
                                             const FloatType* AXOM_NOT_USED(nz))
  {
    ray[0] = x0[idx];
    ray[1] = y0[idx];

    ray[2] = nx[idx];
    ray[3] = ny[idx];
  }

  template <typename BoundingBoxType>
  AXOM_HOST_DEVICE static inline void getBoundingBox(
    BoundingBoxType& box,
    IndexType idx,
    const FloatType* xmin,
    const FloatType* xmax,
    const FloatType* ymin,
    const FloatType* ymax,
    const FloatType* AXOM_NOT_USED(zmin),
    const FloatType* AXOM_NOT_USED(zmax))
  {
    using PointType = typename BoundingBoxType::PointType;

    box = BoundingBoxType(PointType {xmin[idx], ymin[idx]},
                          PointType {xmax[idx], ymax[idx]});
  }
};

//------------------------------------------------------------------------------
// 3D Specialization
//------------------------------------------------------------------------------
template <typename FloatType>
class QueryAccessor<3, FloatType>
{
public:
  template <typename PointType>
  AXOM_HOST_DEVICE static inline void getPoint(PointType& point,
                                               IndexType idx,
                                               const FloatType* x,
                                               const FloatType* y,
                                               const FloatType* z)
  {
    point[0] = x[idx];
    point[1] = y[idx];
    point[2] = z[idx];
  }

  template <typename RayType>
  AXOM_HOST_DEVICE static inline void getRay(RayType& ray,
                                             IndexType idx,
                                             const FloatType* x0,
                                             const FloatType* nx,
                                             const FloatType* y0,
                                             const FloatType* ny,
                                             const FloatType* z0,
                                             const FloatType* nz)
  {
    ray[0] = x0[idx];
    ray[1] = y0[idx];
    ray[2] = z0[idx];

    ray[3] = nx[idx];
    ray[4] = ny[idx];
    ray[5] = nz[idx];
  }

  template <typename BoundingBoxType>
  AXOM_HOST_DEVICE static inline void getBoundingBox(BoundingBoxType& box,
                                                     IndexType idx,
                                                     const FloatType* xmin,
                                                     const FloatType* xmax,
                                                     const FloatType* ymin,
                                                     const FloatType* ymax,
                                                     const FloatType* zmin,
                                                     const FloatType* zmax)
  {
    using PointType = typename BoundingBoxType::PointType;

    box = BoundingBoxType(PointType {xmin[idx], ymin[idx], zmin[idx]},
                          PointType {xmax[idx], ymax[idx], zmax[idx]});
  }
};

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_BVH_QUERYACCESSOR_HPP_ */
