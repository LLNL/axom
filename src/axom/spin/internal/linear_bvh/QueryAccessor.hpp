// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_QUERYACCESSOR_HPP_
#define AXOM_SPIN_BVH_QUERYACCESSOR_HPP_

#include "axom/core/Macros.hpp" // for axom macros
#include "axom/core/Types.hpp"  // for axom::IndexType

// C/C++ includes
#include <type_traits> // for std::is_floating_point()

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
template < int NDIMS, typename FloatType >
class QueryAccessor
{
public:

  AXOM_STATIC_ASSERT_MSG( ( (NDIMS==2) || (NDIMS==3) ),
                          "TraversalPredicates are defined for 2D and 3D only" );
  AXOM_STATIC_ASSERT_MSG( std::is_floating_point< FloatType >::value,
                          "A valid FloatingType must be used, e.g., double or float" );

  /// \name Query Point Access methods
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
  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline void getPoint( PointType& point,
                               IndexType idx,
                               const FloatType* x,
                               const FloatType* y,
                               const FloatType* z );

  /// @}
};

//------------------------------------------------------------------------------
// 2D Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class QueryAccessor< 2, FloatType >
{
public:

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline void getPoint( PointType& point,
                               IndexType idx,
                               const FloatType* x,
                               const FloatType* y,
                               const FloatType* AXOM_NOT_USED(z) )
  {
    point[ 0 ] = x[ idx ];
    point[ 1 ] = y[ idx ];
  }

};

//------------------------------------------------------------------------------
// 3D Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class QueryAccessor< 3, FloatType >
{
public:

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline void getPoint( PointType& point,
                               IndexType idx,
                               const FloatType* x,
                               const FloatType* y,
                               const FloatType* z )
  {
    point[ 0 ] = x[ idx ];
    point[ 1 ] = y[ idx ];
    point[ 2 ] = z[ idx ];
  }

};


} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */



#endif /* AXOM_SPIN_BVH_QUERYACCESSOR_HPP_ */
