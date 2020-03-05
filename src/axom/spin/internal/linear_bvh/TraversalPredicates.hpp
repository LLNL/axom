// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVH_TRAVERSALPREDICATES_HPP_
#define AXOM_SPIN_BVH_TRAVERSALPREDICATES_HPP_

#include "axom/spin/internal/linear_bvh/vec.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{

template < typename FloatType >
using VecType = internal::linear_bvh::Vec< FloatType, 4 >;

/*!
 * \brief TraversalPredicates is a singleton class that defines different
 *  predicates for traversing a BVH instance that operate directly on the
 *  BVH data layout.
 *
 *  \see BVHData
 *
 *  \tparam NDIMS the spatial dimension
 *  \tparam FloatType the floating point type precision, e.g., double or float.
 */
template < int NDIMS, typename FloatType >
class TraversalPredicates
{

public:

  AXOM_STATIC_ASSERT_MSG( ( (NDIMS==2) || (NDIMS==3) ),
                          "TraversalPredicates are defined for 2D and 3D only" );
  AXOM_STATIC_ASSERT_MSG( std::is_floating_point< FloatType >::value,
                          "A valid FloatingType must be used, e.g., double or float" );

  ///\name Predicates for Point Queries
  /// @{

  /*!
   * \brief Checks if the supplied point is within the left bin.
   *
   * \param [in] point the coordinates of the point in query.
   * \param [in] s1 the 1st segment of the BVH that stores the left bin.
   * \param [in] s2 the 2nd segment of the BVH that stores the left bin.
   *
   * \return status true if the point is inside the left bin, else, false.
   *
   */
  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInLeftBin( const PointType& point,
                                     const VecType< FloatType >& s1,
                                     const VecType< FloatType >& s2 ) noexcept;

  /*!
   * \brief Checks if the supplied point is within the right bin.
   *
   * \param [in] point the coordinates of the point in query.
   * \param [in] s2 the 2nd segment of the BVH that stores the right bin.
   * \param [in] s3 the 3rd segment of the BVH that stores the right bin.
   *
   * \return status true if the point is inside the right bin, else, false.
   */
  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInRightBin( const PointType& point,
                                      const VecType< FloatType >& s2,
                                      const VecType< FloatType >& s3 ) noexcept;

  /// @}

};

constexpr int DIMENSION_2 = 2;
constexpr int DIMENSION_3 = 3;

//------------------------------------------------------------------------------
// 2D Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class TraversalPredicates< DIMENSION_2, FloatType >
{

public:

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInLeftBin( const PointType& point,
                                     const VecType< FloatType >& s1,
                                     const VecType< FloatType >& s2  ) noexcept
  {
    // NOTE: See BVHData.hpp for how the BVH bin is organized in the segments.
    bool in_left = true;

    if ( point[0] < s1[0] ) in_left = false;
    if ( point[1] < s1[1] ) in_left = false;

    if ( point[0] > s1[3] ) in_left = false;
    if ( point[1] > s2[0] ) in_left = false;

    return in_left;
  }

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInRightBin( const PointType& point,
                                      const VecType< FloatType >& s2,
                                      const VecType< FloatType >& s3  ) noexcept
  {
    // NOTE: See BVHData.hpp for how the BVH bin is organized in the segments.
    bool in_right = true;

    if ( point[0] < s2[2] ) in_right = false;
    if ( point[1] < s2[3] ) in_right = false;

    if ( point[0] > s3[1] ) in_right = false;
    if ( point[1] > s3[2] ) in_right = false;

    return in_right;
  }

};

//------------------------------------------------------------------------------
// 3D Specialization
//------------------------------------------------------------------------------
template < typename FloatType >
class TraversalPredicates< DIMENSION_3, FloatType >
{

public:

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInLeftBin( const PointType& point,
                                     const VecType< FloatType >& s1,
                                     const VecType< FloatType >& s2  ) noexcept
  {

    // NOTE: See BVHData.hpp for how the BVH bin is organized in the segments.
    bool in_left = true;

    if ( point[0] < s1[0] ) in_left = false;
    if ( point[1] < s1[1] ) in_left = false;
    if ( point[2] < s1[2] ) in_left = false;

    if ( point[0] > s1[3] ) in_left = false;
    if ( point[1] > s2[0] ) in_left = false;
    if ( point[2] > s2[1] ) in_left = false;

    return in_left;
  }

  template < typename PointType >
  AXOM_HOST_DEVICE
  static inline bool pointInRightBin( const PointType& point,
                                      const VecType< FloatType >& s2,
                                      const VecType< FloatType >& s3 ) noexcept
  {
    // NOTE: See BVHData.hpp for how the BVH bin is organized in the segments.
    bool in_right = true;

    if ( point[0] < s2[2] ) in_right = false;
    if ( point[1] < s2[3] ) in_right = false;
    if ( point[2] < s3[0] ) in_right = false;

    if ( point[0] > s3[1] ) in_right = false;
    if ( point[1] > s3[2] ) in_right = false;
    if ( point[2] > s3[3] ) in_right = false;

    return in_right;
  }

};

} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */


#endif /* AXOM_SPIN_BVH_TRAVERSALPREDICATES_HPP_ */
