/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef PRIMAL_SPHERE_HPP_
#define PRIMAL_SPHERE_HPP_

#include "axom/Macros.hpp"           // for AXOM macros
#include "axom/Types.hpp"            // for AXOM_NULLPTR
#include "axom_utils/Utilities.hpp"  // for utilities::isNearlyEqual()
#include "primal/OrientedSide.hpp"   // for OrientedSide enum definition

#include "slic/slic.hpp"             // for SLIC macros

namespace axom
{
namespace primal
{

/*!
 * \class Sphere
 *
 * \brief Defines an oriented Sphere in 2-D (i.e., a circle) or 3-D given by
 *  its center, \f$ \mathcal{X} \f$ and radius \f$ \mathcal{R} \f$. The Sphere
 *  object provides associated operations on a sphere, such as, signed distance
 *  and orientation.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template < typename T, int NDIMS >
class Sphere
{
public:

/// \name Constructors
/// @{

  /*!
   * \brief Constructs a Sphere centered at origin with the given radius.
   * \param [in] radius the radius of the Sphere (optional).
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  explicit Sphere( T radius=1.0 );

  /*!
   * \brief Constructs a Sphere with the given center and radius.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Sphere (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   *
   * \pre center != AXOM_NULLPTR
   */
  explicit Sphere( const T* center, T radius=1.0 );

/// @}

  /*!
   * \brief Destructor.
   */
  ~Sphere();

  /*!
   * \brief Returns the radius of the Sphere.
   * \return r the radius of the Sphere.
   */
  inline T getRadius( ) const { return m_radius; };

  /*!
   * \brief Returns the center of the Sphere.
   *
   * \return c pointer to array that holds the center of the Sphere.
   * \note c points to an array that is NDIMS long.
   * \post c != AXOM_NULLPTR
   */
  inline const T* getCenter( ) const { return m_center; };

  /*!
   * \brief Computes the signed distance of a point to the Sphere's boundary.
   *
   * \param [in] q pointer to buffer consisting of query point coordinates.
   * \return d the computed signed distance of the point, q, to the sphere.
   *
   * \note The signed distance of a point, q, is:
   *  <ul>
   *   <li> negative inside the sphere </li>
   *   <li> positive outside the sphere </li>
   *   <li> zero on the boundary </li>
   *  </ul>
   *
   * \pre q != AXOM_NULLPTR
   * \pre q must be a pointer to an array that is at least NDIMS long
   */
  inline T computeSignedDistance( const T* q ) const;

  /*!
   * \brief Computes the orientation of a point with respect to the Sphere.
   *
   * \param [in] q pointer to user-supplied point q.
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orient the orientation of q with respect to the sphere.
   *
   *  \note This method returns one of the following values:
   *   <ul>
   *    <li> <b>ON_BOUNDARY</b>      : if `q` is on the sphere's boundary </li>
   *    <li> <b>ON_POSITIVE_SIDE</b> : if `q` is outside the sphere </li>
   *    <li> <b>ON_NEGATIVE_SIDE</b> : if `q` is inside the sphere </li>
   *  </ul>
   *
   * \see OrientedSide for the list of possible return values.
   *
   * \pre q != AXOM_NULLPTR
   * \pre q must be a pointer to an array that is at least NDIMS long
   *
   */
  inline int getOrientation( const T* q, double TOL=1.e-9 ) const;

private:
  T m_center[ NDIMS ]; /*!< sphere center */
  T m_radius;          /*!< sphere radius */

  DISABLE_COPY_AND_ASSIGNMENT( Sphere );
  DISABLE_MOVE_AND_ASSIGNMENT( Sphere );
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
// Sphere Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Sphere< T,NDIMS >::Sphere( T radius ) : m_radius(radius)
{
  AXOM_STATIC_ASSERT_MSG( (NDIMS==2) || (NDIMS==3),
                            "A Sphere object may be defined in 2-D or 3-D" );

  for ( int i=0; i < NDIMS; ++i )
  {
    m_center[ i ] = static_cast< T >( 0.0 );
  }
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Sphere< T,NDIMS >::Sphere( const T* center, T radius ) : m_radius(radius)
{
  AXOM_STATIC_ASSERT_MSG( (NDIMS==2) || (NDIMS==3),
                            "A Sphere object may be defined in 2-D or 3-D" );

  SLIC_ASSERT( center != AXOM_NULLPTR );
  for ( int i=0; i < NDIMS; ++i )
  {
    m_center[ i ] = center[ i ];
  }

}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Sphere< T,NDIMS >::~Sphere() { }

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline T Sphere< T,NDIMS >::computeSignedDistance( const T* q ) const
{
  SLIC_ASSERT( q != AXOM_NULLPTR );

  T d = 0.0;
  for ( int i=0 ; i < NDIMS ; ++i )
  {
    const T dx = q[ i ]-m_center[ i ];
    d += (dx*dx);
  }

  return( std::sqrt( d )-m_radius );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline int Sphere< T,NDIMS >::getOrientation( const T* q, double TOL ) const
{
  SLIC_ASSERT( q != AXOM_NULLPTR );

  T signed_distance = this->computeSignedDistance( q );

  int orient = -1;

  if ( axom::utilities::isNearlyEqual( signed_distance, 0.0, TOL) )
  {

    orient = ON_BOUNDARY;

  }
  else if ( signed_distance < 0.0f )
  {

    // inside
    orient = ON_NEGATIVE_SIDE;

  }
  else
  {

    // outside
    orient = ON_POSITIVE_SIDE;

  }

  return orient;
}

} /* namespace primal */
} /* namespace axom */

#endif /* PRIMAL_SPHERE_HPP_ */
