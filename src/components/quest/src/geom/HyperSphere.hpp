/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef HYPERSPHERE_HPP_
#define HYPERSPHERE_HPP_

#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "quest/Orientation.hpp"
#include "slic/slic.hpp"

#include "common/Utilities.hpp"

// C/C++ includes
#include <cstddef>

namespace quest
{

/*!
 *******************************************************************************
 * \class HyperSphere
 *
 * \brief The HyperSphere class provides the means to represent a hypersphere,
 *  \f$ \mathcal{H} \in \mathcal{R}^d \f$, given a center, \f$ \mathcal{C} \f$
 *  and radius \f$ R \f$ namely a circle in \f$ \mathcal{R}^2 \f$ and sphere in
 *  \f$ \mathcal{R}^3 \f$ Further, the class provides the means to perform
 *  distance and inside/outside queries.
 *******************************************************************************
 */
template< typename T, int NDIMS >
class HyperSphere
{
public:

   /*!
    ****************************************************************************
    * \brief Constructs HyperSphere centered at origin with the given radius.
    * \param [in] radius radius of the HyperSphere. Default is 1.0
    ****************************************************************************
    */
   HyperSphere( T radius=1.0 );

   /*!
    ****************************************************************************
    * \brief Constructs HyperSphere with given center and radius
    * \param [in] center user-supplied center.
    * \param [in] radius user-supplied radius. Default is 1.0.
    ****************************************************************************
    */
   HyperSphere( T* center, T radius=1.0 );

   /*!
    ****************************************************************************
    * \brief Copy constructor.
    * \param [in] other The hypersphere to copy
    ****************************************************************************
    */
   HyperSphere( const HyperSphere<T,NDIMS>& other ) { *this = other; };

   /*!
    ****************************************************************************
    * \brief Destructor.
    ****************************************************************************
    */
   ~HyperSphere();

   /*!
    ****************************************************************************
    * \brief Assignment operator.
    * \param [in] rhs  HyperSphere instance on right-hand-side.
    ****************************************************************************
    */
   HyperSphere<T,NDIMS>& operator=(const HyperSphere<T,NDIMS>& rhs);

   /*!
    ****************************************************************************
    * \brief Returns the radius of the HyperSphere.
    * \return r the radius of the HyperSphere.
    ****************************************************************************
    */
   T radius() const { return m_radius; };

   /*!
    ****************************************************************************
    * \brief Returns the center of the HyperSphere
    * \return c pointer to array that holds the center of the HyperSphere.
    * \note The length of the array is NDIMS.
    ****************************************************************************
    */
   T* center() { return m_center; };
   const T* center() const { return m_center; };

   /*!
    ***************************************************************************
    * \brief Computes signed distance of a point to the HyperSphere boundary.
    * \param [in] q pointer to user-supplied point q.
    * \return dist signed distance
    * \pre q != ATK_NULLPTR
    * \pre q must be a pointer to an array that is at least NDIMS long
    ***************************************************************************
    */
   T getSignedDistance( T* q);

   /*!
    ***************************************************************************
    * \brief Computes orientation of a point with respect to the HyperSphere.
    * \param [in] q pointer to user-supplied point q.
    * \return orient orientation of q with respect to the sphere.
    * \pre q != ATK_NULLPTR
    * \pre q must be a pointer to an array that is at least NDIMS long
    ***************************************************************************
    */
   int getOrientation( T* q );

private:

   T m_center[ NDIMS ];
   T m_radius;

};

/// \name Pre-defined HyperSpheres for convenience
/// @{

typedef HyperSphere<double,3> Sphere;
typedef HyperSphere<double,2> Circle;
/// @}

} /* namespace quest */

//------------------------------------------------------------------------------
// HyperSphere Implementation
//------------------------------------------------------------------------------
namespace quest {

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
HyperSphere< T,NDIMS >::HyperSphere( T radius ) : m_radius(radius)
{
   std::fill( m_center, m_center+NDIMS, 0.0 );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
HyperSphere< T,NDIMS >::HyperSphere( T* center, T radius ) : m_radius(radius)
{
   SLIC_ASSERT( center != ATK_NULLPTR );
   memcpy( m_center, center, NDIMS*sizeof(T) );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
HyperSphere< T,NDIMS >::~HyperSphere()
{

}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline HyperSphere< T,NDIMS >& HyperSphere< T,NDIMS >::operator=(
        const HyperSphere< T,NDIMS >& rhs)
{
   if ( this == &rhs ) {
     return *this;
   }

   memcpy( m_center, rhs.m_center, NDIMS*sizeof(T) );
   m_radius = rhs.m_radius;
   return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
T HyperSphere< T,NDIMS >::getSignedDistance( T* q )
{
   SLIC_ASSERT( q != ATK_NULLPTR );

   T d = 0.0;
   for ( int i=0; i < NDIMS; ++i ) {
      const T dx = q[ i ]-m_center[ i ];
      d += (dx*dx);
   }

   return( std::sqrt( d )-m_radius );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
int HyperSphere< T,NDIMS >::getOrientation( T* q )
{
   SLIC_ASSERT( q != ATK_NULLPTR );

   T signed_distance = this->getSignedDistance( q );

   int orient = -1;

   if ( asctoolkit::utilities::isNearlyEqual( signed_distance, 0.0, 1.0e-9) ) {

       orient = ON_BOUNDARY;

   } else if ( signed_distance < 0.0f ) {

       // inside
       orient = ON_NEGATIVE_SIDE;

   } else {

       // outside
       orient = ON_POSITIVE_SIDE;

   }

   return orient;
}

}
#endif /* HYPERSPHERE_HPP_ */
