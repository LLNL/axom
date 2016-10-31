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
 * \file
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the "signed" squared distance between two geometric entities.
 *******************************************************************************
 */

#ifndef SQUAREDDISTANCE_HPP_
#define SQUAREDDISTANCE_HPP_

#include "quest/BoundingBox.hpp"
#include "quest/Point.hpp"
#include "quest/Segment.hpp"
#include "quest/Triangle.hpp"
#include "quest/Vector.hpp"
#include "quest/closest_point.hpp"

#include "slic/slic.hpp"


namespace quest {


/*!
 *******************************************************************************
 * \brief Computes the squared distance from point A to point B.
 * \param [in] A source point
 * \param [in] B end point.
 * \return d the distance from point A to point B.
 *******************************************************************************
 */
template < typename T, int ndims >
inline
double squared_distance( const Point< T,ndims >& A, const Point< T,ndims >& B )
{
  Vector< T,ndims > v( A, B );
  return( v.squared_norm() );
}

/*!
 *******************************************************************************
 * \brief Computes the minimum squared distance from a query point, P, to a
 *  given axis-aligned bounding box B.
 * \param [in] P the query point.
 * \param [in] B the axis-aligned bounding box.
 * \return d the signed distance from P to the closest point on B.
 *******************************************************************************
 */
template < typename T, int ndims >
inline
double squared_distance( const Point< T,ndims >& P,
                         const BoundingBox< T,ndims >& B )
{
   if ( B.contains( P ) ) {
      /* short-circuit */
      return 0.0f;
   }

   // compute closest point to the box
   Point< T,ndims > cp;
   for ( int i=0; i < ndims; ++i ) {

       cp[ i ] = P[ i ];
       if ( cp[i] < B.getMin()[i] ) {
           cp[ i ] = B.getMin()[i];
       }

       if ( cp[i] > B.getMax()[i] ) {
           cp[ i ] = B.getMax()[i];
       }
   }

   // return squared signed distance to the closest point
   return quest::squared_distance( P, cp );
}

/*!
 *******************************************************************************
 * \brief Computes the minimum squared distance from a query point, P, to the
 *  given segment, S.
 * \param [in] P the query point.
 * \param [in] S the input segment.
 * \return d the minimum distance from P on the
 *******************************************************************************
 */
template < typename T, int ndims >
inline
double squared_distance( const Point< T,ndims >& P,
                         const Segment< T,ndims >& S )
{
   Vector< T,ndims > ab( S.source(), S.target() );
   Vector< T,ndims > ac( S.source(), P );

   const T e = Vector< T,ndims >::dot_product( ac, ab );

   // outside segment, on the side of a
   // Testing if closest point is A
   if ( e <= 0.0f ) {

      return ac.squared_norm();

   }

   // outside segment, on the side of b
   // Testing if closest point is B
   const T f = ab.squared_norm();
   if ( e >= f ) {

      Vector< T,ndims > bc( S.target(), P );
      return bc.squared_norm();

   }

   // P projects onto the segment
   // Otherwise, we are in between A,B, therefore we project inside A,B.
   const T dist = ac.squared_norm() - ( e*e/f );
   return dist;
}

/*!
 *******************************************************************************
 * \brief Computes the minimum squared distance from a query point, P, to the
 *  closest point on the given triangle.
 * \param [in] P the query point.
 * \param [in] tri the supplied triangle.
 * \return d the distance from Q to the closest point on the triangle T.
 *******************************************************************************
 */
template < typename T, int ndims >
inline
double squared_distance( const Point< T,ndims >& P,
                         const Triangle< T,ndims >& tri )
{
  Point< T,ndims > cpt = closest_point( P, tri );
  return squared_distance( P, cpt );
}


}


#endif /* SQUAREDDISTANCE_HPP_ */
