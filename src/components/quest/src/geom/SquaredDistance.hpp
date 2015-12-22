/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the "signed" squared distance between two geometric entities.
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef SQUAREDDISTANCE_HPP_
#define SQUAREDDISTANCE_HPP_

#include "quest/BoundingBox.hpp"
#include "quest/Point.hpp"
#include "quest/Triangle.hpp"
#include "quest/Vector.hpp"


#include "slic/slic.hpp"


namespace quest
{


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
 *  closest point on the given triangle.
 * \param [in] P the query point.
 * \param [in] tri the supplied triangle.
 * \return d the distance from Q to the closest point on the triangle T.
 * \note The algorithm is based on "Real Time Collision Detection, Chapter 5.1.5
 *  Closest Point on Triangle to Point".
 *******************************************************************************
 */
template < typename T, int ndims >
inline
double squared_distance( const Point< T,ndims >& P,
                         const Triangle< T,ndims >& tri )
{
  // Check if P in vertex region outside A
  Vector< T, ndims > ab( tri.A(),tri.B() );
  Vector< T, ndims > ac( tri.A(),tri.C() );
  Vector< T, ndims > ap( tri.A(),P );
  T d1 = Vector< T,ndims >::dot_product( ab, ap );
  T d2 = Vector< T,ndims >::dot_product( ac, ap );
  if ( d1 <= 0.0f && d2 <= 0.0f ) {

      // A is the closest point
      return quest::squared_distance( P, tri.A() );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside B
  Vector< T,ndims > bp( tri.B(), P );
  T d3 = Vector< T,ndims >::dot_product( ab, bp );
  T d4 = Vector< T,ndims >::dot_product( ac, bp );
  if ( d3 >= 0.0f && d4 <= d3 ) {

      // B is the closest point
      return quest::squared_distance( P, tri.B() );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of AB
  T vc = d1*d4 - d3*d2;
  if ( vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f ) {

    T v = d1 / ( d1-d3 );
    Vector< T,ndims > v_ab = ab*v;

    double x = tri.A()[0] + v_ab[0];
    double y = tri.A()[1] + v_ab[1];
    double z = (ndims==3)? tri.A()[2] + v_ab[2] : 0.0;

    Point< T,ndims > pt = Point<T,ndims>::make_point( x,y,z );
    return quest::squared_distance( P, pt );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside C
  Vector< T,ndims > cp( tri.C(), P );
  T d5 = Vector< T,ndims >::dot_product(ab,cp);
  T d6 = Vector< T,ndims >::dot_product(ac,cp);
  if ( d6 >= 0.0f && d5 <= d6 ) {

     // C is the closest point
     return quest::squared_distance( P, tri.C() );
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of AC
  T vb = d5*d2 - d1*d6;
  if ( vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f ) {

    T w = d2 / (d2-d6);
    Vector< T, ndims > w_ac = ac*w;

    double x = tri.A()[0] + w_ac[0];
    double y = tri.A()[1] + w_ac[1];
    double z = (ndims==3)? tri.A()[2] + w_ac[2] : 0.0;

    Point< T,ndims > pt = Point< T,ndims >::make_point( x,y,z );
    return quest::squared_distance( P, pt );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of BC
  T va = d3*d6 - d5*d4;
  if ( va <= 0.0f && (d4-d3) >= 0.0f && (d5-d6) >= 0.0f ) {

    T w = (d4-d3)/( (d4-d3)+(d5-d6) );
    Vector< T,ndims > bc( tri.B(), tri.C() );
    Vector< T,ndims > w_bc = bc*w;

    double x = tri.B()[0] + w_bc[0];
    double y = tri.B()[1] + w_bc[1];
    double z = (ndims==3)? tri.B()[2] + w_bc[2] : 0.0;

    Point< T,ndims > pt = Point< T,ndims >::make_point( x,y,z );
    return quest::squared_distance( P, pt );

  } // END if

  //----------------------------------------------------------------------------
  // P is inside face region
  T denom = 1.0f / (va + vb + vc );
  T v     = vb * denom;
  T w     = vc * denom;
  Vector< T,ndims > N = (ab*v) + (ac*w);

  double x = tri.A()[0] + N[0];
  double y = tri.A()[1] + N[1];
  double z = (ndims==3)? tri.A()[2] + N[2] : 0.0;

  Point< T,ndims > pt = Point< T,ndims >::make_point( x,y,z );

  return quest::squared_distance( P, pt );
}


}


#endif /* SQUAREDDISTANCE_HPP_ */
