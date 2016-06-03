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
 * \brief Consists of a set of methods that compute the closest point on a
 *  geometric primitive B from another geometric primitive A.
 *
 *******************************************************************************
 */

#ifndef CLOSEST_POINT_HPP_
#define CLOSEST_POINT_HPP_

#include "quest/Point.hpp"
#include "quest/Triangle.hpp"

namespace quest
{

/*!
 *******************************************************************************
 * \brief Computes the closest point from a point, P, to a given triangle.
 *
 * \param [in] P the query point
 * \param [in] tri user-supplied triangle.
 *
 * \return cp the closest point from a point P and a triangle.
 *
 * \pre NDIMS==2 || NDIMS==3
 *
 * \note Implementation is based on "Real Time Collision Detection,
 *  Chapter 5.1.5 Closest Point on Triangle to Point".
 *******************************************************************************
 */
template < typename T, int NDIMS >
inline Point< T,NDIMS > closest_point( const Point< T,NDIMS >& P,
                                       const Triangle< T,NDIMS >& tri )
{
  // Check if P in vertex region outside A
  Vector< T, NDIMS > ab( tri.A(),tri.B() );
  Vector< T, NDIMS > ac( tri.A(),tri.C() );
  Vector< T, NDIMS > ap( tri.A(),P );
  T d1 = Vector< T,NDIMS >::dot_product( ab, ap );
  T d2 = Vector< T,NDIMS >::dot_product( ac, ap );
  if ( d1 <= 0.0f && d2 <= 0.0f ) {

      // A is the closest point
      return ( tri.A() );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside B
  Vector< T,NDIMS > bp( tri.B(), P );
  T d3 = Vector< T,NDIMS >::dot_product( ab, bp );
  T d4 = Vector< T,NDIMS >::dot_product( ac, bp );
  if ( d3 >= 0.0f && d4 <= d3 ) {

      // B is the closest point
      return ( tri.B() );

  } // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of AB
  T vc = d1*d4 - d3*d2;
  if ( vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f ) {

    T v = d1 / ( d1-d3 );
    Vector< T,NDIMS > v_ab = ab*v;

    double x = tri.A()[0] + v_ab[0];
    double y = tri.A()[1] + v_ab[1];
    double z = (NDIMS==3)? tri.A()[2] + v_ab[2] : 0.0;

    return ( Point<T,NDIMS>::make_point( x,y,z ) );
  } // END if

  //----------------------------------------------------------------------------
  // Check if P in vertex region outside C
  Vector< T,NDIMS > cp( tri.C(), P );
  T d5 = Vector< T,NDIMS >::dot_product(ab,cp);
  T d6 = Vector< T,NDIMS >::dot_product(ac,cp);
  if ( d6 >= 0.0f && d5 <= d6 ) {

     // C is the closest point
     return ( tri.C() );
  }

  //----------------------------------------------------------------------------
  // Check if P in edge region of AC
  T vb = d5*d2 - d1*d6;
  if ( vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f ) {

    T w = d2 / (d2-d6);
    Vector< T, NDIMS > w_ac = ac*w;

    double x = tri.A()[0] + w_ac[0];
    double y = tri.A()[1] + w_ac[1];
    double z = (NDIMS==3)? tri.A()[2] + w_ac[2] : 0.0;

    return ( Point< T,NDIMS >::make_point( x,y,z ) );
  } // END if

  //----------------------------------------------------------------------------
  // Check if P in edge region of BC
  T va = d3*d6 - d5*d4;
  if ( va <= 0.0f && (d4-d3) >= 0.0f && (d5-d6) >= 0.0f ) {

    T w = (d4-d3)/( (d4-d3)+(d5-d6) );
    Vector< T,NDIMS > bc( tri.B(), tri.C() );
    Vector< T,NDIMS > w_bc = bc*w;

    double x = tri.B()[0] + w_bc[0];
    double y = tri.B()[1] + w_bc[1];
    double z = (NDIMS==3)? tri.B()[2] + w_bc[2] : 0.0;

    return ( Point< T,NDIMS >::make_point( x,y,z ) );
  } // END if

  //----------------------------------------------------------------------------
  // P is inside face region
  T denom = 1.0f / (va + vb + vc );
  T v     = vb * denom;
  T w     = vc * denom;
  Vector< T,NDIMS > N = (ab*v) + (ac*w);

  double x = tri.A()[0] + N[0];
  double y = tri.A()[1] + N[1];
  double z = (NDIMS==3)? tri.A()[2] + N[2] : 0.0;

  return ( Point< T,NDIMS >::make_point( x,y,z ) );
}

}



#endif /* CLOSEST_POINT_HPP_ */
