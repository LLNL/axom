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

/*!
 *
 * \file vector_utilities.hpp
 *
 * \brief Header file containing utility functions to be used for vector
 *  calculations where "vectors" are arrays.
 *
 */

#ifndef AXOM_NUMERICS_VECTOR_UTILITIES_HPP_
#define AXOM_NUMERICS_VECTOR_UTILITIES_HPP_

#include "axom/Types.hpp"                  // for AXOM_NULLPTR

#include "axom_utils/Determinants.hpp"     // for numerics::determinant()
#include "axom_utils/Utilities.hpp"        // for isNearlyEqual()

// C/C++ includes
#include <cassert> // for assert()
#include <cmath>   // for sqrt()

namespace axom
{
namespace numerics
{

/*!
 * \brief Generates a vector consisting of a sequence of N uniformly spaced
 *  values over the interval [x0,x1].
 *
 * \param [in] x0 scalar, the start value of the sequence.
 * \param [in] x1 scalar, the end value of the sequence.
 * \param [out] v user-supplied buffer where to store the sequence of numbers
 * \param [in] N the size of the computed vector sequence.
 *
 * \return status true if successful, otherwise, false.
 *
 * \note The output vector, v, must be able to hold at least N values.
 * \note if x0 < x1, the sequence of values will be in ascending order.
 * \note if x0 > x1, the sequence of values will be in descending order.
 *
 * \pre v != AXOM_NULLPTR
 * \pre N > 1
 *
 */
template < typename T >
inline bool linspace( const T& x0, const T& x1, T* v, int N );

/*!
 * \brief Computes the vector cross-product of two vectors, u and v.
 *
 * \param [in] u array pointer to the vector u
 * \param [in] v array pointer to the vector v
 * \param [out] w array pointer where to store the cross-product
 *
 * \note The u, v, and w are expected to be 3D vectors and are expected to
 *  be pointing to array of at least size 3.
 *
 * \pre u != AXOM_NULLPTR
 * \pre v != AXOM_NULLPTR
 * \pre w != AXOM_NULLPTR
 *
 */
template < typename T >
inline void cross_product( const T* u, const T* v, T* w );

/*!
 * \brief Computes the dot product of the arrays u and v.
 *
 * \tparam T data type
 * \param [in] u pointer to an array of size dim
 * \param [in] v pointer to an array of size dim
 * \param [in] dim the dimension of the arrays at u and v
 * \return the dot product of the arrays u and v
 *
 * \pre dim >= 1
 * \pre u != AXOM_NULLPTR
 * \pre v != AXOM_NULLPTR
 * \pre u has at least dim entries
 * \pre v has at least dim entries
 */
template < typename T >
inline T dot_product( const T* u, const T* v, int dim);

/*!
 * \brief Makes u orthogonal to v.
 *
 * \tparam T data type
 * \param [in, out] u vector to be made orthogonal to other; saves in-place
 * \param [in] v vector that u is made orthogonal to
 * \param [in] dim dimension of vectors
 * \param [in] tol tolerance; if the norm of v is less than tol we do nothing
 *
 * \pre dim >= 1
 * \pre u != AXOM_NULLPTR
 * \pre v != AXOM_NULLPTR
 * \pre T is a floating point type
 */
template < typename T >
void make_orthogonal(T* u, T* v, int dim, double tol=1E-16);

/*!
 * \brief Performs Gram-Schmidt orthonormalization in-place on a 2D array
 * of shape size,dim where it treats the rows as the individual vectors.
 *
 * \tparam T data type
 * \param [in, out] basis vectors to be made orthonormal; saves them in-place
 * \param [in] size number of vectors
 * \param [in] dim dimension of vectors
 * \param [in] eps If one of the vectors, after being made orthogonal to the
 *  others, has norm less than eps, then the orthonormalization is declared a
 *  failure. Note that this may well destory the data pointed to by basis.
 * \return true if the orthonormalization is successful, false otherwise
 *
 * \pre dim >= 1
 * \pre 1 <= size <= dim
 * \pre basis != AXOM_NULLPTR
 * \pre T is a floating point type
 */
template < typename T >
bool orthonormalize(T* basis, int size, int dim, double eps = 1E-16);

/*!
 * \brief Normalizes the passed in array.
 *
 * \tparam T data type
 * \param [in] v pointer the array
 * \param [in] dim the dimension of v
 * \param [in] eps fuzz parameter
 * \note If squared norm of v is less than eps, the normalization
 *  fails and we do nothing to v
 *
 * \return true if normalization is successful, false otherwise.
 *
 * \pre dim >= 1
 * \pre v != AXOM_NULLPTR
 * \pre T is a floating point type
 */
template < typename T >
inline bool normalize(T* v, int dim, double eps = 1e-16);

} /* end namespace numerics */
} /* end namespace axom */


//------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace numerics
{

//------------------------------------------------------------------------------
template < typename T >
inline bool linspace( const T& x0, const T& x1, T* v, int N )
{
  AXOM_STATIC_ASSERT_MSG( std::is_floating_point< T >::value,
                             "pre: T is a floating point type" );
  assert( "pre: v pointer is null" && (v != AXOM_NULLPTR) );

  if ( N <= 1 )
  {
    /* short-circuit */
    return false;
  }

  const T h = (x1-x0) / static_cast< T >( N-1 );

  for ( int i=0; i < N; ++i )
  {
    v[ i ] = x0 + i*h;
  }

  return true;
}

//------------------------------------------------------------------------------
template < typename T >
inline void cross_product( const T* u, const T* v, T* w )
{
  assert( "pre: u pointer is null" && (u != AXOM_NULLPTR) );
  assert( "pre: v pointer is null" && (v != AXOM_NULLPTR) );
  assert( "pre: w pointer is null" && (w != AXOM_NULLPTR) );

  w[ 0 ] = numerics::determinant( u[1], u[2], v[1], v[2] );

  // NOTE: transpose u,v to negate
  w[ 1 ] = numerics::determinant( v[0], v[2], u[0], u[2] );
  w[ 2 ] = numerics::determinant( u[0], u[1], v[0], v[1] );
}

//------------------------------------------------------------------------------
template < typename T >
inline T dot_product( const T* u, const T* v, int dim)
{
  assert("pre: u pointer is null" && (u != AXOM_NULLPTR));
  assert("pre: v pointer is null" && (v != AXOM_NULLPTR));
  assert("pre: dim >= 1" && (dim >= 1));

  T res = u[0]*v[0];
  for (int i = 1 ; i < dim ; ++i)
  {
    res += u[i]*v[i];
  }

  return res;
}

//------------------------------------------------------------------------------
template < typename T >
void make_orthogonal(T* u, T* v, int dim, double tol)
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point< T >::value,
                         "pre: T is a floating point type");
  assert("pre: u pointer is null" && (u != AXOM_NULLPTR));
  assert("pre: v pointer is null" && (v != AXOM_NULLPTR));
  assert("pre: dim >= 1" && (dim >= 1));

  double norm = static_cast< double >(dot_product(v, v, dim));

  if (norm < tol)
    return;

  T tnorm = static_cast< T >(norm);

  T dot = dot_product(u, v, dim);

  for (int l = 0 ; l < dim ; ++l)
    u[l] -= ((dot*v[l])/tnorm);
}

//------------------------------------------------------------------------------
template < typename T >
bool orthonormalize(T* basis, int size, int dim, double eps)
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point< T >::value,
                         "pre: T is a floating point type");
  assert("pre: basis pointer is null" && (basis != AXOM_NULLPTR));
  assert("pre: dim >= 1" && (dim >= 1));
  assert("pre: size >= 1" && (size >= 1));
  assert("pre: size <= dim" && (size <= dim));

  for (int i = 0 ; i < size ; ++i)
  {
    T* curr = &basis[i*dim];

    // make curr orthogonal to previous ones
    for (int j = 0 ; j < i ; ++j)
    {
      T* other = &basis[j*dim];

      make_orthogonal(curr, other, dim);
    }

    bool res = normalize(curr, dim, eps);

    if (!res)
    {
      return false;
    }
  }

  // success
  return true;
}

//------------------------------------------------------------------------------
template < typename T >
inline bool normalize(T* v, int dim, double eps)
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point< T >::value,
                         "pre: T is a floating point type");
  assert("pre: v pointer is null" && (v != AXOM_NULLPTR));
  assert("pre: dim >= 1" && (dim >= 1));

  const double norm = static_cast< double >( dot_product(v, v, dim) );

  if (utilities::isNearlyEqual< double >(norm, 0., eps))
  {
    return false;
  }

  const T tnorm = 1.0 / static_cast< T >(std::sqrt(norm));
  for (int l = 0 ; l < dim ; ++l)
  {
    v[l] *= tnorm;
  }

  // success
  return true;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif /* AXOM_NUMERICS_VECTOR_UTILITIES_HPP_ */
