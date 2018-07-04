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
 * \file matvecops.hpp
 *
 * \brief Provides Matrix/Vector operations.
 *
 */

#ifndef AXOM_NUMERICS_MATVECOPS_HPP_
#define AXOM_NUMERICS_MATVECOPS_HPP_

#include "axom/Types.hpp"                  // for AXOM_NULLPTR

#include "axom_utils/Determinants.hpp"     // for numerics::determinant()
#include "axom_utils/Matrix.hpp"           // for numerics::Matrix
#include "axom_utils/Utilities.hpp"        // for isNearlyEqual()

// C/C++ includes
#include <cassert> // for assert()
#include <cmath>   // for sqrt()

namespace axom
{
namespace numerics
{

/// \name Vector Operators
/// @{

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

/// @}

/// \name Matrix Operators
/// @{

/*!
 * \brief Computes \f$ \mathcal{C} = \mathcal{A} + \mathcal{B} \f$
 *
 * \param [in]  A \f$ M \times N \f$ matrix on the left-hand side.
 * \param [in]  B \f$ M \times N \f$ matrix on the right-hand side.
 * \param [out] C \f$ M \times N \f$ output matrix.
 *
 * \post \f$ c_{ij} = a_{ij} + b{ij} \f$, \f$\forall c_{ij} \in \mathcal{C}\f$
 *
 * \note Matrix addition is undefined for matrices that have different
 *  dimensions. If the input matrices, \f$ \mathcal{A} \f$, \f$ \mathcal{B} \f$
 *  or \f$ \mathcal{C} \f$ have different dimensions, a \f$ 1 \times 1 \f$ null
 *  matrix is returned in \f$ \mathcal{C} \f$
 *
 * \return status true if addition is successful, otherwise, false.
 */
template < typename T >
inline bool matrix_add( const Matrix< T >& A,
                        const Matrix< T >& B,
                        Matrix< T >& C );

/*!
 * \brief Computes \f$ \mathcal{C} = \mathcal{A} - \mathcal{B} \f$
 *
 * \param [in]  A \f$ M \times N \f$ matrix on the left-hand side.
 * \param [in]  B \f$ M \times N \f$ matrix on the right-hand side.
 * \param [out] C \f$ M \times N \f$ output matrix.
 *
 * \post \f$ c_{ij} = a_{ij} + b{ij} \f$, \f$\forall c_{ij} \in \mathcal{C}\f$
 *
 * \note Matrix subtraction is undefined for matrices that have different
 *  dimensions. If the input matrices, \f$ \mathcal{A} \f$, \f$ \mathcal{B} \f$
 *  or \f$ \mathcal{C} \f$ have different dimensions, a \f$ 1 \times 1 \f$ null
 *  matrix is returned in \f$ \mathcal{C} \f$
 *
 * \return status true if addition is successful, otherwise, false.
 */
template < typename T >
inline bool matrix_subtract( const Matrix< T >& A,
                             const Matrix< T >& B,
                             Matrix< T >& C );
/*!
 * \brief Computes the matrix-matrix product of \f$ \mathcal{A} \f$ and
 *  \f$ \mathcal{B} \f$ and stores the result in \f$ \mathcal{C} \f$
 *
 * \param [in] A  \f$ M \times N \f$ matrix on the left hand-side.
 * \param [in] B  \f$ N \times K \f$ matrix on the right hand-side.
 * \param [out] C \f$ M \times K \f$ output matrix
 *
 * \pre The inner dimensions of matrices A, B must match
 * \pre Output matrix should be an \f$ M \times K \f$ matrix
 *
 * \note Matrix multiplication is undefined for matrices with different inner
 *  dimension. If the inner dimensions are not matching, the code returns a
 *  \f$ 1 \times 1 \f$ null matrix in \f$ \mathcal{C} \f$
 *
 * * \note Matrix \f$ \mathcal{C} \f$ should be a different Matrix instance than
 *  \f$ \mathcal{A} \f$ , \f$ \mathcal{B} \f$
 *
 * \return status true if the multiplication is successful, otherwise, false.
 */
template < typename T >
inline bool matrix_multiply( const Matrix< T >&A,
                             const Matrix< T >& B,
                             Matrix< T >& C );

/*!
 * \brief Computes a scalar-matrix produect of a matrix \f$ \mathcal{A} \f$ and
 *  a scalar \f$ c \f$ and stores the result in \f$ \mathcal{A} \f$
 *
 * \param [in,out] A an \f$ M \times N \f$ matrix
 * \param [in] c scalar value to multiply the matrix coefficients
 *
 * \post \f$ a_{ij} = c \cdot a_{ij} \f$, \f$ \forall a_{ij} \in mathcal{A} \f$
 */
template < typename T >
inline void matrix_scalar_multiply( Matrix< T >& A, const T& c );

/*!
 * \brief Computes the matrix-vector product of a matrix \f$\mathcal{A}\f$ and
 *  a vector \f$\mathbf{x}\f$ and store the result in the user-supplied output
 *  vector.
 *
 * \param [in] A an \f$ M \times N \f$ matrix
 * \param [in] vec pointer to user-supplied vector storing the vector
 * \param [out] output pointer to the user-supplied output vector
 *
 * \pre vec != AXOM_NULLPTR
 * \pre vec must be of dimension \f$ N \f$
 * \pre output != AXOM_NULLPTR
 * \pre output must be of dimension \f$ M \f$
 * \pre vec != output
 *
 * \post \f$ b_i = \sum\limits_{j=0}^N a_{ij} \cdot x_j \f$,
 *  \f$\forall i \in [0,M-1] \f$
 */
template < typename T >
inline void matrix_vector_multiply( const Matrix< T >& A,
                                    const T* vec,
                                    T* output );

/*!
 * \brief Computes the matrix transpose of a given matrix \f$ \mathcal{A} \f$
 *
 * \param [in] A an \f$ m \times n \f$ matrix
 * \param [out] M an \f$ n \times m \f$ matrix where to store the transpose.
 *
 * \note If the supplied matrix does not have the correct dimensions, the code
 *  returns a \f$ 1 \times 1 \f$ null matrix in \f$ \mathcal{M} \f$
 *
 * \note Matrix \f$ \mathcal{M} \f$ should be a different Matrix instance than
 *  \f$ \mathcal{A} \f$
 *
 * \return status true if the matrix transpose is successful, otherwise, false.
 */
template < typename T >
inline bool matrix_transpose( const Matrix< T >& A, Matrix< T >& M );

/// @}

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
// IMPLEMENTATION OF MATRIX OPERATIONS
//------------------------------------------------------------------------------
template < typename T >
inline bool matrix_add( const Matrix< T >& A,
                        const Matrix< T >& B,
                        Matrix< T >& C )
{
  if ( A.getNumRows()    != B.getNumRows() ||
       A.getNumColumns() != B.getNumColumns() ||
       C.getNumRows()    != B.getNumRows() ||
       C.getNumColumns() != B.getNumColumns() )
  {

    // matrix dimensions do not match
    C = Matrix< T >::zeros(1,1);
    return false;
  }

  typedef typename Matrix< T >::IndexType IndexType;

  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();
  const int N     = nrows*ncols;

  T* target = C.data();

  const T* sourceA = A.data();
  const T* sourceB = B.data();
  for ( IndexType i=0 ; i < N ; ++i )
  {
    target[ i ] = sourceA[ i ] + sourceB[ i ];
  }

  return true;
}

//------------------------------------------------------------------------------
template < typename T >
bool matrix_subtract( const Matrix< T >& A,
                      const Matrix< T >& B,
                      Matrix< T >& C )
{
  if ( A.getNumRows()    != B.getNumRows() ||
       A.getNumColumns() != B.getNumColumns() ||
       C.getNumRows()    != B.getNumRows() ||
       C.getNumColumns() != B.getNumColumns() )
  {

    // matrix dimensions do not match
    C = Matrix< T >::zeros(1,1);
    return false;
  }

  typedef typename Matrix< T >::IndexType IndexType;

  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();
  const int N     = nrows*ncols;

  T* target = C.data();

  const T* sourceB = B.data();
  const T* sourceA = A.data();
  for ( IndexType i=0 ; i < N ; ++i )
  {
    target[ i ] = sourceA[ i ] - sourceB[ i ];
  }

  return true;
}

//------------------------------------------------------------------------------
template < typename T >
inline bool matrix_multiply( const Matrix< T >& A,
                             const Matrix< T >& B,
                             Matrix< T >& C )
{
  if ( A.getNumColumns() != B.getNumRows()   ||
       C.getNumRows() != A.getNumRows()      ||
       C.getNumColumns() != B.getNumColumns()   )
  {

    C = Matrix< T >::zeros( 1,1 );
    return false;
  }

  typedef typename Matrix< T >::IndexType IndexType;

  const int nk    = A.getNumColumns();
  const int nrows = C.getNumRows();
  const int ncols = C.getNumColumns();

  for ( IndexType i=0 ; i < nrows ; ++i )
  {
    for ( IndexType j=0 ; j < ncols ; ++j )
    {

      C( i,j ) = static_cast< T >( 0 );
      for ( IndexType k=0 ; k < nk ; ++k )
      {
        C( i,j ) += A(i,k) * B(k,j);
      }   // END for all internal

    }  // END for all columns
  } // END for all rows

  return true;
}

//------------------------------------------------------------------------------
template < typename T >
inline void matrix_scalar_multiply( Matrix< T >& A, const T& c )
{
  // matrix-scalar multiplication
  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();
  const int N     = nrows*ncols;

  typedef typename Matrix< T >::IndexType IndexType;

  T* target = A.data();
  for ( IndexType i=0 ; i < N ; ++i )
  {
    target[ i ] *= c;
  }

}

//------------------------------------------------------------------------------
template < typename T >
inline void matrix_vector_multiply( const Matrix< T >& A, const T* x, T* b )
{
  // matrix-vector multiplication
  const int nrows = A.getNumRows();
  const int ncols = A.getNumColumns();

  typedef typename Matrix< T >::IndexType IndexType;

  for ( IndexType i=0 ; i < nrows ; ++i )
  {

    double sum = 0.0;
    for ( IndexType j=0 ; j < ncols ; ++j )
    {
      sum += A( i,j )*x[ j ];
    } // END for all columns

    b[ i ] = static_cast< T >( sum );

  } // END for all rows

}

//-----------------------------------------------------------------------------
template < typename T >
inline bool matrix_transpose( const Matrix< T >& A, Matrix< T >& M )
{
  if ( M.getNumRows()    != A.getNumColumns() ||
       M.getNumColumns() != A.getNumRows()  )
  {
    M = Matrix< T >::zeros(1,1);
    return false;
  }

  const int ncols = A.getNumColumns();
  const int nrows = A.getNumRows();
  using IndexType = typename Matrix< T >::IndexType;

  for ( IndexType i=0 ; i < nrows ; ++i )
  {
    for ( IndexType j=0 ; j < ncols ; ++j )
    {
      M( j,i ) = A( i,j );
    }  // END for all columns in A
  } // END for all rows in A

  return true;
}


//------------------------------------------------------------------------------
// IMPLEMENTATION OF VECTOR OPERATIONS
//------------------------------------------------------------------------------

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

  for ( int i=0 ; i < N ; ++i )
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
