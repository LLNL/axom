/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Sep 4, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef VECTOR_HXX_
#define VECTOR_HXX_

#include "quest/Point.hpp"
#include "quest/Determinants.hpp" // For math::determinant()

// C/C++ includes
#include <cmath>

namespace quest {

/*!
 *******************************************************************************
 * \class Vector
 *
 * \brief Represents a vector, \f$ v \in \mathcal{R}^d. It provides access
 *  methods for setting and querying the vector components as well as vector
 *  math operators, e.g., adding, subtracting, dot_product and cross_product.
 *
 * \see Point
 *******************************************************************************
 */
template < typename T, int ndims >
class Vector : public Point< T, ndims >
{
public:

  /*!
   *****************************************************************************
   * \brief Default constructor
   *****************************************************************************
   */
  Vector();

  /*!
   *****************************************************************************
   * \brief Constructs a vector from point A to point B.
   * \param [in] A origin point of the vector.
   * \param [in] B destination point of the vector.
   * \pre A.dimension() == B.dimension()
   * \pre A.dimension() == ndims
   *****************************************************************************
   */
  Vector( const Point< T,ndims >& A, const Point< T,ndims >& B );

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
   ~Vector() {}

  /*!
   *****************************************************************************
   * \brief Multiplies the scalar with this Vector instance.
   * \param [in] scalar the scalar value to multiply with this vector.
   * \return c the resulting vector after the multiplication.
   *****************************************************************************
   */
  Vector< T,ndims > operator*(T scalar);

  /*!
   *****************************************************************************
   * \brief Computes the vector sum \f$\vec{r}=\vec{u}+\vec{v}\f$
   * \param [in] v the vector on the right hand side.
   * \return r the resulting vector after the sum operation.
   *****************************************************************************
   */
  Vector< T,ndims > operator+( const Vector<T,ndims>& v );

  /*!
   *****************************************************************************
   * \brief Computes the vector subtraction \f$\vec{r}=\vec{u}+\vec{v}\f$
   * \param [in] v the vector on the right hand side.
   * \return r the resulting vector after the subtraction operation.
   *****************************************************************************
   */
  Vector< T,ndims > operator-( const Vector<T,ndims>& v );

  /*!
   *****************************************************************************
   * \brief Computes the squared \f$ l^2\f$ norm of this vector instance.
   * \return n the squared norm.
   * \see Vector::norm()
   *****************************************************************************
   */
  double squared_norm() const;

  /*!
   *****************************************************************************
   * \brief Computes the \f$ l^2 \f$ norm of this vector instance.
   * \return n the norm of the vector, a.k.a., magnitude or length.
   *****************************************************************************
   */
  double norm() const;

  /*!
   *****************************************************************************
   * \brief Normalizes this vector instance.
   * \post this->norm() == 1.0f
   *****************************************************************************
   */
  void normalize();

  /*!
   *****************************************************************************
   * \brief Computes the dot product of two vectors u, v.
   * \param [in] u the vector on the right-hand side.
   * \param [in] v the vector on the left-hand side.
   * \return dotprod the computed dot product.
   * \pre u.dimension() == v.dimension().
   *****************************************************************************
   */
  static T dot_product( const Vector< T,ndims >& u,
                        const Vector< T,ndims >& v );

  /*!
   *****************************************************************************
   * \brief Computes the 3-D cross product of vector u and v.
   * \param [in] u the vector on the right-hand side.
   * \param [in] v the vector on the left-hand side.
   * \return C the resulting vector from A x B.
   *****************************************************************************
   */
  static Vector< T,3 > cross_product( const Vector< T,3 >& u,
                                      const Vector< T,3 >& v );

};

/// \name Pre-defined Vector types
/// @{

typedef Vector< double, 2 > Vector2D;
typedef Vector< double, 3 > Vector3D;

/// @}

} /* namespace  */

//------------------------------------------------------------------------------
//  Vector implementation
//------------------------------------------------------------------------------
namespace quest {

//------------------------------------------------------------------------------
template < typename T, int ndims >
Vector< T,ndims >::Vector()
{
  SLIC_ASSERT( ndims >= 1 );
  std::fill( this->m_components, this->m_components+ndims, 1.0 );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
Vector< T,ndims >::Vector( const Point< T,ndims >& A,
                           const Point< T,ndims >& B )
{
  SLIC_ASSERT( A.dimension() == B.dimension( ) );

  for ( int i=0; i < ndims; ++i ) {
      this->m_components[ i ] = B[ i ]-A[ i ];
  }

}


//------------------------------------------------------------------------------
template < typename T, int ndims >
inline Vector< T,ndims > Vector< T,ndims >::operator*( T scalar )
{
  Vector< T, ndims > result;
  for( int i=0; i < ndims; ++i ) {
    result[ i ] = this->m_components[i] * scalar;
  }
  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline Vector< T, ndims > Vector< T,ndims >::operator+(const Vector<T,ndims>& v)
{
  Vector< T, ndims > result;
  for ( int i=0; i < ndims; ++i ) {
    result[ i ] = this->m_components[ i ] + v[ i ];
  }

  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline Vector< T, ndims > Vector< T,ndims >::operator-(const Vector<T,ndims>& v)
{
  Vector< T, ndims > result;
  for ( int i=0; i < ndims; ++i ) {
    result[ i ] = this->m_components[ i ] - v[ i ];
  }

  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline double Vector< T, ndims >::squared_norm() const
{
  double sum = 0.0;
  for ( int i=0; i < ndims; ++i ) {
      sum += this->m_components[ i ] * this->m_components[ i ];
  }
  return( sum );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline double Vector< T, ndims >::norm() const
{
  return std::sqrt( squared_norm() );
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline void Vector< T, ndims >::normalize()
{
  const double length = this->norm();
  for ( int i=0; i < ndims; ++i ) {
      this->m_components[ i ] /= length;
  }
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline T Vector< T,ndims >::dot_product( const Vector< T,ndims >& u,
                                       const Vector< T,ndims >& v )
{
  SLIC_ASSERT( u.dimension() == v.dimension() );
  double dotprod = 0.0;

  for ( int i=0; i < ndims; ++i ) {
     dotprod += u[ i ]*v[ i ];
  } // END for all ndims

  return dotprod;
}

//------------------------------------------------------------------------------
template < typename T, int ndims >
inline Vector< T,3 > Vector< T,ndims >::cross_product(
        const Vector< T,3 >& u, const Vector< T,3 >& v )
{
  Vector< T,3 > c;
  c[ 0 ] = math::determinant( u[1],u[2], v[1],v[2] );
  c[ 1 ] = math::determinant( u[0],u[2], v[0],v[2] ) * ( -1.0f );
  c[ 2 ] = math::determinant( u[0],u[1], v[0],v[1] );
  return( c );
}


}
#endif /* VECTOR_HXX_ */
