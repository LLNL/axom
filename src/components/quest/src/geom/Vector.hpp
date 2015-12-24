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
template < typename T, int DIM >
class Vector : public Point< T, DIM >
{
public:
    typedef Point<T,DIM> PointType;

public:

  /*!
   *****************************************************************************
   * \brief Fill vector with single value. Acts as default constructor
   * Sets first sz components of the vector to val (default 0).
   * \param [in] val The value to set the coordinates to.  Defaults to zero
   * \param [in] sz The number of coordinates to set to val.
   * The rest will be set to zero.  Defaults is DIM.
   *****************************************************************************
   */
  explicit Vector(T val = T(), int sz = DIM) : PointType(val, sz) {}

  /*!
   *****************************************************************************
   * \brief Creates a vector from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz The number of coordinates to take from the array.  Defaults
   * to DIM.
   * It sz is greater than DIM, we only take the first DIM values.
   *****************************************************************************
   */
  Vector(T* vals, int sz = DIM) : PointType(vals, sz) {}

  /*!
   *****************************************************************************
   * \brief Constructor to create vector from a Point
   * \param [in] pt The point containing the vector's coordinates.
   * \note Equivalent to Vector( Point::zero(), pt)
   *****************************************************************************
   */
  Vector(const Point<T,DIM> & pt) : PointType(pt) {}

  /*!
   *****************************************************************************
   * \brief Constructs a vector from point A to point B.
   * \param [in] A origin point of the vector.
   * \param [in] B destination point of the vector.
   * \pre A.dimension() == B.dimension()
   * \pre A.dimension() == ndims
   *****************************************************************************
   */
  Vector( const Point< T,DIM >& A, const Point< T,DIM >& B );

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
  Vector< T,DIM > operator*(T scalar);

  /*!
   *****************************************************************************
   * \brief Computes the vector sum \f$\vec{r}=\vec{u}+\vec{v}\f$
   * \param [in] v the vector on the right hand side.
   * \return r the resulting vector after the sum operation.
   *****************************************************************************
   */
  Vector< T,DIM > operator+( const Vector<T,DIM>& v );

  /*!
   *****************************************************************************
   * \brief Computes the vector subtraction \f$\vec{r}=\vec{u}+\vec{v}\f$
   * \param [in] v the vector on the right hand side.
   * \return r the resulting vector after the subtraction operation.
   *****************************************************************************
   */
  Vector< T,DIM > operator-( const Vector<T,DIM>& v );

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
   * \note The zero vector becomes the unit vector (1,0,0,...) when normalized.
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
  static T dot_product( const Vector< T,DIM >& u,
                        const Vector< T,DIM >& v );

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
template < typename T, int DIM >
Vector< T,DIM >::Vector( const Point< T,DIM >& A,
                           const Point< T,DIM >& B )
{
  SLIC_ASSERT( A.dimension() == B.dimension( ) );

  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] = B[ i ]-A[ i ];
  }

}


//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T,DIM > Vector< T,DIM >::operator*( T scalar )
{
  Vector< T, DIM > result;
  for( int i=0; i < DIM; ++i ) {
    result[ i ] = this->m_components[i] * scalar;
  }
  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T, DIM > Vector< T,DIM >::operator+(const Vector<T,DIM>& v)
{
  Vector< T, DIM > result;
  for ( int i=0; i < DIM; ++i ) {
    result[ i ] = this->m_components[ i ] + v[ i ];
  }

  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T, DIM > Vector< T,DIM >::operator-(const Vector<T,DIM>& v)
{
  Vector< T, DIM > result;
  for ( int i=0; i < DIM; ++i ) {
    result[ i ] = this->m_components[ i ] - v[ i ];
  }

  return( result );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline double Vector< T, DIM >::squared_norm() const
{
   return dot_product( *this, *this);
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline double Vector< T, DIM >::norm() const
{
  return std::sqrt( squared_norm() );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline void Vector< T, DIM >::normalize()
{
  static const double EPS  = 1.0e-100;

  const double length = this->norm();
  if(length <= EPS)
  {
      // normalizing the zero vector gives the unit vector
      this->m_components[ 0 ] = static_cast<T>(1.);
      std::fill( this->m_components+1, this->m_components+DIM, T());
  }
  else
  {
    for ( int i=0; i < DIM; ++i )
      this->m_components[ i ] /= length;
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T Vector< T,DIM >::dot_product( const Vector< T,DIM >& u,
                                       const Vector< T,DIM >& v )
{
  SLIC_ASSERT( u.dimension() == v.dimension() );
  double dotprod = 0.0;

  for ( int i=0; i < DIM; ++i ) {
     dotprod += u[ i ]*v[ i ];
  } // END for all DIM

  return static_cast<T>(dotprod);
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T,3 > Vector< T,DIM >::cross_product(
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
