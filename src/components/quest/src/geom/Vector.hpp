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


// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Vector;

/*!
 * \brief Forward declaration for vector addition
 */
template<typename T, int DIM> Vector<T,DIM> operator+(const Vector<T, DIM> & vec1, const Vector<T, DIM> & vec2);

/*!
 * \brief Forward declaration for vector subtraction
 */
template<typename T, int DIM> Vector<T,DIM> operator-(const Vector<T, DIM> & vec1, const Vector<T, DIM> & vec2);

/*!
 * \brief Forward declaration for vector (unary) negation
 */
template<typename T, int DIM> Vector<T,DIM> operator-(const Vector<T, DIM> & vec1, const Vector<T, DIM> & vec2);

/*!
 * \brief Forward declaration for scalar multiplication of vector; Scalar on rhs.
 */
template<typename T, int DIM> Vector<T,DIM> operator*(const Vector<T, DIM> & vec, const T scalar);

/*!
 * \brief Forward declaration for scalar multiplication of vector; Scalar on lhs.
 */
template<typename T, int DIM> Vector<T,DIM> operator*(const T scalar, const Vector<T, DIM> & vec);

/*!
 * \brief Forward declaration for scalar multiplication of vector; Scalar on rhs.
 */
template<typename T, int DIM> Vector<T,DIM> operator/(const Vector<T, DIM> & vec, const T scalar);

/*!
 * \brief Overloaded output operator for vectors
 */
template<typename T, int DIM> std::ostream& operator<<(std::ostream & os, const Vector<T,DIM> & vec);


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
    * \brief Adds the vector to the Vector instance \f$\vec{u} +=\vec{v}\f$
    * \param [in] v the vector to add.
    * \return A reference to the Vector instance after vector addition.
    *****************************************************************************
    */
   Vector< T,DIM >& operator+=( const Vector<T,DIM>& v );

   /*!
    *****************************************************************************
    * \brief Adds the vector ot the Vector instance \f$\vec{u} -=\vec{v}\f$
    * \param [in] v the vector to subtract.
    * \return A reference to the Vector instance after vector subtraction.
    *****************************************************************************
    */
   Vector< T,DIM >& operator-=( const Vector<T,DIM>& v );

  /*!
   *****************************************************************************
   * \brief Scalar multiplication on the Vector instance.
   * \param [in] scalar the scalar value to multiply with this vector.
   * \return A reference to the vector instance after scalar multiplication.
   *****************************************************************************
   */
  Vector< T,DIM>& operator*=(T scalar);

  /*!
   *****************************************************************************
   * \brief Scalar division on the Vector instance.
   * \param [in] scalar the scalar value to divide with this vector.
   * \pre scalar != 0
   * \return A reference to the vector instance after scalar division.
   *****************************************************************************
   */
  Vector< T,DIM>& operator/=(T scalar);

  /*!
   *****************************************************************************
   * \brief Dot product of the Vector instance with another vector v
   * \param [in] v the other vector in the dot product
   * \return The dot product of the two vectors.
   *****************************************************************************
   */
  T dot(const Vector<T,DIM>& v ) const;


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
   * \brief Component-wise negation of the vector.
   *****************************************************************************
   */
  void negate();

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
   * \brief Simple formatted print of a Vector instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const;



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


  /*!
   *****************************************************************************
   * \brief Creates a unit vector of dimension DIM.
   * \return A unit vector whose first coordinate is 1 and all others are zero.
   *****************************************************************************
   */
  static Vector unit_vector() { return Vector(1,1);}
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
inline Vector< T,DIM >& Vector< T,DIM >::operator*=( T scalar )
{
    for ( int i=0; i < DIM; ++i )
        this->m_components[ i ] *= scalar;

    return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T,DIM >& Vector< T,DIM >::operator/=( T scalar )
{
    SLIC_ASSERT(scalar != 0.);

    return operator*=( 1./scalar );
}


//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T, DIM >& Vector< T,DIM >::operator+=(const Vector<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] +=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Vector< T, DIM >& Vector< T,DIM >::operator-=(const Vector<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] -= v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline double Vector< T, DIM >::squared_norm() const
{
   return this->dot(*this);
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
  static const double EPS  = 1.0e-50;

  const double len_sq = this->squared_norm();
  if(len_sq >= EPS)
  {
    operator/=( std::sqrt(len_sq) );
  }
  else
  {
    *this = Vector::unit_vector();
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline void Vector< T, DIM >::negate()
{
  for ( int i=0; i < DIM; ++i )
    this->m_components[ i ] = -this->m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T Vector< T,DIM >::dot(const Vector<T,DIM>& vec) const
{
    return dot_product(*this,vec);
}


//------------------------------------------------------------------------------
template < typename T, int DIM >
std::ostream& Vector< T, DIM >::print(std::ostream& os) const
{
    os <<"[";
    for(int dim=0; dim < DIM -1; ++ dim)
        os << this->m_components[dim] << ",";
    os << this->m_components[DIM-1] << "]";

    return os;
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
  c[ 1 ] = math::determinant( u[0],u[2], v[0],v[2] ) * ( -1.0 );
  c[ 2 ] = math::determinant( u[0],u[1], v[0],v[1] );
  return( c );
}

///  Free functions involving vectors

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator*(const Vector<T,DIM>& vec, const T scalar)
{
  Vector< T, DIM > result(vec);
  result *=scalar;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator*(const T scalar, const Vector<T,DIM>& vec)
{
  Vector< T, DIM > result(vec);
  result *=scalar;

  return result;
}


//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator+(const Vector<T,DIM>& vec1, const Vector<T,DIM>& vec2)
{
  Vector< T, DIM > result(vec1);
  result +=vec2;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator/(const Vector<T,DIM>& vec, const T scalar)
{
  Vector< T, DIM > result(vec);
  result /=scalar;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator-(const Vector<T,DIM>& vec1, const Vector<T,DIM>& vec2)
{
  Vector< T, DIM > result(vec1);
  result -=vec2;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline Vector<T,DIM> operator-(const Vector<T,DIM>& vec1)
{
  Vector< T, DIM > result(vec1);
  result.negate();

  return result;
}

template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Vector<T,DIM> & vec)
{
    vec.print(os);
    return os;
}



}
#endif /* VECTOR_HXX_ */
