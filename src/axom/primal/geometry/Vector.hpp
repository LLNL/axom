// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

// axom_utils includes
#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/numerics/matvecops.hpp"

// Primal includes
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"

// C/C++ includes
#include <cmath> // for sqrt()

namespace axom
{
namespace primal
{

// Forward declare the templated classes and operator functions
template < typename T,int NDIMS >
class Vector;

/// \name Forward Declared Overloaded Operators
/// @{

/*!
 * \brief Adds vectors A, B and stores the result into a new vector C.
 * \param [in] A vector on the left-hand side.
 * \param [in] B vector on the right-hand side.
 * \return C resulting vector, \f$ C_i = A_i + B_i \forall i \f$
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator+( const Vector< T,NDIMS >& A,
                             const Vector< T,NDIMS >& B  );

/*!
 * \brief Subtracts vectors A, B and stores the result into a new vector C
 * \param [in] A vector on the left-hand side.
 * \param [in] B vector on the right-hand side.
 * \return C resulting vector, \f$ C_i = A_i - B_i \forall i \f$
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator-( const Vector< T,NDIMS >& A,
                             const Vector< T,NDIMS >& B  );

/*!
 * \brief Unary negation of a vector instance.
 * \param [in] vec1 vector instance to negate.
 * \return C resulting vector from unary negation.
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator-( const Vector< T, NDIMS > & vec1 );

/*!
 * \brief Scalar multiplication of vector; Scalar on rhs.
 * \param [in] vec vector instance.
 * \param [in] scalar user-supplied scalar.
 * \return C resulting vector, \f$ \ni: C_i = scalar*vec_i, \forall i\f$
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator*(const Vector< T, NDIMS >& vec, const T scalar);

/*!
 * \brief Scalar multiplication of vector; Scalar on lhs.
 * \param [in] scalar user-supplied scalar.
 * \param [in] vec vector instance.
 * \return C resulting vector, \f$ \ni: C_i = scalar*vec_i, \forall i\f$
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator*(const T scalar, const Vector< T, NDIMS > & vec);

/*!
 * \brief Scalar division of vector; Scalar on rhs.
 * \param [in] vec vector instance
 * \param [in]n scalar user-supplied scalar.
 * \return C resulting vector, \f$ \ni: C_i = vec_i / scalar, \forall i\f$
 * \pre scalar != 0.0
 */
template < typename T,int NDIMS >
Vector< T,NDIMS > operator/(const Vector< T, NDIMS > & vec, const T scalar);

/*!
 * \brief Overloaded output operator for vectors
 * \param [in] os C++ output stream
 * \param [in] vec vector instance
 */
template < typename T,int NDIMS >
std::ostream& operator<<(std::ostream & os, const Vector< T,NDIMS > & vec);

/// @}

/*!
 * \class Vector
 *
 * \brief Represents a vector, \f$ v \in \mathcal{R}^d \f$. It provides access
 *  methods for setting and querying the vector components as well as vector
 *  math operators, e.g., adding, subtracting, dot_product and cross_product.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \see NumericArray
 * \see Point
 */
template < typename T,int NDIMS >
class Vector
{
public:
  typedef Point< T,NDIMS > PointType;

public:

  /*!
   * \brief Fill vector with single value. Acts as default constructor
   * Sets first sz components of the vector to val (default 0).
   * \param [in] val The value to set the coordinates to.  Defaults to zero
   * \param [in] sz The number of coordinates to set to val.
   * The rest will be set to zero.  Defaults is NDIMS.
   */
  AXOM_HOST_DEVICE
  explicit Vector( T val = T(), int sz = NDIMS) : m_components(val, sz) { }

  /*!
   * \brief Constructor from a numeric array
   * \param [in] arr The numeric array to copy from
   */
  AXOM_HOST_DEVICE
  Vector( const NumericArray< T,NDIMS >& arr) : m_components(arr) { }

  /*!
   * \brief Creates a vector from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz The number of coordinates to take from the array.  Defaults
   * to NDIMS.
   * If sz is greater than NDIMS, we only take the first NDIMS values.
   */
  AXOM_HOST_DEVICE
  Vector(const T* vals, int sz = NDIMS) : m_components(vals, sz) {}

  /*!
   * \brief Constructor to create vector from a Point
   * \param [in] pt The point containing the vector's coordinates.
   * \note Equivalent to Vector( Point::zero(), pt)
   */
  AXOM_HOST_DEVICE
  Vector(const Point< T,NDIMS > & pt) : m_components( pt.array() ) {}

  /*!
   * \brief Copy constructor
   * \param [in] other The vector to copy
   */
  AXOM_HOST_DEVICE
  Vector(const Vector< T,NDIMS > & other) : m_components( other.array() ) {}

  /*!
   * \brief Constructs a vector from point A to point B.
   * \param [in] A origin point of the vector.
   * \param [in] B destination point of the vector.
   * \pre A.dimension() == B.dimension()
   * \pre A.dimension() == ndims
   */
  AXOM_HOST_DEVICE
  Vector( const Point< T,NDIMS >& A, const Point< T,NDIMS >& B ) :
    m_components(B.array() - A.array()) {}

  /*!
   * \brief Destructor.
   */
  AXOM_HOST_DEVICE
  ~Vector() {}

  /*!
   * \brief Returns the dimension of this vector instance.
   * \return d the dimension (size) of the vector
   * \post d >= 1.
   */
  AXOM_HOST_DEVICE
  int dimension() const { return NDIMS; };

  /*!
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   */
  AXOM_HOST_DEVICE
  const T& operator[](int i) const { return m_components[i]; }
  
  AXOM_HOST_DEVICE
  T& operator[](int i)             { return m_components[i]; }

  /*!
   * \brief Returns a reference to the underlying NumericArray.
   */
  AXOM_HOST_DEVICE 
  const NumericArray< T,NDIMS > & array() const { return m_components; }
  
  AXOM_HOST_DEVICE 
  NumericArray< T,NDIMS >& array()              { return m_components; }

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE
  const T* data() const { return m_components.data(); }
  
  AXOM_HOST_DEVICE
  T* data()             { return m_components.data(); }

  /*!
   * \brief Equality comparison operator for vectors.
   */
  friend bool operator==(const Vector& lhs, const Vector& rhs)
  { return lhs.m_components == rhs.m_components; }

  /*!
   * \brief Inequality operator for points
   */
  friend bool operator!=(const Vector& lhs, const Vector& rhs)
  { return !(lhs == rhs); }

  /*!
   * \brief Adds the vector to the Vector instance \f$\vec{u} +=\vec{v}\f$
   * \param [in] v the vector to add.
   * \return A reference to the Vector instance after vector addition.
   */
  Vector< T,NDIMS >& operator+=( const Vector< T,NDIMS >& v );

  /*!
   * \brief Adds the vector ot the Vector instance \f$\vec{u} -=\vec{v}\f$
   * \param [in] v the vector to subtract.
   * \return A reference to the Vector instance after vector subtraction.
   */
  Vector< T,NDIMS >& operator-=( const Vector< T,NDIMS >& v );

  /*!
   * \brief Scalar multiplication on the Vector instance.
   * \param [in] scalar the scalar value to multiply with this vector.
   * \return A reference to the vector instance after scalar multiplication.
   */
  Vector< T,NDIMS >& operator*=(T scalar);

  /*!
   * \brief Scalar division on the Vector instance.
   * \param [in] scalar the scalar value to divide with this vector.
   * \pre scalar != 0
   * \return A reference to the vector instance after scalar division.
   */
  AXOM_HOST_DEVICE
  Vector< T,NDIMS >& operator/=(T scalar);

  /*!
   * \brief Dot product of the Vector instance with another vector v
   * \param [in] v the other vector in the dot product
   * \return The dot product of the two vectors.
   */
  AXOM_HOST_DEVICE
  T dot(const Vector< T,NDIMS >& v ) const;

  /*!
   * \brief Computes the squared \f$ l^2\f$ norm of this vector instance.
   * \return n the squared norm.
   * \see Vector::norm()
   */
  AXOM_HOST_DEVICE
  double squared_norm() const;

  /*!
   * \brief Computes the \f$ l^2 \f$ norm of this vector instance.
   * \return n the norm of the vector, a.k.a., magnitude or length.
   */
  AXOM_HOST_DEVICE
  double norm() const;

  /*!
   * \brief Component-wise negation of the vector.
   */
  void negate();

  /*!
   * \brief Creates a new unit vector in the direction of the vector instance.
   * \note The unit vector of the zero vector is (1,0,0,...) when normalized.
   * \post this->norm() == 1.0f
   */
  AXOM_HOST_DEVICE
  Vector unitVector() const;

  /*!
   * \brief Simple formatted print of a Vector instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /*!
   * \brief Computes the dot product of two vectors u, v.
   * \param [in] u the vector on the right-hand side.
   * \param [in] v the vector on the left-hand side.
   * \return dotprod the computed dot product.
   * \pre u.dimension() == v.dimension().
   */
  AXOM_HOST_DEVICE
  static T dot_product( const Vector< T,NDIMS >& u,
                        const Vector< T,NDIMS >& v );

  /*!
   * \brief Computes the cross product of vector u and v, treating them as 3D.
   * \param [in] u the vector on the right-hand side.
   * \param [in] v the vector on the left-hand side.
   * \return C the resulting vector from A x B.
   */
  AXOM_HOST_DEVICE
  static Vector< T,3 > cross_product( const Vector< T,2 >& u,
                                      const Vector< T,2 >& v );

  /*!
   * \brief Computes the 3-D cross product of vector u and v.
   * \param [in] u the vector on the right-hand side.
   * \param [in] v the vector on the left-hand side.
   * \return C the resulting vector from A x B.
   */
  AXOM_HOST_DEVICE
  static Vector< T,3 > cross_product( const Vector< T,3 >& u,
                                      const Vector< T,3 >& v );

  /*!
   * \brief Utility function to constructs a Vector with the given coordinates.
   * \param [in] x the x--coordinate of the vector.
   * \param [in] y the y--coordinate of the vector.
   * \param [in] z the z--coordinate of the vector. Default is 0.0.
   * \return v a Vector instance with the given coordinates.
   */
  static Vector make_vector( const T& x, const T& y, const T& z=0.0 );

private:
  NumericArray< T,NDIMS > m_components;
};

/// \name Pre-defined Vector types
/// @{

typedef Vector< double, 2 > Vector2D;
typedef Vector< double, 3 > Vector3D;

/// @}

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  Vector implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T,NDIMS >& Vector< T,NDIMS >::operator*=( T scalar )
{
  m_components *= scalar;
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T,NDIMS >& Vector< T,NDIMS >::operator/=( T scalar )
{
  m_components /= scalar;
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T, NDIMS >&
Vector< T,NDIMS >::operator+=(const Vector< T,NDIMS >& v)
{
  m_components += v.array();
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T, NDIMS >&
Vector< T,NDIMS >::operator-=(const Vector< T,NDIMS >& v)
{
  m_components -= v.array();
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline double Vector< T, NDIMS >::squared_norm() const
{
  return dot(*this);
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline double Vector< T, NDIMS >::norm() const
{
  return std::sqrt( squared_norm() );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T, NDIMS > Vector< T, NDIMS >::unitVector() const
{
  constexpr double EPS = 1.0e-50;

  Vector v( *this);

  const double len_sq = squared_norm();
  if ( len_sq >= EPS )
  {
    v /= ( std::sqrt(len_sq) );
  }
  else
  {
    // Create a vector whose first coordinate is 1 and all others are 0
    v = Vector( static_cast< T >( 1. ), 1 );
  }

  return v;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline void Vector< T, NDIMS >::negate()
{
  for ( int i=0 ; i < NDIMS ; ++i )
  {
    m_components[ i ] = -m_components[ i ];
  }
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline T Vector< T,NDIMS >::dot(const Vector< T,NDIMS >& vec) const
{
  return dot_product(*this,vec);
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
std::ostream& Vector< T, NDIMS >::print(std::ostream& os) const
{
  os <<"<";
  for ( int dim=0 ; dim < NDIMS -1 ; ++dim )
  {
    os << static_cast< typename NonChar< T >::type >( m_components[dim] )
       << ",";
  }
  os << static_cast< typename NonChar< T >::type >( m_components[NDIMS-1] )
     << ">";

  return os;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline T Vector< T,NDIMS >::dot_product( const Vector< T,NDIMS >& u,
                                         const Vector< T,NDIMS >& v )
{
  SLIC_ASSERT( u.dimension() == v.dimension() );
  return static_cast< T >( numerics::dot_product( u.data(), v.data(), NDIMS ) );
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,3 > Vector< T,NDIMS >::cross_product(
  const Vector< T,2 >& u, const Vector< T,2 >& v )
{
  Vector< T,3 > c;
  c[ 0 ] = 0;
  c[ 1 ] = 0;
  c[ 2 ] = numerics::determinant( u[0],u[1], v[0],v[1] );
  return( c );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Vector< T,3 > Vector< T,NDIMS >::cross_product(
  const Vector< T,3 >& u, const Vector< T,3 >& v )
{
  Vector< T,3 > c;
  numerics::cross_product( u.data(), v.data(), c.data() );
  return( c );
}

///  Free functions involving vectors

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator*(const Vector< T,NDIMS >& vec, const T scalar)
{
  Vector< T, NDIMS > result(vec);
  result *=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator*(const T scalar, const Vector< T,NDIMS >& vec)
{
  Vector< T, NDIMS > result(vec);
  result *=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator+( const Vector< T,NDIMS >& vec1,
                                    const Vector< T,NDIMS >& vec2)
{
  Vector< T, NDIMS > result(vec1);
  result +=vec2;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator/(const Vector< T,NDIMS >& vec, const T scalar)
{
  Vector< T,NDIMS > result(vec);
  result /=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator-( const Vector< T,NDIMS >& vec1,
                                    const Vector< T,NDIMS >& vec2)
{
  Vector< T, NDIMS > result(vec1);
  result -=vec2;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > operator-(const Vector< T,NDIMS >& vec1)
{
  Vector< T, NDIMS > result(vec1);
  result.negate();
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
std::ostream& operator<<( std::ostream & os,
                          const Vector< T,NDIMS > & vec)
{
  vec.print(os);
  return os;
}

//------------------------------------------------------------------------------
template < typename T,int NDIMS >
inline Vector< T,NDIMS > Vector< T,NDIMS >::make_vector( const T& x,
                                                         const T& y,
                                                         const T& z )
{
  T tmp_array[3] = { x, y, z};
  return Vector(tmp_array, NDIMS);
}

} /* namespace primal */

} /* namespace axom */

#endif /* VECTOR_HXX_ */
