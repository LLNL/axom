// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef NUMERIC_ARRAY_HPP_
#define NUMERIC_ARRAY_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <cstring>    // For memcpy()
#include <algorithm>  // For std:: copy and fill
#include <ostream>    // For print() and operator <<

namespace axom
{
namespace primal
{

// Forward declare the templated classes and operator functions
template < typename T, int SIZE >
class NumericArray;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Checks if two numeric arrays are component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs==rhs, otherwise, false.
 */
template < typename T,int SIZE >
bool operator==( const NumericArray< T,SIZE >& lhs,
                 const NumericArray< T,SIZE >& rhs );

/*!
 * \brief Checks if two numeric arrays are *not* component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs!=rhs, otherwise, false.
 */
template < typename T,int SIZE >
bool operator!=( const NumericArray< T,SIZE >& lhs,
                 const NumericArray< T,SIZE >& rhs);

/*!
 * \brief Performs component-wise addition of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array from the component-wise addition.
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator+( const NumericArray< T,SIZE >& lhs,
                                  const NumericArray< T,SIZE >& rhs  );

/*!
 * \brief Performs component-wise subtraction of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \result C resulting numeric array from component-wise subtraction.
 */
template < typename T,int SIZE >
AXOM_HOST_DEVICE
NumericArray< T,SIZE > operator-( const NumericArray< T,SIZE >& lhs,
                                  const NumericArray< T,SIZE >& rhs  );

/*!
 * \brief Unary negation of a numeric array instance.
 * \param [in] arr numeric array instance on the left-hand side.
 * \result C resulting numeric array from unary negation.
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator-( const NumericArray< T,SIZE >& arr );

/*!
 * \brief Scalar multiplication a numeric array; Scalar on rhs.
 * \param [in] arr numeric array instance.
 * \param [in] scalar user-supplied scalar.
 * \return C resutling numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator*( const NumericArray< T, SIZE > & arr,
                                  double scalar );

/*!
 * \brief Scalar multiplication a numeric array; Scalar on lhs.
 * \param [in] scalar user-supplied scalar.
 * \param [in] arr numeric array instance.
 * \return C resulting numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator*( double scalar,
                                  const NumericArray< T, SIZE > & arr );

/*!
 * \brief Component-wise multiplication of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i * rhs_i, \forall i\f$
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator*( const NumericArray< T, SIZE > & lhs,
                                  const NumericArray< T, SIZE > & rhs  );

/*!
 * \brief Component-wise division of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i / rhs_i, \forall i\f$
 * \pre \f$ rhs_i != 0.0, \forall i \f$
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator/( const NumericArray< T, SIZE > & lhs,
                                  const NumericArray< T, SIZE > & rhs  );

/*!
 * \brief Scalar division of NumericArray; Scalar on rhs
 * \param [in] arr numeric array instance
 * \param [in] scalar user-supplied scalar
 * \return C resulting numeric array, \f$ \ni: C_i = arr_i/scalar, \forall i\f$
 * \pre scalar != 0.0
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > operator/( const NumericArray< T, SIZE >& arr,
                                  double scalar );

/*!
 * \brief Coordinate-wise absolute value on the NumericArray
 * \param [in] arr numeric array instance
 * \pre std::abs is defined for template type T
 * \return A NumericArray whose coordinates are the absolute value of arr
 */
template < typename T,int SIZE >
NumericArray< T,SIZE > abs( const NumericArray< T, SIZE >& arr);

/*!
 * \brief Overloaded output operator for numeric arrays
 * \param [in] os C++ output stream
 * \param [in] arr numeric array instance.
 */
template < typename T,int SIZE >
std::ostream& operator<<( std::ostream & os,
                          const NumericArray< T,SIZE > & arr );

///@}

/**
 * \brief Type trait to avoid outputting chars when a value is expected
 *  This avoids unintentionally outputting system beeps
 */
template < typename T >
struct NonChar
{
  typedef T type;     /** The non-char type to return */
};

template < >
struct NonChar< char >
{
  /** A non-char signed type to which we can cast a char for output */
  typedef int type;
};

template < >
struct NonChar < unsigned char >
{
  /** A non-char unsigned type to which we can cast a char for output */
  typedef unsigned int type;
};

/*!
 * \class NumericArray
 *
 * \brief A simple statically sized array of data with component-wise operators.
 *
 * \tparam T the numeric type of the elements in the array, e.g., float, double.
 * \tparam SIZE the size of the array
 */
template < typename T,int SIZE >
class NumericArray
{
public:
  enum
  {
    NBYTES = SIZE*sizeof(T)
  };

public:

  // -- TODO: Add static_assert that T has numeric type --

  /*!
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to. Defaults to zero.
   * \param [in] sz The number of components to set to val.
   * The rest will be set to zero.  Defaults is SIZE.
   * If sz is greater than SIZE, we set all coordinates to val
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  explicit NumericArray( T val = T(), int sz = SIZE);

  /*!
   * \brief Creates a numeric array from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz number of coordinates. Defaults to SIZE.
   * \note If sz is greater than SIZE, we only take the first SIZE values.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  NumericArray(const T* vals, int sz = SIZE);

  /*!
   * \brief Copy constructor.
   * \param [in] other The numeric array to copy
   */
  AXOM_HOST_DEVICE
  NumericArray( const NumericArray& other ) { *this = other; };

  /*!
   * \brief Destructor.
   */
  AXOM_HOST_DEVICE
  ~NumericArray() { }

  /*!
   * \brief Returns the dimension of this numeric array instance.
   * \return d the dimension (size) of the array
   * \post d >= 1.
   */
  static int size() { return SIZE; };

  /*!
   * \brief Assignment operator.
   * \param [in] rhs a numeric array instance on the right hand side.
   */
  AXOM_HOST_DEVICE
  NumericArray& operator=(const NumericArray& rhs);

  /*!
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return \f$ p_i \f$ the value at the given component index.
   * \pre \f$  0 \le i < SIZE \f$
   */
  AXOM_HOST_DEVICE 
  const T& operator[](int i) const;
  
  AXOM_HOST_DEVICE 
  T& operator[](int i);

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE 
  const T* data() const;
  
  AXOM_HOST_DEVICE 
  T* data();

  /*!
   *
   * \brief Copy the coordinate data to the provided array
   * \param [in] arr The array to which we are copying.
   * \pre The user needs to make sure that the provided array has been allocated
   * and has sufficient space for SIZE coordinates.
   */
  void to_array(T* arr) const;

  /*!
   * \brief Simple formatted print of a numeric array instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /*!
   * \brief Component-wise addition assignment operator.
   * \param [in] arr the array to add.
   * Adds the numeric array arr to this instance (component-wise).
   * \return A reference to the NumericArray instance after addition.
   */
  NumericArray< T,SIZE >& operator+=( const NumericArray< T,SIZE >& arr );

  /*!
   * \brief Component-wise subtraction assignment operator.
   * \param [in] arr the array to subtract.
   * Subtracts the numeric array arr from this instance (component-wise).
   * \return A reference to the NumericArray instance after subtraction.
   */
  AXOM_HOST_DEVICE
  NumericArray< T,SIZE >& operator-=( const NumericArray< T,SIZE >& arr );

  /*!
   * \brief Scalar multiplication on the NumericArray instance.
   * \param [in] scalar the scalar value with which to multiply.
   * Each element of the numeric array is multiplied by scalar
   * \return A reference to the NumericArray instance after scalar
   * multiplication.
   */
  AXOM_HOST_DEVICE
  NumericArray< T,SIZE >& operator*=(double scalar);

  /*!
   * \brief Scalar division on the NumericArray instance.
   * \param [in] scalar the scalar value with which to divide .
   * \pre scalar != 0
   * Each element of the numeric array is divided by scalar
   * \return A reference to the NumericArray instance after scalar division.
   */
  AXOM_HOST_DEVICE
  NumericArray< T,SIZE >& operator/=(double scalar);

  /*!
   * \brief Component-wise multiplication assignment operator.
   * \param [in] arr the array to multiply (component-wise).
   * Multiplies the numeric array arr with this instance (component-wise).
   * \return A reference to the NumericArray instance after cwise
   * multiplication.
   */
  AXOM_HOST_DEVICE
  NumericArray< T,SIZE >& operator*=( const NumericArray< T,SIZE >& arr );

  /*!
   * \brief Component-wise division assignment operator.
   * \param [in] arr the array to divide (component-wise).
   * Divides the numeric array arr with this instance (component-wise).
   * \pre forall i, arr[i] != 0
   * \return A reference to the NumericArray instance after cwise division.
   */
  NumericArray< T,SIZE >& operator/=( const NumericArray< T,SIZE >& arr );

  /*!
   * \brief Ensures that the highest value of the coordinates is at most
   *  upperVal.
   *
   * \param [in] upperVal The highest possible value
   * \post forall i, arr[i] <= upperVal
   * \return A reference to the NumericArray instance after clamping upper
   */
  NumericArray< T,SIZE >& clampUpper( const T& upperVal);

  /*!
   * \brief Ensures that the lowest value of the coordinates is at least
   *  lowerVal.
   *
   * \param [in] lowerVal The lowest possible value
   *
   * \post forall i, arr[i] >= lowerVal
   *
   * \return A reference to the NumericArray instance after clamping lower
   */
  NumericArray< T,SIZE >& clampLower( const T& lowerVal);

  /*!
   * \brief Ensures that each coordinate's value is in range
   *  [lowerVal,upperVal].
   *
   * \param [in] lowerVal The lowest possible value
   * \param [in] upperVal The highest possible value
   *
   * \pre lowerVal <= upperVal
   * \post forall i, lowerVal <= arr[i] <= upperVal
   *
   * \return A reference to the NumericArray instance after clamping
   */
  NumericArray< T,SIZE >& clamp( const T& lowerVal, const T& upperVal);

  /*!
   * \brief Find the max component.
   * \return The value of the largest component.
   */
  T max() const;

  /*!
   * \brief Find the min component.
   * \return The value of the smallest component.
   */
  T min() const;

  /*!
   * \brief Find the index of the max component.
   * \return The index of the largest component ( \f$ 0 \le ret < SIZE \f$)
   */
  int argMax() const;

  /*!
   * \brief Find the index of the min component.
   * \return The index of the smallest component ( \f$ 0 \le ret < SIZE \f$)
   */
  int argMin() const;

private:
  AXOM_HOST_DEVICE
  void verifyIndex(int AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT(idx >= 0 && idx < SIZE);
  }

protected:
  T m_components[ SIZE ];    /*! The encapsulated array */
};

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  NumericArray implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{

//------------------------------------------------------------------------------
template < typename T, int SIZE >
NumericArray< T,SIZE >::NumericArray(T val, int sz)
{
  // NOTE (KW): This should be a static assert in the class
  SLIC_ASSERT( SIZE >= 1 );

  // Fill first nvals coordinates with val ( 0 <= nvals <= SIZE )
  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);
  for (int i = 0; i < nvals; i++)
  {
    m_components[i] = val;    
  }

  // Fill any remaining coordinates with zero
  for (int j = nvals; j < SIZE; j++)
  {
    m_components[j] = T();    
  }
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
NumericArray< T, SIZE >::NumericArray(const T* vals, int sz)
{
  SLIC_ASSERT( SIZE >= 1 );

  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);

  // Copy first nvals coordinates from vals array ( 0 <= nvals <= SIZE )
  for (int i = 0; i < nvals; i++)
  {
    m_components[i] = vals[i];    
  }

  // Fill any remaining coordinates with zero
  for (int j = nvals; j < SIZE; j++) 
  {
    m_components[j] = T();    
  }
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline  AXOM_HOST_DEVICE 
NumericArray< T,SIZE >&
NumericArray< T,SIZE >::operator=( const NumericArray< T,SIZE >& rhs )
{

  if ( this == &rhs )
  {
    return *this;
  }

  // copy all the data
  memcpy( m_components, rhs.m_components, NBYTES);
  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline T& NumericArray< T,SIZE >::operator[](int i)
{
  verifyIndex(i);
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline const T& NumericArray< T,SIZE >::operator[](int i) const
{
  verifyIndex(i);
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline const T* NumericArray< T,SIZE >::data() const
{
  return m_components;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline T* NumericArray< T,SIZE >::data()
{
  return m_components;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
void NumericArray< T,SIZE >::to_array(T* arr) const
{
  SLIC_ASSERT( arr != nullptr);
  memcpy( arr, m_components, NBYTES );
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
std::ostream& NumericArray< T, SIZE >::print(std::ostream& os) const
{
  os <<"[ ";
  for (int dim=0 ; dim < SIZE -1 ; ++dim)
  {
    os << static_cast< typename NonChar< T >::type >( m_components[dim] )
       << " ";
  }

  os << static_cast< typename NonChar< T >::type >(m_components[SIZE-1])
     << "]";

  return os;
}

//------------------------------------------------------------------------------
// Member function arithmetic operators (component-wise)
//------------------------------------------------------------------------------

template < typename T,int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::operator*=( double scalar )
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] = static_cast<T>(m_components[ i ] * scalar);
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::operator/=( double scalar )
{
  SLIC_ASSERT(scalar != 0.);
  return operator*=( 1./scalar );
}

//------------------------------------------------------------------------------
template < typename T, int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::operator*=(const NumericArray< T,SIZE >& v)
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] *=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray<  T,SIZE >&
NumericArray< T,SIZE >::operator/=( const NumericArray< T,SIZE >& v )
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    SLIC_ASSERT( v[ i ] != 0.);
    m_components[ i ] /=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::operator+=(const NumericArray< T,SIZE >& v)
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] +=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int SIZE >
inline NumericArray< T, SIZE >&
NumericArray< T,SIZE >::operator-=(const NumericArray< T,SIZE >& v)
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] -= v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T, SIZE >&
NumericArray< T,SIZE >::clamp( const T& lowerVal, const T& upperVal )
{
  SLIC_ASSERT( lowerVal <= upperVal);

  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] =
      axom::utilities::clampVal(m_components[ i ],lowerVal, upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::clampLower( const T& lowerVal)
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] = std::max( m_components[ i ], lowerVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE >&
NumericArray< T,SIZE >::clampUpper( const T& upperVal)
{
  for ( int i=0 ; i < SIZE ; ++i )
  {
    m_components[ i ] = std::min( m_components[ i ], upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline T NumericArray< T,SIZE >::max() const
{
  T result = this->m_components[0];
  for ( int i=1 ; i < SIZE ; ++i )
  {

    T tmp = m_components[i];

    if ( tmp > result)
    {
      result = tmp;
    }

  }

  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline T NumericArray< T,SIZE >::min() const
{
  T result = this->m_components[0];
  for ( int i=1 ; i < SIZE ; ++i )
  {
    T tmp = this->m_components[i];

    if ( tmp < result)
    {
      result = tmp;
    }

  }

  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline int NumericArray< T,SIZE >::argMax() const
{
  int idx = 0;
  for ( int i=1 ; i < SIZE ; ++i )
  {
    if ( m_components[i] > m_components[idx])
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline int NumericArray< T,SIZE >::argMin() const
{
  int idx = 0;
  for ( int i=1 ; i < SIZE ; ++i )
  {
    if ( m_components[i] < m_components[idx] )
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template < typename T,int SIZE >
bool operator==( const NumericArray< T,SIZE >& lhs,
                 const NumericArray< T,SIZE >& rhs)
{
  for ( int dim=0 ; dim < SIZE ; ++dim )
  {
    if ( lhs[dim] != rhs[dim] )
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
bool operator!=( const NumericArray< T,SIZE >& lhs,
                 const NumericArray< T,SIZE >& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
std::ostream& operator<<(std::ostream & os, const NumericArray< T,SIZE > & arr)
{
  arr.print(os);
  return os;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator*( const NumericArray< T,SIZE >& arr,
                                         double scalar)
{
  NumericArray< T,SIZE > result(arr);
  result *=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator*( double scalar,
                                         const NumericArray< T,SIZE >& arr)
{
  NumericArray< T, SIZE > result(arr);
  result *=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator+( const NumericArray< T,SIZE >& lhs,
                                         const NumericArray< T,SIZE >& rhs)
{
  NumericArray< T, SIZE > result(lhs);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator*( const NumericArray< T,SIZE >& lhs,
                                         const NumericArray< T,SIZE >& rhs)
{
  NumericArray< T,SIZE > result(lhs);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator/( const NumericArray< T,SIZE >& lhs,
                                         const NumericArray< T,SIZE >& rhs)
{
  NumericArray< T,SIZE > result(lhs);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator/( const NumericArray< T,SIZE >& arr,
                                         double scalar)
{
  NumericArray< T, SIZE > result(arr);
  result /=scalar;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator-( const NumericArray< T,SIZE >& lhs,
                                         const NumericArray< T,SIZE >& rhs)
{
  NumericArray< T,SIZE > result(lhs);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > operator-(const NumericArray< T,SIZE >& arr)
{
  NumericArray< T, SIZE > result;
  result -= arr;
  return result;
}

//------------------------------------------------------------------------------
template < typename T,int SIZE >
inline NumericArray< T,SIZE > abs(const NumericArray< T,SIZE >& arr)
{
  NumericArray< T, SIZE > result(arr);

  for (int i=0 ; i<SIZE ; ++i)
  {
    result[i] = axom::utilities::abs(result[i]);
  }

  return result;
}

} /* namespace primal*/

} /* namespace axom */

#endif /* NUMERIC_ARRAY_HXX_ */
