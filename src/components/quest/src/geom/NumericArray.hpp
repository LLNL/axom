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
 *******************************************************************************
 */

#ifndef NUMERIC_ARRAY_HXX_
#define NUMERIC_ARRAY_HXX_

#include "slic/slic.hpp"

// C/C++ includes
#include <cstring>      // For memcpy()
#include <algorithm>    // For std:: copy and fill
#include <ostream>      // For print() and operator <<


namespace {
    /*!
     *****************************************************************************
     * \brief Utility function that clamps an input val to a given range.
     * \param [in] val  The value to clamp
     * \param [in] lower The lower range
     * \param [in] upper The upper range
     * \return The clamped value.
     * \post lower <= returned value <= upper.
     *****************************************************************************
     */
    template<typename T>
    T clamp(T val, T lower, T upper)
    {
        SLIC_ASSERT( lower <= upper);
        return std::min( std::max( val, lower), upper);
    }
}


namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class NumericArray;

/*!
* \brief Forward declaration of NumericArray's equality operator
* Two numeric arrays are considered equal when they are component-wise equal
*/
template<typename T, int DIM> bool operator==(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);

/*!
* \brief Forward declaration of NumericArray's inequality operator
*/
template<typename T, int DIM> bool operator!=(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);



/*!
 * \brief Forward declaration for NumericArray's addition
 */
template<typename T, int DIM> NumericArray<T,DIM> operator+(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);

/*!
 * \brief Forward declaration for NumericArray's subtraction
 */
template<typename T, int DIM> NumericArray<T,DIM> operator-(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);

/*!
 * \brief Forward declaration for NumericArray's (unary) negation
 */
template<typename T, int DIM> NumericArray<T,DIM> operator-(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);

/*!
 * \brief Forward declaration for scalar multiplication of NumericArray's; Scalar on rhs.
 */
template<typename T, int DIM> NumericArray<T,DIM> operator*(const NumericArray<T, DIM> & arr, const T scalar);

/*!
 * \brief Forward declaration for scalar multiplication of NumericArray; Scalar on lhs.
 */
template<typename T, int DIM> NumericArray<T,DIM> operator*(const T scalar, const NumericArray<T, DIM> & arr);

/*!
 * \brief Forward declaration for component-wise multiplication of NumericArrays
 */
template<typename T, int DIM> NumericArray<T,DIM> operator*(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);

/*!
 * \brief Forward declaration for component-wise division of NumericArrays
 */
template<typename T, int DIM> NumericArray<T,DIM> operator/(const NumericArray<T, DIM> & lhs, const NumericArray<T, DIM> & rhs);


/*!
 * \brief Forward declaration for scalar division of NumericArray; Scalar on rhs.
 */
template<typename T, int DIM> NumericArray<T,DIM> operator/(const NumericArray<T, DIM> & arr, const T scalar);


/*!
 * \brief Overloaded output operator for numeric arrays
 */
template<typename T, int DIM> std::ostream& operator<<(std::ostream & os, const NumericArray<T,DIM> & arr);


/*!
 *******************************************************************************
 * \class NumericArray
 *
 * \brief A simple statically sized array of data with component-wise operators
 *******************************************************************************
 */
template < typename T, int DIM >
class NumericArray
{
public:
    enum { NDIMS = DIM
         , NBYTES = DIM * sizeof(T)
         };

public:

    // -- TODO: Add static_assert that T has numeric type --


  /*!
   *****************************************************************************
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to.  Defaults to zero ( T())
   * \param [in] sz The number of components to set to val.
   * The rest will be set to zero.  Defaults is DIM.
   * If sz is greater than DIM, we set all coordinates to val
   *****************************************************************************
   */
  explicit NumericArray(T val = T(), int sz = DIM);

  /*!
   *****************************************************************************
   * \brief Creates a numeric array from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz The number of coordinates to take from the array.  Defaults
   * to DIM.
   * It sz is greater than DIM, we only take the first DIM values.
   *****************************************************************************
   */
  NumericArray(T* vals, int sz = DIM);

  /*!
   *****************************************************************************
   * \brief Copy constructor.
   * \param [in] rhs
   *****************************************************************************
   */
  NumericArray( const NumericArray& rhs ) { *this = rhs; };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
   ~NumericArray() {}

  /*!
   *****************************************************************************
   * \brief Returns the dimension of this numeric array instance.
   * \return d the dimension (size) of the array
   * \post d >= 1.
   *****************************************************************************
   */
  int dimension() const { return DIM; };

  /*!
   *****************************************************************************
   * \brief Assignment operator.
   * \param [in] rhs a numeric array instance on the right hand side.
   *****************************************************************************
   */
  NumericArray& operator=(const NumericArray& rhs);

  /*!
   *****************************************************************************
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   *****************************************************************************
   */
  const T& operator[](int i) const;
  T& operator[](int i);


  /*!
   *****************************************************************************
   * \brief Returns a pointer to the underlying data.
   *****************************************************************************
   */
  const T* data() const;
  T* data();

  /*!
   *****************************************************************************
   * \brief Copy the coordinate data to the provided array
   * \param arr The array to which we are copying.
   * \pre The user needs to make sure that the provided array has been allocated
   * and has sufficient space for DIM coordinates.
   *****************************************************************************
   */
  void to_array(T* arr) const;

  /*!
   *****************************************************************************
   * \brief Simple formatted print of a numeric array instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const;



  /*!
   *****************************************************************************
   * \brief Component-wise addition assignment operator.
   * \param [in] arr the array to add.
   * Adds the numeric array @arr to this instance (component-wise).
   * \return A reference to the NumericArray instance after addition.
   *****************************************************************************
   */
  NumericArray< T,DIM >& operator+=( const NumericArray<T,DIM>& arr );

  /*!
   *****************************************************************************
   * \brief Component-wise subtraction assignment operator.
   * \param [in] arr the array to subtract.
   * Subtracts the numeric array @arr from this instance (component-wise).
   * \return A reference to the NumericArray instance after subtraction.
   *****************************************************************************
   */
  NumericArray< T,DIM >& operator-=( const NumericArray<T,DIM>& arr );

 /*!
  *****************************************************************************
  * \brief Scalar multiplication on the NumericArray instance.
  * \param [in] scalar the scalar value with which to multiply.
  * Each element of the numeric array is multiplied by @scalar
  * \return A reference to the NumericArray instance after scalar multiplication.
  *****************************************************************************
  */
 NumericArray< T,DIM>& operator*=(T scalar);

 /*!
  *****************************************************************************
  * \brief Scalar division on the NumericArray instance.
  * \param [in] scalar the scalar value with which to divide .
  * \pre scalar != 0
  * Each element of the numeric array is divided by @scalar
  * \return A reference to the NumericArray instance after scalar division.
  *****************************************************************************
  */
 NumericArray< T,DIM>& operator/=(T scalar);

 /*!
  *****************************************************************************
  * \brief Component-wise multiplication assignment operator.
  * \param [in] arr the array to multiply (component-wise).
  * Multiplies the numeric array @arr with this instance (component-wise).
  * \return A reference to the NumericArray instance after cwise multiplication.
  *****************************************************************************
  */
 NumericArray< T,DIM >& operator*=( const NumericArray<T,DIM>& arr );

 /*!
  *****************************************************************************
  * \brief Component-wise division assignment operator.
  * \param [in] arr the array to divide (component-wise).
  * Divides the numeric array @arr with this instance (component-wise).
  * \pre forall i, arr[i] != 0
  * \return A reference to the NumericArray instance after cwise division.
  *****************************************************************************
  */
 NumericArray< T,DIM >& operator/=( const NumericArray<T,DIM>& arr );

 /*!
  *****************************************************************************
  * \brief Find the max component.
  * \return The value of the largest component.
  *****************************************************************************
  */
 T max() const;

 /*!
  *****************************************************************************
  * \brief Find the min component.
  * \return The value of the smallest component.
  *****************************************************************************
  */
 T min() const;

 /*!
  *****************************************************************************
  * \brief Find the index of the max component.
  * \return The index of the largest component ( 0 <= ret < DIM)
  *****************************************************************************
  */
 int argMax() const;

 /*!
  *****************************************************************************
  * \brief Find the index of the min component.
  * \return The index of the smallest component ( 0 <= ret < DIM)
  *****************************************************************************
  */
 int argMin() const;

private:
  void verifyIndex(int idx) const { SLIC_ASSERT(idx >= 0 && idx < DIM); }

protected:
  T m_components[ DIM ];
};

} /* namespace quest */

//------------------------------------------------------------------------------
//  Point implementation
//------------------------------------------------------------------------------
namespace quest {

//------------------------------------------------------------------------------
template < typename T, int DIM >
NumericArray< T, DIM >::NumericArray(T val, int sz)
{
  // NOTE (KW): This should be a static assert in the class
  SLIC_ASSERT( DIM >= 1 );

  // Fill first nvals coordinates with val ( 0 <= nvals <= DIM )
  const int nvals = clamp(sz, 0, DIM);
  std::fill( m_components, m_components+nvals, val );

  // Fill any remaining coordinates with zero
  if(nvals < DIM)
  {
      std::fill( m_components+nvals, m_components+DIM, T() );
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
NumericArray< T, DIM >::NumericArray(T* vals, int sz)
{
  SLIC_ASSERT( DIM >= 1 );

  const int nvals = clamp(sz, 0, DIM);

  // Copy first nvals coordinates from vals array ( 0 <= nvals <= DIM )
  std::copy( vals, vals+nvals, m_components);

  // Fill any remaining coordinates with zero
  if(nvals < DIM)
  {
      std::fill( m_components+nvals, m_components+DIM, T());
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray<T,DIM>& NumericArray< T,DIM >::operator=(const NumericArray<T,DIM>& rhs )
{

  if( this == &rhs ) {
    return *this;
  }

  // copy all the data
  memcpy( m_components, rhs.m_components, NBYTES);
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T& NumericArray< T, DIM >::operator[](int i)
{
    verifyIndex(i);
    return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T& NumericArray< T, DIM >::operator[](int i) const
{
  verifyIndex(i);
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T* NumericArray< T, DIM >::data() const
{
  return m_components;
}
//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T* NumericArray< T, DIM >::data()
{
    return m_components;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
void NumericArray< T, DIM >::to_array(T* arr) const
{
    SLIC_ASSERT( arr != ATK_NULLPTR);
    memcpy( arr, m_components, NBYTES );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
std::ostream& NumericArray< T, DIM >::print(std::ostream& os) const
{
    os <<"[";
    for(int dim=0; dim < DIM -1; ++ dim)
        os << m_components[dim] << ",";
    os << m_components[DIM-1] << "]";

    return os;
}


//------------------------------------------------------------------------------
// Member function arithmetic operators (component-wise)
//------------------------------------------------------------------------------

template < typename T, int DIM >
inline NumericArray< T,DIM >& NumericArray< T,DIM >::operator*=( T scalar )
{
    for ( int i=0; i < DIM; ++i )
        this->m_components[ i ] *= scalar;

    return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray< T,DIM >& NumericArray< T,DIM >::operator/=( T scalar )
{
    SLIC_ASSERT(scalar != 0.);

    return operator*=( 1./scalar );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray< T, DIM >& NumericArray< T,DIM >::operator*=(const NumericArray<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] *=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray< T, DIM >& NumericArray< T,DIM >::operator/=(const NumericArray<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      SLIC_ASSERT( v[i] != 0.);
      this->m_components[ i ] /=  v[ i ];
  }

  return *this;
}


//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray< T, DIM >& NumericArray< T,DIM >::operator+=(const NumericArray<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] +=  v[ i ];
  }

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline NumericArray< T, DIM >& NumericArray< T,DIM >::operator-=(const NumericArray<T,DIM>& v)
{
  for ( int i=0; i < DIM; ++i ) {
      this->m_components[ i ] -= v[ i ];
  }

  return *this;
}


//------------------------------------------------------------------------------
template<typename T, int DIM>
inline T NumericArray< T,DIM >::max() const
{
  T result = this->m_components[0];
  for ( int i=1; i < DIM; ++i )
  {
    T tmp = this->m_components[i];
    if( (tmp) > result)
        result = tmp;
  }

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline T NumericArray< T,DIM >::min() const
{
  T result = this->m_components[0];
  for ( int i=1; i < DIM; ++i )
  {
    T tmp = this->m_components[i];
    if( (tmp) < result)
        result = tmp;
  }

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline int NumericArray< T,DIM >::argMax() const
{
  int idx = 0;
  for ( int i=1; i < DIM; ++i )
  {
    if( this->m_components[i] > this->m_components[idx])
        idx = i;
  }

  return idx;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline int NumericArray< T,DIM >::argMin() const
{
  int idx = 0;
  for ( int i=1; i < DIM; ++i )
  {
    if( this->m_components[i] < this->m_components[idx])
        idx = i;
  }

  return idx;
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
bool operator==(const NumericArray<T, DIM>& lhs, const NumericArray<T, DIM>& rhs)
{
    for(int dim=0;dim<DIM;++dim)
    {
        if( lhs[dim] != rhs[dim])
            return false;
    }
    return true;
}

//------------------------------------------------------------------------------

template<typename T, int DIM>
bool operator!=(const NumericArray<T, DIM>& lhs, const NumericArray<T, DIM>& rhs)
{
    return !(lhs == rhs);
}


template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const NumericArray<T,DIM> & arr)
{
    arr.print(os);
    return os;
}


///  Free functions involving vectors

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator*(const NumericArray<T,DIM>& arr, const T scalar)
{
  NumericArray< T, DIM > result(arr);
  result *=scalar;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator*(const T scalar, const NumericArray<T,DIM>& arr)
{
  NumericArray< T, DIM > result(arr);
  result *=scalar;

  return result;
}


//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator+(const NumericArray<T,DIM>& lhs, const NumericArray<T,DIM>& rhs)
{
  NumericArray< T, DIM > result(lhs);
  result +=rhs;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator*(const NumericArray<T,DIM>& lhs, const NumericArray<T,DIM>& rhs)
{
  NumericArray< T, DIM > result(lhs);
  result *=rhs;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator/(const NumericArray<T,DIM>& lhs, const NumericArray<T,DIM>& rhs)
{
  NumericArray< T, DIM > result(lhs);
  result /=rhs;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator/(const NumericArray<T,DIM>& arr, const T scalar)
{
  NumericArray< T, DIM > result(arr);
  result /=scalar;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator-(const NumericArray<T,DIM>& lhs, const NumericArray<T,DIM>& rhs)
{
  NumericArray< T, DIM > result(lhs);
  result -=rhs;

  return result;
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
inline NumericArray<T,DIM> operator-(const NumericArray<T,DIM>& arr)
{
  NumericArray< T, DIM > result;
  result -= arr;

  return result;
}


} /* namespace quest*/

#endif /* NUMERIC_ARRAY_HXX_ */
