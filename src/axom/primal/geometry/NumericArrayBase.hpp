// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_NUMERIC_ARRAY_BASE_HPP_
#define AXOM_PRIMAL_NUMERIC_ARRAY_BASE_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <algorithm>
#include <ostream>
#include <initializer_list>
#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int SIZE, typename ArrayType>
class NumericArrayBase;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Checks if two numeric arrays are component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs==rhs, otherwise, false.
 */
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE bool operator==(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                 const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Checks if two numeric arrays are *not* component-wise equal.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return status true if lhs!=rhs, otherwise, false.
 */
template <typename T, int SIZE, typename ArrayType>
bool operator!=(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Performs component-wise addition of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array from the component-wise addition.
 */
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE NumericArrayBase<T, SIZE, ArrayType> operator+(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                                                const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Performs component-wise subtraction of two numeric arrays.
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \result C resulting numeric array from component-wise subtraction.
 */
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE NumericArrayBase<T, SIZE, ArrayType> operator-(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                                                const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Unary negation of a numeric array instance.
 * \param [in] arr numeric array instance on the left-hand side.
 * \result C resulting numeric array from unary negation.
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> operator-(const NumericArrayBase<T, SIZE, ArrayType>& arr);

/*!
 * \brief Scalar multiplication a numeric array; Scalar on rhs.
 * \param [in] arr numeric array instance.
 * \param [in] scalar user-supplied scalar.
 * \return C resutling numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> operator*(const NumericArrayBase<T, SIZE, ArrayType>& arr, double scalar);

/*!
 * \brief Scalar multiplication a numeric array; Scalar on lhs.
 * \param [in] scalar user-supplied scalar.
 * \param [in] arr numeric array instance.
 * \return C resulting numeric array, \f$ \ni: C_i = scalar*arr_i, \forall i\f$
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> operator*(double scalar, const NumericArrayBase<T, SIZE, ArrayType>& arr);

/*!
 * \brief Component-wise multiplication of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i * rhs_i, \forall i\f$
 */
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE NumericArrayBase<T, SIZE, ArrayType> operator*(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                                                const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Component-wise division of NumericArrays
 * \param [in] lhs numeric array instance on the left-hand side.
 * \param [in] rhs numeric array instance on the right-hand side.
 * \return C resulting numeric array, \f$ \ni: C_i = lhs_i / rhs_i, \forall i\f$
 * \pre \f$ rhs_i != 0.0, \forall i \f$
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> operator/(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                               const NumericArrayBase<T, SIZE, ArrayType>& rhs);

/*!
 * \brief Scalar division of NumericArray; Scalar on rhs
 * \param [in] arr numeric array instance
 * \param [in] scalar user-supplied scalar
 * \return C resulting numeric array, \f$ \ni: C_i = arr_i/scalar, \forall i\f$
 * \pre scalar != 0.0
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> operator/(const NumericArrayBase<T, SIZE, ArrayType>& arr, double scalar);

/*!
 * \brief Coordinate-wise absolute value on the NumericArray
 * \param [in] arr numeric array instance
 * \pre std::abs is defined for template type T
 * \return A NumericArray whose coordinates are the absolute value of arr
 */
template <typename T, int SIZE, typename ArrayType>
NumericArrayBase<T, SIZE, ArrayType> abs(const NumericArrayBase<T, SIZE, ArrayType>& arr);

/*!
 * \brief Overloaded output operator for numeric arrays
 * \param [in] os C++ output stream
 * \param [in] arr numeric array instance.
 */
template <typename T, int SIZE, typename ArrayType>
std::ostream& operator<<(std::ostream& os, const NumericArrayBase<T, SIZE, ArrayType>& arr);

///@}

/**
 * \brief Type trait to avoid outputting chars when a value is expected
 *  This avoids unintentionally outputting system beeps
 */
template <typename T>
struct NonChar
{
  typedef T type; /** The non-char type to return */
};

template <>
struct NonChar<char>
{
  /** A non-char signed type to which we can cast a char for output */
  typedef int type;
};

template <>
struct NonChar<unsigned char>
{
  /** A non-char unsigned type to which we can cast a char for output */
  typedef unsigned int type;
};

/*!
 * \accelerated
 * \class NumericArrayBase
 *
 * \brief A simple statically sized array of data with component-wise operators.
 *
 * \tparam T the numeric type of the elements in the array, e.g., float, double.
 * \tparam SIZE the size of the array
 */
template <typename T, int SIZE, typename ArrayType>
class NumericArrayBase  // NOLINT
{
public:
  /*!
   * \brief Returns the dimension of this numeric array instance.
   * \return d the dimension (size) of the array
   * \post d >= 1.
   */
  static int size() { return SIZE; };

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
   * \brief Copy the coordinate data to the provided array
   * \param [in] arr The array to which we are copying.
   * \pre The user needs to make sure that the provided array has been allocated
   * and has sufficient space for SIZE coordinates.
   */
  AXOM_HOST_DEVICE
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
  AXOM_HOST_DEVICE
  NumericArrayBase& operator+=(const NumericArrayBase& arr);

  /*!
   * \brief Component-wise subtraction assignment operator.
   * \param [in] arr the array to subtract.
   * Subtracts the numeric array arr from this instance (component-wise).
   * \return A reference to the NumericArray instance after subtraction.
   */
  AXOM_HOST_DEVICE
  NumericArrayBase& operator-=(const NumericArrayBase& arr);

  /*!
   * \brief Scalar multiplication on the NumericArray instance.
   * \param [in] scalar the scalar value with which to multiply.
   * Each element of the numeric array is multiplied by scalar
   * \return A reference to the NumericArray instance after scalar
   * multiplication.
   */
  AXOM_HOST_DEVICE
  NumericArrayBase& operator*=(double scalar);

  /*!
   * \brief Scalar division on the NumericArray instance.
   * \param [in] scalar the scalar value with which to divide .
   * \pre scalar != 0
   * Each element of the numeric array is divided by scalar
   * \return A reference to the NumericArray instance after scalar division.
   */
  AXOM_HOST_DEVICE
  NumericArrayBase& operator/=(double scalar);

  /*!
   * \brief Component-wise multiplication assignment operator.
   * \param [in] arr the array to multiply (component-wise).
   * Multiplies the numeric array arr with this instance (component-wise).
   * \return A reference to the NumericArray instance after cwise
   * multiplication.
   */
  AXOM_HOST_DEVICE
  NumericArrayBase& operator*=(const NumericArrayBase& arr);

  /*!
   * \brief Component-wise division assignment operator.
   * \param [in] arr the array to divide (component-wise).
   * Divides the numeric array arr with this instance (component-wise).
   * \pre forall i, arr[i] != 0
   * \return A reference to the NumericArray instance after cwise division.
   */
  NumericArrayBase& operator/=(const NumericArrayBase& arr);

  /*!
   * \brief Ensures that the highest value of the coordinates is at most
   *  upperVal.
   *
   * \param [in] upperVal The highest possible value
   * \post forall i, arr[i] <= upperVal
   * \return A reference to the NumericArray instance after clamping upper
   */
  NumericArrayBase& clampUpper(const T& upperVal);

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
  NumericArrayBase& clampLower(const T& lowerVal);

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
  NumericArrayBase& clamp(const T& lowerVal, const T& upperVal);

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

  /// \brief Computes the sum of the components
  T sum() const;

protected:
  AXOM_HOST_DEVICE
  void verifyIndex(int AXOM_DEBUG_PARAM(idx)) const
  {
    SLIC_ASSERT(idx >= 0 && idx < SIZE);
  }

private:
  /// \brief Returns a reference to the Derived CRTP object - see https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/
  AXOM_HOST_DEVICE ArrayType& asDerived()
  {
    return static_cast<ArrayType&>(*this);
  }
  /// \overload
  AXOM_HOST_DEVICE const ArrayType& asDerived() const
  {
    return static_cast<const ArrayType&>(*this);
  }
};

}  // namespace primal
}  // namespace axom

//------------------------------------------------------------------------------
//  NumericArrayBase implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline T& NumericArrayBase<T, SIZE, ArrayType>::operator[](int i)
{
  verifyIndex(i);
  return asDerived().component(i);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline const T& NumericArrayBase<T, SIZE, ArrayType>::operator[](int i) const
{
  verifyIndex(i);
  return asDerived().component(i);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE void NumericArrayBase<T, SIZE, ArrayType>::to_array(T* arr) const
{
  SLIC_ASSERT(arr != nullptr);
  for(int dim = 0; dim < SIZE; ++dim)
  {
    arr[dim] = asDerived().component(dim);
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
std::ostream& NumericArrayBase<T, SIZE, ArrayType>::print(std::ostream& os) const
{
  os << "[ ";
  for(int dim = 0; dim < SIZE - 1; ++dim)
  {
    os << static_cast<typename NonChar<T>::type>(asDerived().component(dim)) << " ";
  }

  os << static_cast<typename NonChar<T>::type>(asDerived().component(SIZE - 1)) << "]";

  return os;
}

//------------------------------------------------------------------------------
// Member function arithmetic operators (component-wise)
//------------------------------------------------------------------------------

template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator*=(
  double scalar)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) = static_cast<T>(asDerived().component(i) * scalar);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator/=(
  double scalar)
{
  SLIC_ASSERT(scalar != 0.);
  return operator*=(1. / scalar);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator*=(
  const NumericArrayBase<T, SIZE, ArrayType>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) *= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator/=(
  const NumericArrayBase<T, SIZE, ArrayType>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    SLIC_ASSERT(v[i] != 0.);
    asDerived().component(i) /= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator+=(
  const NumericArrayBase<T, SIZE, ArrayType>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) += v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::operator-=(
  const NumericArrayBase<T, SIZE, ArrayType>& v)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) -= v[i];
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::clamp(const T& lowerVal,
                                                                                         const T& upperVal)
{
  SLIC_ASSERT(lowerVal <= upperVal);

  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) =
      axom::utilities::clampVal(asDerived().component(i), lowerVal, upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::clampLower(const T& lowerVal)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) = std::max(asDerived().component(i), lowerVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType>& NumericArrayBase<T, SIZE, ArrayType>::clampUpper(const T& upperVal)
{
  for(int i = 0; i < SIZE; ++i)
  {
    asDerived().component(i) = std::min(asDerived().component(i), upperVal);
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline T NumericArrayBase<T, SIZE, ArrayType>::max() const
{
  T result = this->asDerived().component(0);
  for(int i = 1; i < SIZE; ++i)
  {
    T tmp = asDerived().component(i);

    if(tmp > result)
    {
      result = tmp;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline T NumericArrayBase<T, SIZE, ArrayType>::min() const
{
  T result = this->asDerived().component(0);
  for(int i = 1; i < SIZE; ++i)
  {
    T tmp = this->asDerived().component(i);

    if(tmp < result)
    {
      result = tmp;
    }
  }

  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline int NumericArrayBase<T, SIZE, ArrayType>::argMax() const
{
  int idx = 0;
  for(int i = 1; i < SIZE; ++i)
  {
    if(asDerived().component(i) > asDerived().component(idx))
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline int NumericArrayBase<T, SIZE, ArrayType>::argMin() const
{
  int idx = 0;
  for(int i = 1; i < SIZE; ++i)
  {
    if(asDerived().component(i) < asDerived().component(idx))
    {
      idx = i;
    }
  }

  return idx;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline T NumericArrayBase<T, SIZE, ArrayType>::sum() const
{
  T result {};
  for(int i = 0; i < SIZE; ++i)
  {
    result += this->asDerived().component(i);
  }

  return result;
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE bool operator==(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                 const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  for(int dim = 0; dim < SIZE; ++dim)
  {
    if(lhs[dim] != rhs[dim])
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
bool operator!=(const NumericArrayBase<T, SIZE, ArrayType>& lhs, const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
std::ostream& operator<<(std::ostream& os, const NumericArrayBase<T, SIZE, ArrayType>& arr)
{
  arr.print(os);
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> operator*(const NumericArrayBase<T, SIZE, ArrayType>& arr,
                                                      double scalar)
{
  NumericArrayBase<T, SIZE, ArrayType> result(arr);
  result *= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> operator*(double scalar,
                                                      const NumericArrayBase<T, SIZE, ArrayType>& arr)
{
  NumericArrayBase<T, SIZE, ArrayType> result(arr);
  result *= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType> operator+(
  const NumericArrayBase<T, SIZE, ArrayType>& lhs,
  const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  NumericArrayBase<T, SIZE, ArrayType> result(lhs);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType> operator*(
  const NumericArrayBase<T, SIZE, ArrayType>& lhs,
  const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  NumericArrayBase<T, SIZE, ArrayType> result(lhs);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> operator/(const NumericArrayBase<T, SIZE, ArrayType>& lhs,
                                                      const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  NumericArrayBase<T, SIZE, ArrayType> result(lhs);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> operator/(const NumericArrayBase<T, SIZE, ArrayType>& arr,
                                                      double scalar)
{
  NumericArrayBase<T, SIZE, ArrayType> result(arr);
  result /= scalar;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
AXOM_HOST_DEVICE inline NumericArrayBase<T, SIZE, ArrayType> operator-(
  const NumericArrayBase<T, SIZE, ArrayType>& lhs,
  const NumericArrayBase<T, SIZE, ArrayType>& rhs)
{
  NumericArrayBase<T, SIZE, ArrayType> result(lhs);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> operator-(const NumericArrayBase<T, SIZE, ArrayType>& arr)
{
  NumericArrayBase<T, SIZE, ArrayType> result;
  result -= arr;
  return result;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE, typename ArrayType>
inline NumericArrayBase<T, SIZE, ArrayType> abs(const NumericArrayBase<T, SIZE, ArrayType>& arr)
{
  NumericArrayBase<T, SIZE, ArrayType> result(arr);

  for(int i = 0; i < SIZE; ++i)
  {
    result[i] = axom::utilities::abs(result[i]);
  }

  return result;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_NUMERIC_ARRAY_BASE_HPP_
