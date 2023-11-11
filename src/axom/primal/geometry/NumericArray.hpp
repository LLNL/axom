// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_NUMERIC_ARRAY_HPP_
#define AXOM_PRIMAL_NUMERIC_ARRAY_HPP_

#include "axom/primal/geometry/NumericArrayBase.hpp"

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
/*!
 * \accelerated
 * \class NumericArray
 *
 * \brief A simple statically sized array of data with component-wise operators.
 *
 * \tparam T the numeric type of the elements in the array, e.g., float, double.
 * \tparam SIZE the size of the array
 */
template <typename T, int SIZE>
class NumericArray : public NumericArrayBase<T, SIZE, NumericArray<T, SIZE>>
{
public:
  /*!
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to. Defaults to zero.
   * \param [in] sz The number of components to set to val.
   * The rest will be set to zero.  Defaults is SIZE.
   * If sz is greater than SIZE, we set all coordinates to val
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  explicit NumericArray(T val = T(), int sz = SIZE);

  /*!
   * \brief Creates a numeric array from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz number of coordinates. Defaults to SIZE.
   * \note If sz is greater than SIZE, we only take the first SIZE values.
   */
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  explicit NumericArray(const T* vals, int sz = SIZE);

  /*!
   * \brief Creates a numeric array from an initializer list
   * \param [in] values an initializer list containing the values of the
   * array. If the size is not the same as the size of this array, this
   * behaves the same way as the constructor which takes a pointer and size.
   */
  NumericArray(std::initializer_list<T> values)
    : NumericArray {values.begin(), static_cast<int>(values.size())}
  { }

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE
  const T* data() const;

  AXOM_HOST_DEVICE
  T* data();

  AXOM_HOST_DEVICE
  T& component(int i) { return m_components[i]; }
  
  AXOM_HOST_DEVICE
  const T& component(int i) const { return m_components[i]; }

protected:
  T m_components[SIZE];  /// The encapsulated array
};

}  // namespace primal
}  // namespace axom

//------------------------------------------------------------------------------
//  NumericArray implementation
//------------------------------------------------------------------------------

namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE>::NumericArray(T val, int sz)
{
  // NOTE (KW): This should be a static assert in the class
  SLIC_ASSERT(SIZE >= 1);

  // Fill first nvals coordinates with val ( 0 <= nvals <= SIZE )
  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);
  for(int i = 0; i < nvals; i++)
  {
    m_components[i] = val;
  }

  // Fill any remaining coordinates with zero
  for(int j = nvals; j < SIZE; j++)
  {
    m_components[j] = T();
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArray<T, SIZE>::NumericArray(const T* vals, int sz)
{
  SLIC_ASSERT(SIZE >= 1);

  const int nvals = axom::utilities::clampVal(sz, 0, SIZE);

  // Copy first nvals coordinates from vals array ( 0 <= nvals <= SIZE )
  for(int i = 0; i < nvals; i++)
  {
    m_components[i] = vals[i];
  }

  // Fill any remaining coordinates with zero
  for(int j = nvals; j < SIZE; j++)
  {
    m_components[j] = T();
  }
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline const T* NumericArray<T, SIZE>::data() const
{
  return m_components;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T* NumericArray<T, SIZE>::data()
{
  return m_components;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_NUMERIC_ARRAY_HPP_
