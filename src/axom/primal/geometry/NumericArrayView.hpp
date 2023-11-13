// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_NUMERIC_ARRAY_VIEW_HPP_
#define AXOM_PRIMAL_NUMERIC_ARRAY_VIEW_HPP_

#include "axom/primal/geometry/NumericArrayBase.hpp"

#include "axom/core/Macros.hpp"
#include "axom/core/ArrayView.hpp"
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
 * \class NumericArrayView
 *
 * \brief A simple statically sized array of data with component-wise operators.
 *
 * \tparam T the numeric type of the elements in the array, e.g., float, double.
 * \tparam SIZE the size of the array
 */
template <typename T, int SIZE>
class NumericArrayView : public NumericArrayBase<T, SIZE, NumericArrayView<T, SIZE>>
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
  explicit NumericArrayView(T* data, int stride = 1);

  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE
  explicit NumericArrayView(ArrayView<T>& view);

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  AXOM_HOST_DEVICE
  const T* data() const;

  AXOM_HOST_DEVICE
  T* data();

  AXOM_HOST_DEVICE
  T& component(int i) { return m_components[i*m_stride]; }
  
  AXOM_HOST_DEVICE
  const T& component(int i) const { return m_components[i*m_stride]; }

protected:
  T* m_components = nullptr;  /// The encapsulated array
  int m_stride;
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
AXOM_HOST_DEVICE NumericArrayView<T, SIZE>::NumericArrayView(T* data, int stride)
  : m_components {data}, m_stride {stride}
{ 
  // NOTE (KW): This should be a static assert in the class
  SLIC_ASSERT(SIZE >= 1);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE NumericArrayView<T, SIZE>::NumericArrayView(ArrayView<T>& view)
  : m_components {view.data()}, m_stride {view.minStride()}
{
  SLIC_ASSERT(SIZE >= 1);
  SLIC_ASSERT(view.size() >= SIZE);
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline const T* NumericArrayView<T, SIZE>::data() const
{
  return m_components;
}

//------------------------------------------------------------------------------
template <typename T, int SIZE>
AXOM_HOST_DEVICE inline T* NumericArrayView<T, SIZE>::data()
{
  return m_components;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::NumericArrayView using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::NumericArrayView<T, NDIMS>>
  : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_NUMERIC_ARRAY_VIEW_HPP_
