// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYVIEW_HPP_
#define AXOM_ARRAYVIEW_HPP_

#include "axom/core/memory_management.hpp"  // for memory allocation functions
#include "axom/core/ArrayBase.hpp"
#include "axom/core/ArrayIteratorBase.hpp"

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM>
class ArrayView;

/// \name ArrayView to wrap a pointer and provide indexing semantics
/// @{

/*!
 * \class ArrayView
 *
 * \brief Provides a view over a generic array container.
 * 
 * The ArrayView expresses a non-owning relationship over a pointer
 *
 * \tparam T the type of the values to hold.
 * \tparam DIM The dimension of the array.
 *
 */
template <typename T, int DIM = 1>
class ArrayView : public ArrayBase<T, DIM, ArrayView<T, DIM>>
{
public:
  using value_type = T;
  static constexpr int dimension = DIM;
  using ArrayViewIterator = ArrayIteratorBase<ArrayView<T, DIM>>;

  /// \brief Default constructor
  ArrayView() = default;

  /*!
   * \brief Generic constructor for an ArrayView of arbitrary dimension with external data
   *
   * \param [in] data the external data this ArrayView will wrap.
   * \param [in] args The parameter pack containing the "shape" of the ArrayView
   *
   * \pre sizeof...(Args) == DIM
   *
   * \post size() == num_elements
   */
  template <typename... Args>
  ArrayView(T* data, Args... args);

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  inline IndexType size() const { return m_num_elements; }

  /*!
   * \brief Returns an ArrayViewIterator to the first element of the Array
   */
  ArrayViewIterator begin()
  {
    assert(m_data != nullptr);
    return ArrayViewIterator(0, this);
  }

  /*!
   * \brief Returns an ArrayViewIterator to the element following the last
   *  element of the Array.
   */
  ArrayViewIterator end()
  {
    assert(m_data != nullptr);
    return ArrayViewIterator(size(), this);
  }

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  inline T* data() { return m_data; }
  inline const T* data() const { return m_data; }

  /// @}

  /*!
   * \brief Get the ID for the umpire allocator
   * 
   * FIXME: This is just a stand-in impl, extend this class to support wrapping of GPU pointers
   */
  int getAllocatorID() const { return axom::getDefaultAllocatorID(); }

private:
  T* m_data = nullptr;
  /// \brief The full number of elements in the array
  ///  i.e., 3 for a 1D Array of size 3, 9 for a 3x3 2D array, etc
  IndexType m_num_elements = 0;
};

/// \brief Helper alias for multi-component arrays
template <typename T>
using MCArrayView = ArrayView<T, 2>;

//------------------------------------------------------------------------------
//                            ArrayView IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM>
template <typename... Args>
ArrayView<T, DIM>::ArrayView(T* data, Args... args)
  : ArrayBase<T, DIM, ArrayView<T, DIM>>(args...)
  , m_data(data)
{
  static_assert(sizeof...(Args) == DIM,
                "Array size must match number of dimensions");
  // Intel hits internal compiler error when casting as part of function call
  IndexType tmp_args[] = {args...};
  m_num_elements = detail::packProduct(tmp_args);
}

} /* namespace axom */

#endif /* AXOM_ARRAYVIEW_HPP_ */
