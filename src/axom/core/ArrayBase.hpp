// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYBASE_HPP_
#define AXOM_ARRAYBASE_HPP_

#include "axom/config.hpp"                    // for compile-time defines
#include "axom/core/Macros.hpp"               // for axom macros
#include "axom/core/utilities/Utilities.hpp"  // for processAbort()
#include "axom/core/Types.hpp"                // for IndexType definition

// C/C++ includes
#include <array>     // for std::array
#include <iostream>  // for std::cerr and std::ostream
#include <numeric>   // for std::inner_product

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM>
class Array;

/*!
 * \brief Policy class for implementing Array behavior that differs
 * between the 1D and multidimensional cases
 * 
 * \tparam T The element/value type
 * \tparam DIM The dimension of the Array
 */
template <typename T, int DIM>
class ArrayBase
{
public:
  using ArrayType = Array<T, DIM>;

  /*!
   * \brief Parameterized constructor that sets up the default strides
   *
   * \param [in] args the parameter pack of sizes in each dimension.
   */
  template <typename... Args>
  ArrayBase(Args... args) : m_dims {static_cast<IndexType>(args)...}
  {
    updateStrides();
  }

  /*!
   * \brief Dimension-aware accessor, returns a reference to the given value.
   *
   * \param [in] args the parameter pack of indices in each dimension.
   *
   * \note equivalent to *(array.data() + idx).
   *
   * \pre sizeof...(Args) == DIM
   * \pre 0 <= args[i] < m_dims[i] for i in [0, DIM)
   */
  template <typename... Args,
            typename SFINAE = typename std::enable_if<sizeof...(Args) == DIM>::type>
  T& operator()(Args... args)
  {
    IndexType indices[] = {static_cast<IndexType>(args)...};
    IndexType idx =
      std::inner_product(indices, indices + DIM, m_strides.begin(), 0);
    // assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  template <typename... Args,
            typename SFINAE = typename std::enable_if<sizeof...(Args) == DIM>::type>
  const T& operator()(Args... args) const
  {
    IndexType indices[] = {static_cast<IndexType>(args)...};
    IndexType idx =
      std::inner_product(indices, indices + DIM, m_strides.begin(), 0);
    // assert(inBounds(idx));
    return asDerived().data()[idx];
  }

  /// \brief Swaps two ArrayBases
  friend void swap(ArrayBase& lhs, ArrayBase& rhs)
  {
    std::swap(lhs.m_dims, rhs.m_dims);
    std::swap(lhs.m_strides, rhs.m_strides);
  }

  /// \brief Returns the dimensions of the Array
  const std::array<IndexType, DIM>& shape() const { return m_dims; }

  /// \brief Returns the strides of the Array
  const std::array<IndexType, DIM>& strides() const { return m_strides; }

  /*!
   * \brief Appends an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   * 
   * \pre The shapes of the calling Array and @a other are the same
   * (excluding the leading dimension), i.e., shape()[1:] == other.shape()[1:]
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void insert(IndexType pos, const ArrayType& other)
  {
#ifdef AXOM_DEBUG
    if(!std::equal(m_dims.begin() + 1, m_dims.end(), other.shape().begin() + 1))
    {
      std::cerr << "Cannot append a multidimensional array of incorrect shape.";
      utilities::processAbort();
    }
#endif

    // First update the dimensions - we're adding only to the leading dimension
    m_dims[0] += other.shape()[0];
    // Then add the raw data to the buffer
    asDerived().insert(pos, other.size(), other.data());
    updateStrides();
  }

protected:
  /*!
   * \brief Returns the minimum "chunk size" that should be allocated
   * For example, 2 would be the chunk size of a 2D array whose second dimension is of size 2.
   * This is used when resizing/reallocating; it wouldn't make sense to have a
   * capacity of 3 in the array described above.
   */
  IndexType blockSize() const { return m_strides[0]; }

  /*!
   * \brief Updates the internal striding information to a row-major format
   * Intended to be called after @p m_dims is updated.
   * In the future, this class will support different striding schemes (e.g., column-major)
   * and/or user-provided striding
   */
  void updateStrides()
  {
    // Row-major
    m_strides[DIM - 1] = 1;
    for(int i = static_cast<int>(DIM) - 2; i >= 0; i--)
    {
      m_strides[i] = m_strides[i + 1] * m_dims[i + 1];
    }
  }

private:
  /// \brief Returns a reference to the Derived CRTP object - see https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/
  ArrayType& asDerived() { return static_cast<ArrayType&>(*this); }
  /// \overload
  const ArrayType& asDerived() const
  {
    return static_cast<const ArrayType&>(*this);
  }

protected:
  /// \brief The sizes (extents?) in each dimension
  std::array<IndexType, DIM> m_dims;
  /// \brief The strides in each dimension
  std::array<IndexType, DIM> m_strides;
};

/// \brief Array implementation specific to 1D Arrays
template <typename T>
class ArrayBase<T, 1>
{
public:
  using ArrayType = Array<T, 1>;
  ArrayBase(IndexType = 0) { }

  /*!
   * \brief Push a value to the back of the array.
   *
   * \param [in] value the value to be added to the back.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void push_back(const T& value);

  /*!
   * \brief Push a value to the back of the array.
   *
   * \param [in] value the value to move to the back.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void push_back(T&& value);

  /*!
   * \brief Inserts new element at the end of the Array.
   *
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   */
  template <typename... Args>
  void emplace_back(Args&&... args);

  /// \brief Returns the dimensions of the Array
  // FIXME: std::array is used for consistency with multidim case, should we just return the scalar?
  // Double curly braces needed for C++11 prior to resolution of CWG issue 1720
  std::array<IndexType, 1> shape() const { return {{asDerived().size()}}; }

  /// \brief Swaps two ArrayBases
  void swap(ArrayBase&) { }

  /*!
   * \brief Appends an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void insert(IndexType pos, const ArrayType& other)
  {
    asDerived().insert(pos, other.size(), other.data());
  }

protected:
  /*!
   * \brief Returns the minimum "chunk size" that should be allocated
   */
  IndexType blockSize() const { return 1; }

private:
  /// \brief Returns a reference to the Derived CRTP object - see https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/
  ArrayType& asDerived() { return static_cast<ArrayType&>(*this); }
  /// \overload
  const ArrayType& asDerived() const
  {
    return static_cast<const ArrayType&>(*this);
  }
};

} /* namespace axom */

#endif /* AXOM_ARRAYBASE_HPP_ */
