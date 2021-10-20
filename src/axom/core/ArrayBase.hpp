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
template <typename T, int DIM, typename ArrayType>
class ArrayBase;

/// \name Overloaded ArrayBase Operator(s)
/// @{

/*! 
 * \brief Overloaded output stream operator. Outputs the Array-like to the
 *  given output stream.
 *
 * \param [in,out] os output stream object.
 * \param [in] arr user-supplied Array-like instance.
 * \return os the updated output stream object.
 */
template <typename T, int DIM, typename ArrayType>
std::ostream& operator<<(std::ostream& os,
                         const ArrayBase<T, DIM, ArrayType>& arr);

/*!
 * \brief Equality comparison operator for Array-likes
 *
 * \param [in] lhs left Array-like to compare
 * \param [in] rhs right Array-like to compare
 * \return true if the Arrays have the same allocator ID, are of equal shape,
 * and have the same elements.
 */
template <typename T, int DIM, typename LArrayType, typename RArrayType>
bool operator==(const ArrayBase<T, DIM, LArrayType>& lhs,
                const ArrayBase<T, DIM, RArrayType>& rhs);

/*!
 * \brief Inequality comparison operator for Arrays
 *
 * \param [in] lhs left Array to compare
 * \param [in] rhs right Array to compare
 * \return true if the Arrays do not have the same allocator ID, are not of
 * equal shape, or do not have the same elements.
 */
template <typename T, int DIM, typename LArrayType, typename RArrayType>
bool operator!=(const ArrayBase<T, DIM, LArrayType>& lhs,
                const ArrayBase<T, DIM, RArrayType>& rhs);

/// @}

/*
 * \brief Policy class for implementing Array behavior that differs
 * between the 1D and multidimensional cases
 * 
 * \tparam T The element/value type
 * \tparam DIM The dimension of the Array
 * \tparam ArrayType The type of the underlying array
 * 
 * \pre ArrayType must provide methods with the following signatures:
 * \code{.cpp}
 * IndexType size() const;
 * T* data();
 * const T* data() const;
 * int getAllocatorID() const;
 * \endcode
 */
template <typename T, int DIM, typename ArrayType>
class ArrayBase
{
public:
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
    assert(inBounds(idx));
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
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }

  /*!
   * \brief Accessor, returns a reference to the given value.
   * For multidimensional arrays, indexes into the (flat) raw data.
   *
   * \param [in] idx the position of the value to return.
   *
   * \note equivalent to *(array.data() + idx).
   *
   * \pre 0 <= idx < m_num_elements
   */
  /// @{
  T& operator[](const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  const T& operator[](const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// @}

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
  template <typename OtherArrayType>
  void insert(IndexType pos, const ArrayBase<T, DIM, OtherArrayType>& other)
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
    asDerived().insert(pos,
                       static_cast<const OtherArrayType&>(other).size(),
                       static_cast<const OtherArrayType&>(other).data());
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

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if idx is within bounds */
  inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < asDerived().size();
  }
  /// @}

protected:
  /// \brief The sizes (extents?) in each dimension
  std::array<IndexType, DIM> m_dims;
  /// \brief The strides in each dimension
  std::array<IndexType, DIM> m_strides;
};

/// \brief Array implementation specific to 1D Arrays
template <typename T, typename ArrayType>
class ArrayBase<T, 1, ArrayType>
{
public:
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

  /*!
   * \brief Accessor, returns a reference to the given value.
   * For multidimensional arrays, indexes into the (flat) raw data.
   *
   * \param [in] idx the position of the value to return.
   *
   * \note equivalent to *(array.data() + idx).
   *
   * \pre 0 <= idx < m_num_elements
   */
  /// @{
  T& operator[](const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  const T& operator[](const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// @}

  /// \brief Swaps two ArrayBases
  void swap(ArrayBase&) { }

  /*!
   * \brief Appends an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <typename OtherArrayType>
  void insert(IndexType pos, const ArrayBase<T, 1, OtherArrayType>& other)
  {
    asDerived().insert(pos,
                       static_cast<const OtherArrayType&>(other).size(),
                       static_cast<const OtherArrayType&>(other).data());
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

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if idx is within bounds */
  inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < asDerived().size();
  }
  /// @}
};

//------------------------------------------------------------------------------
//                            ArrayBase IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/// Free functions implementing ArrayBase's operator(s)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM, typename ArrayType>
inline std::ostream& print(std::ostream& os,
                           const ArrayBase<T, DIM, ArrayType>& array)
{
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  // FIXME: Re-add check for umpire::resource::Constant as well, but this will crash
  // if there exists no allocator for Constant memory. Is there a more fine-grained
  // approach we can use to see what allocators are available before trying to get their IDs?
  if(static_cast<const ArrayType&>(array).getAllocatorID() ==
     axom::getUmpireResourceAllocatorID(umpire::resource::Device))
  {
    std::cerr << "Cannot print Array allocated on the GPU" << std::endl;
    utilities::processAbort();
  }
#endif
  const T* data = static_cast<const ArrayType&>(array).data();
  os << "[ ";
  for(IndexType i = 0; i < static_cast<const ArrayType&>(array).size(); i++)
  {
    os << data[i] << " ";
  }
  os << " ]";

  return os;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, typename ArrayType>
std::ostream& operator<<(std::ostream& os, const ArrayBase<T, DIM, ArrayType>& arr)
{
  print(os, arr);
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, typename LArrayType, typename RArrayType>
bool operator==(const ArrayBase<T, DIM, LArrayType>& lhs,
                const ArrayBase<T, DIM, RArrayType>& rhs)
{
  if(static_cast<const LArrayType&>(lhs).getAllocatorID() !=
     static_cast<const RArrayType&>(rhs).getAllocatorID())
  {
    return false;
  }

  if(lhs.shape() != rhs.shape())
  {
    return false;
  }

  for(int i = 0; i < static_cast<const LArrayType&>(lhs).size(); i++)
  {
    if(!(lhs[i] == rhs[i]))
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, typename LArrayType, typename RArrayType>
bool operator!=(const ArrayBase<T, DIM, LArrayType>& lhs,
                const ArrayBase<T, DIM, RArrayType>& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template <typename T, typename ArrayType>
inline void ArrayBase<T, 1, ArrayType>::push_back(const T& value)
{
  emplace_back(value);
}

//------------------------------------------------------------------------------
template <typename T, typename ArrayType>
inline void ArrayBase<T, 1, ArrayType>::push_back(T&& value)
{
  emplace_back(std::move(value));
}

//------------------------------------------------------------------------------
template <typename T, typename ArrayType>
template <typename... Args>
inline void ArrayBase<T, 1, ArrayType>::emplace_back(Args&&... args)
{
  asDerived().emplace(asDerived().size(), args...);
}

namespace detail
{
// Takes the product of a parameter pack of T
template <typename T, int N>
T packProduct(const T (&arr)[N])
{
  return std::accumulate(arr, arr + N, 1, std::multiplies<T> {});
}

template <typename T, int N>
bool allNonNegative(const T (&arr)[N])
{
  for(int i = 0; i < N; i++)
  {
    if(arr[i] < 0)
    {
      return false;
    }
  }
  return true;
}

}  // namespace detail

} /* namespace axom */

#endif /* AXOM_ARRAYBASE_HPP_ */
