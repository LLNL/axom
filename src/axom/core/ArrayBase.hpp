// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYBASE_HPP_
#define AXOM_ARRAYBASE_HPP_

#include "axom/config.hpp"                    // for compile-time defines
#include "axom/core/Macros.hpp"               // for axom macros
#include "axom/core/memory_management.hpp"    // for memory allocation functions
#include "axom/core/utilities/Utilities.hpp"  // for processAbort()
#include "axom/core/Types.hpp"                // for IndexType definition
#include "axom/core/StackArray.hpp"
#include "axom/core/numerics/matvecops.hpp"  // for dot_product

// C/C++ includes
#include <iostream>  // for std::cerr and std::ostream
#include <numeric>   // for std::accumulate
#ifdef __ibmxl__
  #include <algorithm>  // for std::fill_n due to XL bug
#endif

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM, typename ArrayType>
class ArrayBase;

namespace detail
{
template <typename ArrayType>
struct ArrayTraits;

template <typename T, int DIM>
class ArrayViewProxy;
}  // namespace detail

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
template <typename T1, typename T2, int DIM, typename LArrayType, typename RArrayType>
bool operator==(const ArrayBase<T1, DIM, LArrayType>& lhs,
                const ArrayBase<T2, DIM, RArrayType>& rhs);

/*!
 * \brief Inequality comparison operator for Arrays
 *
 * \param [in] lhs left Array to compare
 * \param [in] rhs right Array to compare
 * \return true if the Arrays do not have the same allocator ID, are not of
 * equal shape, or do not have the same elements.
 */
template <typename T1, typename T2, int DIM, typename LArrayType, typename RArrayType>
bool operator!=(const ArrayBase<T1, DIM, LArrayType>& lhs,
                const ArrayBase<T2, DIM, RArrayType>& rhs);

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
 *
 * \pre A specialization of ArrayTraits for all ArrayTypes must also be provided
 * with the boolean value IsView, which affects the const-ness of returned
 * references.
 */
template <typename T, int DIM, typename ArrayType>
class ArrayBase
{
private:
  constexpr static bool is_array_view = detail::ArrayTraits<ArrayType>::is_view;

public:
  /* If ArrayType is an ArrayView, we use shallow-const semantics, akin to
   * std::span; a const ArrayView will still allow for mutating the underlying
   * pointed-to data.
   *
   * If ArrayType is an Array, we use deep-const semantics, akin to std::vector;
   * a const Array will prevent modifications of the underlying Array data.
   */
  using RealConstT = typename std::conditional<is_array_view, T, const T>::type;

  template <int SliceDim>
  using SliceType = typename std::conditional<
    SliceDim == DIM,
    T&,
    const typename detail::ArrayTraits<ArrayType>::template Slice<SliceDim>>::type;

  template <int SliceDim>
  using ConstSliceType = typename std::conditional<
    SliceDim == DIM,
    RealConstT&,
    const typename detail::ArrayTraits<ArrayType>::template Slice<SliceDim>>::type;

  AXOM_HOST_DEVICE ArrayBase() : m_dims {} { updateStrides(); }

  /*!
   * \brief Parameterized constructor that sets up the default strides
   *
   * \param [in] args the parameter pack of sizes in each dimension.
   */
  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, DIM>& args)
    : m_dims {args}
  {
    updateStrides();
  }

  /*!
   * \brief Copy constructor for arrays of different type
   * Because the element type (T) and dimension (DIM) are still locked down,
   * this function is nominally used for copying ArrayBase metadata from
   * Array <-> ArrayView and/or Array-like objects whose data are in different
   * memory spaces
   */
  template <typename OtherArrayType>
  ArrayBase(
    const ArrayBase<typename std::remove_const<T>::type, DIM, OtherArrayType>& other)
    : m_dims(other.shape())
    , m_strides(other.strides())
  { }

  /// \overload
  template <typename OtherArrayType>
  ArrayBase(
    const ArrayBase<const typename std::remove_const<T>::type, DIM, OtherArrayType>& other)
    : m_dims(other.shape())
    , m_strides(other.strides())
  { }

  /*!
   * \brief Dimension-aware accessor; with N=DIM indices, returns a reference
   *  to the given value at that index. Otherwise, returns a sub-array pointed
   *  to by the given sub-index.
   *
   * \param [in] args the parameter pack of indices in each dimension.
   *
   * \note equivalent to *(array.data() + idx).
   *
   * \pre sizeof...(Args) <= DIM
   * \pre 0 <= args[i] < m_dims[i] for i in [0, sizeof...(Args))
   */
  template <typename... Args>
  AXOM_HOST_DEVICE SliceType<sizeof...(Args)> operator()(Args... args)
  {
    static_assert(sizeof...(Args) <= DIM,
                  "Index dimensions different from array dimensions");
    constexpr int UDim = sizeof...(Args);
    const StackArray<IndexType, UDim> indices {static_cast<IndexType>(args)...};
    return (*this)[indices];
  }
  /// \overload
  template <typename... Args>
  AXOM_HOST_DEVICE ConstSliceType<sizeof...(Args)> operator()(Args... args) const
  {
    static_assert(sizeof...(Args) <= DIM,
                  "Index dimensions different from array dimensions");
    constexpr int UDim = sizeof...(Args);
    const StackArray<IndexType, UDim> indices {static_cast<IndexType>(args)...};
    return (*this)[indices];
  }

  /*!
   * \brief Scalar accessor; returns a sub-array referenced by the given sub-
   *  index, beginning at array(idx, 0...)
   *
   * \param [in] idx the index of the first dimension.
   *
   * \pre 0 <= idx < m_dims[0]
   */
  AXOM_HOST_DEVICE SliceType<1> operator[](const IndexType idx)
  {
    const StackArray<IndexType, 1> slice {idx};
    return (*this)[slice];
  }

  /// \overload
  AXOM_HOST_DEVICE ConstSliceType<1> operator[](const IndexType idx) const
  {
    const StackArray<IndexType, 1> slice {idx};
    return (*this)[slice];
  }

  /*!
   * \brief Dimension-aware accessor; with UDim=DIM indices, returns a reference
   *  to the given value at that index. Otherwise, returns a sub-array pointed
   *  to by the given sub-index.
   *
   * \param [in] args a stack array of indices in each dimension.
   *
   * \note equivalent to *(array.data() + idx).
   *
   * \pre UDim <= DIM
   * \pre 0 <= args[i] < m_dims[i] for i in [0, UDim)
   */
  template <int UDim>
  AXOM_HOST_DEVICE SliceType<UDim> operator[](const StackArray<IndexType, UDim>& idx)
  {
    static_assert(UDim <= DIM,
                  "Index dimensions cannot be larger than array dimensions");
    return sliceImpl(idx);
  }

  /// \overload
  template <int UDim>
  AXOM_HOST_DEVICE ConstSliceType<UDim> operator[](
    const StackArray<IndexType, UDim>& idx) const
  {
    static_assert(UDim <= DIM,
                  "Index dimensions cannot be larger than array dimensions");
    return sliceImpl(idx);
  }

  /// @{

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
  AXOM_HOST_DEVICE T& flatIdx(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIdx(const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// @}

  /// \brief Swaps two ArrayBases
  void swap(ArrayBase& other)
  {
    std::swap(m_dims, other.m_dims);
    std::swap(m_strides, other.m_strides);
  }

  /// \brief Returns the dimensions of the Array
  AXOM_HOST_DEVICE const StackArray<IndexType, DIM>& shape() const
  {
    return m_dims;
  }

  /// \brief Returns the strides of the Array
  AXOM_HOST_DEVICE const StackArray<IndexType, DIM>& strides() const
  {
    return m_strides;
  }

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
  AXOM_HOST_DEVICE void updateStrides()
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
  AXOM_HOST_DEVICE ArrayType& asDerived()
  {
    return static_cast<ArrayType&>(*this);
  }
  /// \overload
  AXOM_HOST_DEVICE const ArrayType& asDerived() const
  {
    return static_cast<const ArrayType&>(*this);
  }

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if idx is within bounds */
  AXOM_HOST_DEVICE inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < asDerived().size();
  }
  /// @}

  /// \name Internal subarray slicing methods
  /// @{

  /*! \brief Returns a subarray given UDim indices */
  template <int UDim>
  AXOM_HOST_DEVICE SliceType<UDim> sliceImpl(const StackArray<IndexType, UDim>& idx)
  {
    const IndexType baseIdx =
      numerics::dot_product((const IndexType*)idx, m_strides.begin(), UDim);
    assert(inBounds(baseIdx));
    StackArray<IndexType, DIM - UDim> new_inds;
    for(int i = 0; i < DIM - UDim; i++)
    {
      new_inds[i] = m_dims[UDim + i];
    }
    detail::ArrayViewProxy<T, DIM - UDim> viewProxy(asDerived().data() + baseIdx,
                                                    asDerived().getAllocatorID(),
                                                    new_inds);
    return SliceType<UDim>(viewProxy);
  }

  /// \overload
  template <int UDim>
  AXOM_HOST_DEVICE ConstSliceType<UDim> sliceImpl(
    const StackArray<IndexType, UDim>& idx) const
  {
    const IndexType baseIdx =
      numerics::dot_product((const IndexType*)idx, m_strides.begin(), UDim);
    assert(inBounds(baseIdx));
    StackArray<IndexType, DIM - UDim> new_inds;
    for(int i = 0; i < DIM - UDim; i++)
    {
      new_inds[i] = m_dims[UDim + i];
    }
    detail::ArrayViewProxy<T, DIM - UDim> viewProxy(asDerived().data() + baseIdx,
                                                    asDerived().getAllocatorID(),
                                                    new_inds);
    return ConstSliceType<UDim>(viewProxy);
  }

  /*! \brief Returns a scalar reference given a full set of indices */
  AXOM_HOST_DEVICE SliceType<DIM> sliceImpl(const StackArray<IndexType, DIM>& idx)
  {
    const IndexType baseIdx =
      numerics::dot_product((const IndexType*)idx, m_strides.begin(), DIM);
    assert(inBounds(baseIdx));
    return asDerived().data()[baseIdx];
  }

  /// \overload
  AXOM_HOST_DEVICE ConstSliceType<DIM> sliceImpl(
    const StackArray<IndexType, DIM>& idx) const
  {
    const IndexType baseIdx =
      numerics::dot_product((const IndexType*)idx, m_strides.begin(), DIM);
    assert(inBounds(baseIdx));
    return asDerived().data()[baseIdx];
  }
  /// @}

protected:
  /// \brief The sizes (extents?) in each dimension
  StackArray<IndexType, DIM> m_dims;
  /// \brief The strides in each dimension
  StackArray<IndexType, DIM> m_strides;
};

/// \brief Array implementation specific to 1D Arrays
template <typename T, typename ArrayType>
class ArrayBase<T, 1, ArrayType>
{
private:
  constexpr static bool is_array_view = detail::ArrayTraits<ArrayType>::is_view;

public:
  /* If ArrayType is an ArrayView, we use shallow-const semantics, akin to
   * std::span; a const ArrayView will still allow for mutating the underlying
   * pointed-to data.
   *
   * If ArrayType is an Array, we use deep-const semantics, akin to std::vector;
   * a const Array will prevent modifications of the underlying Array data.
   */
  using RealConstT = typename std::conditional<is_array_view, T, const T>::type;

  AXOM_HOST_DEVICE ArrayBase(IndexType = 0) { }

  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, 1>&) { }

  // Empy implementation because no member data
  template <typename OtherArrayType>
  ArrayBase(const ArrayBase<typename std::remove_const<T>::type, 1, OtherArrayType>&)
  { }

  // Empy implementation because no member data
  template <typename OtherArrayType>
  ArrayBase(
    const ArrayBase<const typename std::remove_const<T>::type, 1, OtherArrayType>&)
  { }

  /// \brief Returns the dimensions of the Array
  // Double curly braces needed for C++11 prior to resolution of CWG issue 1720
  AXOM_HOST_DEVICE StackArray<IndexType, 1> shape() const
  {
    return {{asDerived().size()}};
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
  AXOM_HOST_DEVICE T& operator[](const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& operator[](const IndexType idx) const
  {
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
  AXOM_HOST_DEVICE T& flatIdx(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIdx(const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// @}

  /// \brief Swaps two ArrayBases
  /// No member data, so this is a no-op
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
  AXOM_HOST_DEVICE ArrayType& asDerived()
  {
    return static_cast<ArrayType&>(*this);
  }
  /// \overload
  AXOM_HOST_DEVICE const ArrayType& asDerived() const
  {
    return static_cast<const ArrayType&>(*this);
  }

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if idx is within bounds */
  AXOM_HOST_DEVICE inline bool inBounds(IndexType idx) const
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
#if defined(AXOM_USE_UMPIRE) && defined(UMPIRE_ENABLE_DEVICE)
  const int alloc_id = static_cast<const ArrayType&>(array).getAllocatorID();
  if(alloc_id == axom::getUmpireResourceAllocatorID(umpire::resource::Device)
  #ifdef UMPIRE_ENABLE_CONST
     || alloc_id == axom::getUmpireResourceAllocatorID(umpire::resource::Constant)
  #endif
  )
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
template <typename T1, typename T2, int DIM, typename LArrayType, typename RArrayType>
bool operator==(const ArrayBase<T1, DIM, LArrayType>& lhs,
                const ArrayBase<T2, DIM, RArrayType>& rhs)
{
  static_assert(std::is_same<typename std::remove_const<T1>::type,
                             typename std::remove_const<T2>::type>::value,
                "Cannot compare Arrays of incompatible type");
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
    if(!(lhs.flatIdx(i) == rhs.flatIdx(i)))
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T1, typename T2, int DIM, typename LArrayType, typename RArrayType>
bool operator!=(const ArrayBase<T1, DIM, LArrayType>& lhs,
                const ArrayBase<T2, DIM, RArrayType>& rhs)
{
  return !(lhs == rhs);
}

namespace detail
{
// Takes the product of a parameter pack of T
template <typename T, int N>
AXOM_HOST_DEVICE T packProduct(const T (&arr)[N])
{
  T prod = 1;
  for(int idx = 0; idx < N; idx++)
  {
    prod *= arr[idx];
  }
  return prod;
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

/// \brief Indirection needed to dodge an MSVC compiler bug
template <typename... Args>
struct all_types_are_integral_impl : std::true_type
{ };

template <typename First, typename... Rest>
struct all_types_are_integral_impl<First, Rest...>
{
  static constexpr bool value = std::is_integral<First>::value &&
    all_types_are_integral_impl<Rest...>::value;
};

/// \brief Checks if all types in a parameter pack are integral
template <typename... Args>
struct all_types_are_integral
{
  static constexpr bool value = all_types_are_integral_impl<Args...>::value;
};

/*!
 * \brief Default-initializes the "new" segment of an array
 *
 * \param [inout] data The data to initialize
 * \param [in] pos The beginning of the subset of \a data that should be initialized
 * \param [in] len The length of the subset of \a data that is to be initialized
 * \param [in] alloc_id The allocator ID with which \a data was allocated
 * 
 * The final parameter is intended to be used via tag dispatch with std::is_default_constructible
 */
template <typename T>
void initializeInPlace(T* data,
                       const IndexType pos,
                       const IndexType len,
                       const int alloc_id,
                       std::true_type)
{
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE) && \
  defined(UMPIRE_ENABLE_DEVICE)
  if(alloc_id == getAllocatorID<MemorySpace::Device>())
  {
    // If we instantiated a fill kernel here it would require
    // that T's default ctor is device-annotated which is too
    // strict of a requirement, so we copy a buffer instead.
    T* tmp_buffer = new T[len]();
    axom::copy(data + pos, tmp_buffer, len * sizeof(T));
    delete[] tmp_buffer;
    return;
  }
#else
  AXOM_UNUSED_VAR(alloc_id);
#endif
  for(int ielem = 0; ielem < len; ielem++)
  {
    new(&data[pos + ielem]) T {};
  }
  // XL deviates from the standard in the above line, so we initialize directly
#ifdef __ibmxl__
  std::fill_n(data + pos, len, T {});
#endif
}

/// \overload
template <typename T>
void initializeInPlace(T*, const IndexType, const IndexType, const int, std::false_type)
{ }

template <typename T, int DIM>
struct ArrayTraits<ArrayViewProxy<T, DIM>>
{
  constexpr static bool is_view = true;

  template <int SliceDim>
  using Slice = ArrayViewProxy<T, DIM - SliceDim>;
};

/*!
 * \brief Proxy class for efficiently constructing slice ArrayViews by using
 *  the generic ArrayBase constructor in ArrayView. This avoids the cost of
 *  looking up the correct allocator ID from the other ArrayView constructors.
 */
template <typename T, int DIM>
class ArrayViewProxy : public ArrayBase<T, DIM, ArrayViewProxy<T, DIM>>
{
public:
  AXOM_HOST_DEVICE ArrayViewProxy(T* data,
                                  int allocatorID,
                                  const StackArray<IndexType, DIM>& shape)
    : ArrayBase<T, DIM, ArrayViewProxy<T, DIM>>(shape)
    , m_data(data)
    , m_allocatorID(allocatorID)
    , m_size(packProduct(shape.m_data))
  { }

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  inline AXOM_HOST_DEVICE IndexType size() const { return m_size; }

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  AXOM_HOST_DEVICE inline T* data() const { return m_data; }

  /// @}

  /*!
   * \brief Get the ID for the umpire allocator
   */
  int getAllocatorID() const { return m_allocatorID; }

private:
  T* const m_data;
  const int m_allocatorID;
  const IndexType m_size;
};

}  // namespace detail

} /* namespace axom */

#endif /* AXOM_ARRAYBASE_HPP_ */
