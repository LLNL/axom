// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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
#include "axom/core/execution/for_all.hpp"   // for for_all, *_EXEC

// C/C++ includes
#include <iostream>  // for std::cerr and std::ostream
#include <numeric>   // for std::accumulate

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM, typename ArrayType>
class ArrayBase;

namespace detail
{
template <typename ArrayType>
struct ArrayTraits;

template <typename T, int DIM, typename BaseArray>
class ArraySubslice;

template <typename T, int SliceDim, typename BaseArray>
struct SubsliceProxy
{
  using type = ArraySubslice<T, SliceDim, BaseArray>;
};

// The below specializations ensure that subslices of subslices refer to the
// original base array type, e.g.:
//   Array<int, 3> arr;
//   auto slice_1 = arr[i]
//      -> returns ArraySubslice<int, 2, Array<int, 3>>
//   auto slice_2 = slice_1[j]
//      -> returns ArraySubslice<int, 1, Array<int,3>> instead of
//         ArraySubslice<int, 1, ArraySubslice<int, 2, Array<int, 3>>>
template <typename T, int SliceDim, int OldSliceDim, typename BaseArray>
struct SubsliceProxy<T, SliceDim, ArraySubslice<T, OldSliceDim, BaseArray>>
{
  using type = ArraySubslice<T, SliceDim, BaseArray>;
};

template <typename T, int SliceDim, int OldSliceDim, typename BaseArray>
struct SubsliceProxy<T, SliceDim, const ArraySubslice<T, OldSliceDim, BaseArray>>
{
  using type = ArraySubslice<T, SliceDim, BaseArray>;
};

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

  template <int IdxDim>
  using SliceType = typename std::conditional<
    DIM == IdxDim,
    T&,
    typename detail::SubsliceProxy<T, DIM - IdxDim, ArrayType>::type>::type;

  template <int IdxDim>
  using ConstSliceType = typename std::conditional<
    DIM == IdxDim,
    RealConstT&,
    typename detail::SubsliceProxy<T, DIM - IdxDim, const ArrayType>::type>::type;

  constexpr static int Dims = DIM;

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
    return asDerived().sliceImpl(idx);
  }

  /// \overload
  template <int UDim>
  AXOM_HOST_DEVICE ConstSliceType<UDim> operator[](
    const StackArray<IndexType, UDim>& idx) const
  {
    static_assert(UDim <= DIM,
                  "Index dimensions cannot be larger than array dimensions");
    return asDerived().sliceImpl(idx);
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
  AXOM_HOST_DEVICE T& flatIndex(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIndex(const IndexType idx) const
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
    return SliceType<UDim>(&asDerived(), idx);
  }

  /// \overload
  template <int UDim>
  AXOM_HOST_DEVICE ConstSliceType<UDim> sliceImpl(
    const StackArray<IndexType, UDim>& idx) const
  {
    return ConstSliceType<UDim>(&asDerived(), idx);
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
  AXOM_HOST_DEVICE T& flatIndex(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIndex(const IndexType idx) const
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
    if(!(lhs.flatIndex(i) == rhs.flatIndex(i)))
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

/// \brief Takes the last N elements from an array.
template <int N, typename T, int DIM>
AXOM_HOST_DEVICE StackArray<T, N> takeLastElems(const StackArray<T, DIM>& arr)
{
  static_assert(N <= DIM,
                "Attempting to take more elements than the array holds");
  StackArray<T, N> ret;
  for(int i = 0; i < N; i++)
  {
    ret[i] = arr[i + DIM - N];
  }
  return ret;
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

#if defined(__GLIBCXX__) && !defined(_GLIBCXX_USE_CXX11_ABI)
// Some type traits we need aren't fully-implemented in libstdc++ for GCC <5.
// We thus check for the macro _GLIBCXX_USE_CXX11_ABI, which should only be
// defined for GCC 5+; if not defined, we'll use the below fallbacks.
template <typename T>
using HasTrivialDefaultCtor = ::std::has_trivial_default_constructor<T>;
template <typename T>
using HasTrivialCopyCtor = ::std::has_trivial_copy_constructor<T>;
#else
template <typename T>
using HasTrivialDefaultCtor = ::std::is_trivially_default_constructible<T>;
template <typename T>
using HasTrivialCopyCtor = ::std::is_trivially_copy_constructible<T>;
#endif

template <typename T, bool DeviceOps>
struct ArrayOpsBase;

template <typename T>
struct ArrayOpsBase<T, false>
{
  using DefaultCtorTag = std::is_default_constructible<T>;

  /*!
   * \brief Helper for default-initializing the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] end The end index (exclusive) of the subset of \a data that is to be initialized
   * \note Specialization for when T is default-constructible.
   */
  template <typename U = void>
  static void init_impl(T* data, IndexType begin, IndexType nelems, std::true_type)
  {
    for(IndexType i = 0; i < nelems; ++i)
    {
      new(data + i + begin) T();
    }
  }

  /*!
   * \overload
   * \note Specialization for when T is not default-constructible.
   */
  template <typename U = void>
  static void init_impl(T*, IndexType, IndexType, std::false_type)
  { }

  /*!
   * \brief Default-initializes the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] end The end index (exclusive) of the subset of \a data that is to be initialized
   * \note Specialization for when T is default-constructible.
   */
  static void init(T* data, IndexType begin, IndexType nelems)
  {
    init_impl(data, begin, nelems, DefaultCtorTag {});
  }

  /*!
   * \brief Fills an uninitialized array with objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index in the array to begin filling elements at
   * \param [in] nelems the number of elements to fill the array with
   * \param [in] value the value to set each array element to
   */
  static void fill(T* array, IndexType begin, IndexType nelems, const T& value)
  {
    std::uninitialized_fill_n(array + begin, nelems, value);
  }

  /*!
   * \brief Fills an uninitialized array with a range of objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index at which to begin placing elements
   * \param [in] nelems the number of elements in the range to fill the array with
   * \param [in] values the values to set each array element to
   */
  static void fill_range(T* array, IndexType begin, IndexType nelems, const T* values)
  {
    std::uninitialized_copy(values, values + nelems, array + begin);
  }

  /*!
   * \brief Constructs a new element in uninitialized memory.
   *
   * \param [inout] array the array to construct in
   * \param [in] i the array index in which to construct the new object
   * \param [in] args the arguments to forward to constructor of the element.
   */
  template <typename... Args>
  static void emplace(T* array, IndexType i, Args&&... args)
  {
    new(array + i) T(std::forward<Args>(args)...);
  }

  /*!
   * \brief Calls the destructor on a range of typed elements in the array.
   *
   * \param [inout] array the array with elements to destroy
   * \param [in] begin the start index of the range of elements to destroy
   * \param [in] value one past the end index of the range of elements to destroy
   */
  static void destroy(T* array, IndexType begin, IndexType nelems)
  {
    if(!std::is_trivially_destructible<T>::value)
    {
      for(IndexType i = 0; i < nelems; i++)
      {
        array[i + begin].~T();
      }
    }
  }

  /*!
   * \brief Moves a range of data in the array.
   *
   * \param [inout] array the array with elements to move
   * \param [in] src_begin the start index of the source range
   * \param [in] src_end the end index of the source range, exclusive
   * \param [in] dst the destination index of the range of elements
   */
  static void move(T* array, IndexType src_begin, IndexType src_end, IndexType dst)
  {
    if(src_begin < dst)
    {
      IndexType dst_last = dst + src_end - src_begin;
      std::move_backward(array + src_begin, array + src_end, array + dst_last);
    }
    else if(src_begin > dst)
    {
      std::move(array + src_begin, array + src_end, array + dst);
    }
  }
};

#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
template <typename T>
struct ArrayOpsBase<T, true>
{
  using ExecSpace = axom::CUDA_EXEC<256>;

  static constexpr bool InitOnDevice = HasTrivialDefaultCtor<T>::value;
  static constexpr bool DestroyOnHost = !std::is_trivially_destructible<T>::value;
  static constexpr bool DefaultCtor = std::is_default_constructible<T>::value;

  struct TrivialDefaultCtorTag
  { };
  struct NontrivialDefaultCtorTag
  { };
  struct NoDefaultCtorTag
  { };

  /*!
   * \brief Helper for default-initialization of a range of elements.
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
   * \note Specialization for when T is nontrivially default-constructible.
   */
  template <typename U = void>
  static void init_impl(T* data,
                        IndexType begin,
                        IndexType nelems,
                        NontrivialDefaultCtorTag)
  {
    // If we instantiated a fill kernel here it would require
    // that T's default ctor is device-annotated which is too
    // strict of a requirement, so we copy a buffer instead.
    void* tmp_buffer = ::operator new(sizeof(T) * nelems);
    T* typed_buffer = static_cast<T*>(tmp_buffer);
    for(IndexType i = 0; i < nelems; ++i)
    {
      // We use placement-new to avoid calling destructors in the delete
      // statement below.
      new(typed_buffer + i) T();
    }
    axom::copy(data + begin, tmp_buffer, nelems * sizeof(T));
    ::operator delete(tmp_buffer);
  }

  /*!
   * \overload
   * \note Specialization for when T is trivially default-constructible.
   */
  template <typename U = void>
  static void init_impl(T* data,
                        IndexType begin,
                        IndexType nelems,
                        TrivialDefaultCtorTag)
  {
    for_all<ExecSpace>(
      begin,
      begin + nelems,
      AXOM_LAMBDA(IndexType i) { new(&data[i]) T(); });
  }

  /*!
   * \overload
   * \note Specialization for when T is not default-constructible.
   */
  template <typename U = void>
  static void init_impl(T*, IndexType, IndexType, NoDefaultCtorTag)
  { }

  /*!
   * \brief Default-initializes the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
   */
  static void init(T* data, IndexType begin, IndexType nelems)
  {
    using InitSelectTag =
      typename std::conditional<InitOnDevice,
                                TrivialDefaultCtorTag,
                                NontrivialDefaultCtorTag>::type;
    using InitTag =
      typename std::conditional<DefaultCtor, InitSelectTag, NoDefaultCtorTag>::type;

    init_impl(data, begin, nelems, InitTag {});
  }

  /*!
   * \brief Helper for filling an uninitialized array with objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index in the array to begin filling elements at
   * \param [in] nelems the number of elements to fill the array with
   * \param [in] value the value to set each array element to
   * \note Specialization for when T is not trivially-copyable.
   */
  template <typename U = void>
  static void fill_impl(T* array,
                        IndexType begin,
                        IndexType nelems,
                        const T& value,
                        std::false_type)
  {
    void* buffer = ::operator new(sizeof(T) * nelems);
    T* typed_buffer = static_cast<T*>(buffer);
    // If we instantiated a fill kernel here it would require
    // that T's copy ctor is device-annotated which is too
    // strict of a requirement, so we copy a buffer instead.
    std::uninitialized_fill_n(typed_buffer, nelems, value);
    axom::copy(array + begin, typed_buffer, sizeof(T) * nelems);
    ::operator delete(buffer);
  }

  /*!
   * \overload
   * \note Specialization for when T is trivially-copyable.
   */
  template <typename U = void>
  static void fill_impl(T* array,
                        IndexType begin,
                        IndexType nelems,
                        const T& value,
                        std::true_type)
  {
    for_all<ExecSpace>(
      nelems,
      AXOM_LAMBDA(IndexType i) { array[i + begin] = value; });
  }

  /*!
   * \brief Fills an uninitialized array with objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index in the array to begin filling elements at
   * \param [in] nelems the number of elements to fill the array with
   * \param [in] value the value to set each array element to
   */
  static void fill(T* array, IndexType begin, IndexType nelems, const T& value)
  {
    fill_impl(array, begin, nelems, value, HasTrivialCopyCtor<T> {});
  }

  /*!
   * \brief Fills an uninitialized array with a range of objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index at which to begin placing elements
   * \param [in] nelems the number of elements in the range to fill the array with
   * \param [in] values the values to set each array element to
   */
  static void fill_range(T* array, IndexType begin, IndexType nelems, const T* values)
  {
    axom::copy(array + begin, values, sizeof(T) * nelems);
  }

  /*!
   * \brief Constructs a new element in uninitialized memory.
   *
   * \param [inout] array the array to construct in
   * \param [in] i the array index in which to construct the new object
   * \param [in] args the arguments to forward to constructor of the element.
   */
  template <typename... Args>
  static void emplace(T* array, IndexType i, Args&&... args)
  {
    // Similar to fill(), except we can allocate stack memory and placement-new
    // the object with a move constructor.
    alignas(T) axom::uint8 host_buf[sizeof(T)];
    T* host_obj = ::new(&host_buf) T(std::forward<Args>(args)...);
    axom::copy(array + i, host_obj, sizeof(T));
  }

  /*!
   * \brief Calls the destructor on a range of typed elements in the array.
   *
   * \param [inout] array the array with elements to destroy
   * \param [in] begin the start index of the range of elements to destroy
   * \param [in] nelems the number of elements to destroy
   */
  static void destroy(T* array, IndexType begin, IndexType nelems)
  {
    if(DestroyOnHost)
    {
      void* buffer = ::operator new(sizeof(T) * nelems);
      T* typed_buffer = static_cast<T*>(buffer);
      axom::copy(typed_buffer, array + begin, sizeof(T) * nelems);
      for(int i = 0; i < nelems; ++i)
      {
        typed_buffer[i].~T();
      }
      axom::copy(array + begin, typed_buffer, sizeof(T) * nelems);
      ::operator delete(buffer);
    }
  }

  /*!
   * \brief Moves a range of data in the array.
   *
   * \param [inout] array the array with elements to move
   * \param [in] src_begin the start index of the source range
   * \param [in] src_end the end index of the source range, exclusive
   * \param [in] dst the destination index of the range of elements
   */
  static void move(T* array, IndexType src_begin, IndexType src_end, IndexType dst)
  {
    // Since this memory is on the device-side, we copy it to a temporary buffer
    // first.
    IndexType nelems = src_end - src_begin;
    T* tmp_buf =
      axom::allocate<T>(nelems, axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(tmp_buf, array + src_begin, nelems * sizeof(T));
    axom::copy(array + dst, tmp_buf, nelems * sizeof(T));
    axom::deallocate(tmp_buf);
  }
};
#endif

template <typename T, MemorySpace SPACE>
struct ArrayOps
{
private:
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
  constexpr static bool IsDevice = (SPACE == MemorySpace::Device);
#else
  constexpr static bool IsDevice = false;
#endif

  using Base = ArrayOpsBase<T, IsDevice>;

public:
  static void init(T* array, IndexType begin, IndexType nelems, int allocId)
  {
    AXOM_UNUSED_VAR(allocId);
    Base::init(array, begin, nelems);
  }

  static void fill(T* array,
                   IndexType begin,
                   IndexType nelems,
                   int allocId,
                   const T& value)
  {
    AXOM_UNUSED_VAR(allocId);
    Base::fill(array, begin, nelems, value);
  }

  static void fill_range(T* array,
                         IndexType begin,
                         IndexType nelems,
                         int allocId,
                         const T* values)
  {
    AXOM_UNUSED_VAR(allocId);
    Base::fill_range(array, begin, nelems, values);
  }

  static void destroy(T* array, IndexType begin, IndexType nelems, int allocId)
  {
    AXOM_UNUSED_VAR(allocId);
    if(nelems == 0)
    {
      return;
    }
    Base::destroy(array, begin, nelems);
  }

  static void move(T* array,
                   IndexType src_begin,
                   IndexType src_end,
                   IndexType dst,
                   int allocId)
  {
    AXOM_UNUSED_VAR(allocId);
    if(src_begin >= src_end)
    {
      return;
    }
    Base::move(array, src_begin, src_end, dst);
  }

  template <typename... Args>
  static void emplace(T* array, IndexType dst, IndexType allocId, Args&&... args)
  {
    AXOM_UNUSED_VAR(allocId);
    Base::emplace(array, dst, std::forward<Args>(args)...);
  }
};

template <typename T>
struct ArrayOps<T, MemorySpace::Dynamic>
{
private:
  using Base = ArrayOpsBase<T, false>;
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
  using BaseDevice = ArrayOpsBase<T, true>;
#endif

public:
  static void init(T* array, IndexType begin, IndexType nelems, int allocId)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::init(array, begin, nelems);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    Base::init(array, begin, nelems);
  }

  static void fill(T* array,
                   IndexType begin,
                   IndexType nelems,
                   int allocId,
                   const T& value)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::fill(array, begin, nelems, value);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    Base::fill(array, begin, nelems, value);
  }

  static void fill_range(T* array,
                         IndexType begin,
                         IndexType nelems,
                         int allocId,
                         const T* values)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::fill_range(array, begin, nelems, values);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    AXOM_UNUSED_VAR(allocId);
    Base::fill_range(array, begin, nelems, values);
  }

  static void destroy(T* array, IndexType begin, IndexType nelems, int allocId)
  {
    if(nelems == 0)
    {
      return;
    }
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::destroy(array, begin, nelems);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    Base::destroy(array, begin, nelems);
  }

  static void move(T* array,
                   IndexType src_begin,
                   IndexType src_end,
                   IndexType dst,
                   int allocId)
  {
    if(src_begin >= src_end)
    {
      return;
    }
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::move(array, src_begin, src_end, dst);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    Base::move(array, src_begin, src_end, dst);
  }

  template <typename... Args>
  static void emplace(T* array, IndexType dst, IndexType allocId, Args&&... args)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      BaseDevice::emplace(array, dst, std::forward<Args>(args)...);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    Base::emplace(array, dst, std::forward<Args>(args)...);
  }
};

template <typename T, int SliceDim, typename BaseArray>
struct ArrayTraits<ArraySubslice<T, SliceDim, BaseArray>>
{
  constexpr static bool is_view = true;
};

/*!
 * \brief Proxy class for efficiently constructing slice ArrayViews by using
 *  the generic ArrayBase constructor in ArrayView. This avoids the cost of
 *  looking up the correct allocator ID from the other ArrayView constructors.
 */
template <typename T, int SliceDim, typename BaseArray>
class ArraySubslice
  : public ArrayBase<T, SliceDim, ArraySubslice<T, SliceDim, BaseArray>>
{
  using BaseClass = ArrayBase<T, SliceDim, ArraySubslice<T, SliceDim, BaseArray>>;
  using RefType = typename std::conditional<std::is_const<BaseArray>::value,
                                            typename BaseClass::RealConstT&,
                                            T&>::type;

  constexpr static int OrigDims = BaseArray::Dims;
  constexpr static int NumIndices = OrigDims - SliceDim;

  template <int UDim>
  using SliceType =
    typename std::conditional<UDim == SliceDim,
                              RefType,
                              ArraySubslice<T, SliceDim - UDim, BaseArray>>::type;

public:
  AXOM_HOST_DEVICE ArraySubslice(BaseArray* array,
                                 const StackArray<IndexType, NumIndices>& idxs)
    : BaseClass(detail::takeLastElems<SliceDim>(array->shape()))
    , m_array(array)
    , m_inds(idxs)
  { }

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  inline AXOM_HOST_DEVICE IndexType size() const
  {
    IndexType size = 1;
    for(int idx = NumIndices; idx < OrigDims; idx++)
    {
      size *= m_array->shape()[idx];
    }
    return size;
  }

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  AXOM_HOST_DEVICE inline T* data() const
  {
    IndexType offset = numerics::dot_product(m_inds.begin(),
                                             m_array->strides().begin(),
                                             NumIndices);
    return m_array->data() + offset;
  }

  /// @}

  /*!
   * \brief Get the ID for the umpire allocator
   */
  int getAllocatorID() const { return m_array->getAllocatorID(); }

protected:
  friend BaseClass;

  template <int UDim>
  AXOM_HOST_DEVICE SliceType<UDim> sliceImpl(
    const StackArray<IndexType, UDim>& idx) const
  {
    StackArray<IndexType, UDim + NumIndices> full_inds;
    for(int i = 0; i < NumIndices; i++)
    {
      full_inds[i] = m_inds[i];
    }
    for(int i = 0; i < UDim; i++)
    {
      full_inds[i + NumIndices] = idx[i];
    }
    return (*m_array)[full_inds];
  }

private:
  BaseArray* const m_array;
  const StackArray<IndexType, NumIndices> m_inds;
};

}  // namespace detail

} /* namespace axom */

#endif /* AXOM_ARRAYBASE_HPP_ */
