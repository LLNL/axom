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
#include "axom/core/execution/for_all.hpp"   // for for_all, *_EXEC

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
}

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
  AXOM_HOST_DEVICE T& operator()(Args... args)
  {
    const IndexType indices[] = {static_cast<IndexType>(args)...};
    const IndexType idx = numerics::dot_product(indices, m_strides.begin(), DIM);
    assert(inBounds(idx));
    return asDerived().data()[idx];
  }
  /// \overload
  template <typename... Args,
            typename SFINAE = typename std::enable_if<sizeof...(Args) == DIM>::type>
  AXOM_HOST_DEVICE RealConstT& operator()(Args... args) const
  {
    const IndexType indices[] = {static_cast<IndexType>(args)...};
    const IndexType idx = numerics::dot_product(indices, m_strides.begin(), DIM);
    assert(inBounds(idx));
    return asDerived().data()[idx];
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

  ArrayBase(IndexType = 0) { }

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
    if(!(lhs[i] == rhs[i]))
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

template <bool has_default_ctor>
struct OpInitBase
{
  template <typename T, typename ExecSpace>
  static void init(T* data, IndexType begin, IndexType end);
};

template <bool fill_on_host>
struct OpFillBase
{
  template <typename T, typename ExecSpace>
  static void fill(T* array, IndexType n, const T& value);
};

template <bool destroy_on_host>
struct OpDestroyBase
{
  template <typename T, typename ExecSpace>
  static void destroy(T* array, IndexType begin, IndexType end);
};

template <bool device>
struct OpMemmoveBase
{
  template <typename T>
  static void move(T* array,
                   IndexType src_begin,
                   IndexType src_end,
                   IndexType dst,
                   int allocId);
};

/*!
 * \brief Default-initializes the "new" segment of an array
 *
 * \param [inout] data The data to initialize
 * \param [in] begin The beginning of the subset of \a data that should be initialized
 * \param [in] end The end index (exclusive) of the subset of \a data that is to be initialized
 */
template <>
template <typename T, typename ExecSpace>
void OpInitBase<true>::init(T* data, IndexType begin, IndexType end)
{
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    IndexType len = end - begin;
    // If we instantiated a fill kernel here it would require
    // that T's default ctor is device-annotated which is too
    // strict of a requirement, so we copy a buffer instead.
    void* tmp_buffer = ::operator new(sizeof(T) * len);
    T* typed_buffer = static_cast<T*>(tmp_buffer);
    for(IndexType i = 0; i < len; ++i)
    {
      // We use placement-new to avoid calling destructors in the delete
      // statement below.
      T* inst = new(typed_buffer + i) T {};
    }
    axom::copy(data + begin, tmp_buffer, len * sizeof(T));
    ::operator delete(tmp_buffer);
    return;
  }
#endif
  for(int ielem = begin; ielem < end; ++ielem)
  {
    new(&data[ielem]) T {};
  }
  // XL deviates from the standard in the above line, so we initialize directly
#ifdef __ibmxl__
  IndexType len = end - begin;
  std::fill_n(data + begin, len, T {});
#endif
}

/*!
 * \brief Specialization for types T that aren't default-constructible
 */
template <>
template <typename T, typename ExecSpace>
void OpInitBase<false>::init(T*, const IndexType, const IndexType)
{ }

/*!
 * \brief Fills an array with objects of type T. This specialization is used
 *  for trivially-copyable objects, as well as non-trivially-copyable objects
 *  when the array is host-accessible.
 *
 * \param [inout] array the array to fill
 * \param [in] n the number of elements to fill the array with
 * \param [in] value the value to set each array element to
 */
template <>
template <typename T, typename ExecSpace>
void OpFillBase<false>::fill(T* array, IndexType n, const T& value)
{
  for_all<ExecSpace>(
    n,
    AXOM_LAMBDA(IndexType i) { array[i] = value; });
}

/*!
 * \brief Fills an array with objects of type T. This specialization is used
 *  for device memory when the type is not trivially-copyable.
 *
 * \param [inout] array the array to fill
 * \param [in] n the number of elements to fill the array with
 * \param [in] value the value to set each array element to
 */
template <>
template <typename T, typename ExecSpace>
void OpFillBase<true>::fill(T* array, IndexType n, const T& value)
{
  void* buffer = ::operator new(sizeof(T) * n);
  T* typed_buffer = static_cast<T*>(buffer);
  // If we instantiated a fill kernel here it would require
  // that T's copy ctor is device-annotated which is too
  // strict of a requirement, so we copy a buffer instead.
  std::uninitialized_fill_n(typed_buffer, n, value);
  axom::copy(array, typed_buffer, sizeof(T) * n);
  ::operator delete(buffer);
}

/*!
 * \brief Calls the destructor on a range of typed elements in the array. This
 *  specialization is used for trivial objects, or non-trivially-destructible
 *  objects when the array is host-accessible.
 *
 * \param [inout] array the array with elements to destroy
 * \param [in] begin the start index of the range of elements to destroy
 * \param [in] value one past the end index of the range of elements to destroy
 */
template <>
template <typename T, typename ExecSpace>
void OpDestroyBase<false>::destroy(T* array, IndexType begin, IndexType end)
{
  if(!std::is_trivial<T>::value)
  {
    for_all<ExecSpace>(
      begin,
      end,
      AXOM_LAMBDA(IndexType i) { array[i].~T(); });
  }
}

/*!
 * \brief Calls the destructor on a range of typed elements in the array. This
 *  specialization is used for non-trivially-destructible objects when the
 *  array is in device space.
 *
 * \param [inout] array the array with elements to destroy
 * \param [in] begin the start index of the range of elements to destroy
 * \param [in] value one past the end index of the range of elements to destroy
 */
template <>
template <typename T, typename ExecSpace>
void OpDestroyBase<true>::destroy(T* array, IndexType begin, IndexType end)
{
  IndexType n = end - begin;
  void* buffer = ::operator new(sizeof(T) * n);
  T* typed_buffer = static_cast<T*>(buffer);
  axom::copy(typed_buffer, array, sizeof(T) * n);
  for(int i = begin; i < end; ++i)
  {
    typed_buffer[i].~T();
  }
  axom::copy(array, typed_buffer, sizeof(T) * n);
  ::operator delete(buffer);
}

/*!
 * \brief Moves a range of data in the array. This is called on host-accessible
 *  memory.
 *
 * \param [inout] array the array with elements to move
 * \param [in] src_begin the start index of the source range
 * \param [in] src_end the end index of the source range, exclusive
 * \param [in] dst the destination index of the range of elements
 */
template <>
template <typename T>
void OpMemmoveBase<false>::move(T* array,
                                IndexType src_begin,
                                IndexType src_end,
                                IndexType dst,
                                int allocId)
{
  AXOM_UNUSED_VAR(allocId);
  std::memmove(array + dst, array + src_begin, (src_end - src_begin) * sizeof(T));
}

/*!
 * \brief Moves a range of data in the array. This is called on device-only
 *  memory.
 *
 * \param [inout] array the array with elements to move
 * \param [in] src_begin the start index of the source range
 * \param [in] src_end the end index of the source range, exclusive
 * \param [in] dst the destination index of the range of elements
 * \param [in] allocId the allocator ID of the array
 */
template <>
template <typename T>
void OpMemmoveBase<true>::move(T* array,
                               IndexType src_begin,
                               IndexType src_end,
                               IndexType dst,
                               int allocId)
{
  // Since this memory is on the device-side, we copy it to a temporary buffer
  // first.
  IndexType nelems = src_end - src_begin;
  T* tmp_buf = axom::allocate<T>(nelems, allocId);
  axom::copy(tmp_buf, array + src_begin, nelems * sizeof(T));
  axom::copy(array + dst, tmp_buf, nelems * sizeof(T));
  axom::deallocate(tmp_buf);
}

template <typename T, MemorySpace SPACE>
struct ArrayOps
{
private:
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
  constexpr static bool IsDevice = (SPACE == MemorySpace::Device);

  using ExecSpace =
    typename std::conditional<IsDevice, axom::CUDA_EXEC<256>, axom::SEQ_EXEC>::type;

#else
  constexpr static bool IsDevice = false;

  using ExecSpace = axom::SEQ_EXEC;
#endif

  using InitBase = OpInitBase<std::is_default_constructible<T>::value>;
  // TODO: the below should be std::is_trivially_copyable and
  // std::is_trivially_destructible; however these aren't available
  // in gcc 4.9.3
  using FillBase = OpFillBase<!std::is_trivial<T>::value && IsDevice>;
  using DestroyBase = OpDestroyBase<!std::is_trivial<T>::value && IsDevice>;
  using MoveBase = OpMemmoveBase<IsDevice>;

public:
  static void init(T* array, IndexType begin, IndexType end, int allocId)
  {
    AXOM_UNUSED_VAR(allocId);
    InitBase::template init<T, ExecSpace>(array, begin, end);
  }

  static void fill(T* array, IndexType n, int allocId, const T& value)
  {
    AXOM_UNUSED_VAR(allocId);
    FillBase::template fill<T, ExecSpace>(array, n, value);
  }

  static void destroy(T* array, IndexType begin, IndexType end, int allocId)
  {
    AXOM_UNUSED_VAR(allocId);
    if(!array || end <= begin)
    {
      return;
    }
    DestroyBase::template destroy<T, ExecSpace>(array, begin, end);
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
    MoveBase::template move<T>(array, src_begin, src_end, dst, allocId);
  }
};

template <typename T>
struct ArrayOps<T, MemorySpace::Dynamic>
{
private:
  using InitBase = OpInitBase<std::is_default_constructible<T>::value>;
  using FillBase = OpFillBase<false>;
  using DestroyBase = OpDestroyBase<false>;
  using MoveBase = OpMemmoveBase<false>;

public:
  static void init(T* array, IndexType begin, IndexType end, int allocId)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      ArrayOps<T, MemorySpace::Device>::init(array, begin, end, allocId);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    InitBase::template init<T, axom::SEQ_EXEC>(array, begin, end);
  }

  static void fill(T* array, IndexType n, int allocId, const T& value)
  {
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      ArrayOps<T, MemorySpace::Device>::fill(array, n, allocId, value);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    FillBase::template fill<T, axom::SEQ_EXEC>(array, n, value);
  }

  static void destroy(T* array, IndexType begin, IndexType end, int allocId)
  {
    if(!array || end <= begin)
    {
      return;
    }
#if defined(__CUDACC__) && defined(AXOM_USE_UMPIRE)
    MemorySpace space = getAllocatorSpace(allocId);

    if(space == MemorySpace::Device)
    {
      ArrayOps<T, MemorySpace::Device>::destroy(array, begin, end, allocId);
      return;
    }
#else
    AXOM_UNUSED_VAR(allocId);
#endif
    DestroyBase::template destroy<T, axom::SEQ_EXEC>(array, begin, end);
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
      ArrayOps<T, MemorySpace::Device>::move(array, src_begin, src_end, dst, allocId);
      return;
    }
#endif
    MoveBase::template move<T>(array, src_begin, src_end, dst, allocId);
  }
};

}  // namespace detail

} /* namespace axom */

#endif /* AXOM_ARRAYBASE_HPP_ */
