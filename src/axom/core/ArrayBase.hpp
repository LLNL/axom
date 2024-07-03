// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAYBASE_HPP_
#define AXOM_ARRAYBASE_HPP_

#include "axom/config.hpp"                    // for compile-time defines
#include "axom/core/Macros.hpp"               // for axom macros
#include "axom/core/MDMapping.hpp"            // for index conversion
#include "axom/core/memory_management.hpp"    // for memory allocation functions
#include "axom/core/utilities/Utilities.hpp"  // for processAbort()
#include "axom/core/Types.hpp"                // for IndexType definition
#include "axom/core/StackArray.hpp"
#include "axom/core/numerics/matvecops.hpp"  // for dot_product
#include "axom/core/execution/for_all.hpp"   // for for_all, *_EXEC

// C/C++ includes
#include <iostream>  // for std::cerr and std::ostream
#include <numeric>   // for std::accumulate
#if defined(__GLIBCXX__) && !defined(_GLIBCXX_USE_CXX11_ABI)
  #error \
    "GNU libstdc++ versions less than 5 do not fully support C++11 and are unsupported by Axom."
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

  //! @brief Construct row-major, unitnitialized array.
  AXOM_HOST_DEVICE ArrayBase() : m_shape(), m_mapping(ArrayStrideOrder::ROW)
  {
    updateStrides();
  }

  /*!
   * \brief Parameterized constructor that sets up the array shape.
   *
   * \param [in] shape Array size in each direction.
   * \param [in] min_stride The minimum stride between two consecutive
   *  elements in row-major order.
   */
  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, DIM>& shape,
                             int min_stride = 1)
    : m_shape {shape}
    , m_mapping(shape, ArrayStrideOrder::ROW, min_stride)
    , m_minStride(m_mapping.fastestStrideLength())
  { }

  /*!
   * \brief Parameterized constructor that sets up the array shape,
   * with an MDMapping to specify data ordering.
   *
   * \param [in] shape Array size in each direction.
   * \param [in] mapping Model mapper, specifying
   *   the array stride order and minimum stride.
   *
   * The object is constructed with the given shape,
   * not the partial shape information in \c mapping.
   */
  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, DIM>& shape,
                             const MDMapping<DIM>& mapping)
    : m_shape {shape}
    , m_mapping(shape, mapping.slowestDirs(), mapping.fastestStrideLength())
    , m_minStride(m_mapping.fastestStrideLength())
  { }

  /*!
   * \brief Parameterized constructor that sets up the array shape and stride.
   *
   * \param [in] shape Array size in each direction.
   * \param [in] stride Array stride for each direction.
   */
  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, DIM>& shape,
                             const StackArray<IndexType, DIM>& stride)
    : m_shape {shape}
  {
    m_mapping.initializeStrides(stride, ArrayStrideOrder::ROW);
    m_minStride = m_mapping.fastestStrideLength();
    validateShapeAndStride(shape, stride);
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
    : m_shape(other.shape())
    , m_mapping(other.mapping())
    , m_minStride(m_mapping.fastestStrideLength())
  { }

  /// \overload
  template <typename OtherArrayType>
  ArrayBase(
    const ArrayBase<const typename std::remove_const<T>::type, DIM, OtherArrayType>& other)
    : m_shape(other.shape())
    , m_mapping(other.mapping())
    , m_minStride(m_mapping.fastestStrideLength())
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
   * \pre 0 <= args[i] < shape()[i] for i in [0, sizeof...(Args))
   */
  template <typename... Args>
  AXOM_HOST_DEVICE SliceType<sizeof...(Args)> operator()(Args... args)
  {
    static_assert(sizeof...(Args) <= DIM,
                  "Index dimensions different from array dimensions");
    constexpr int UDim = sizeof...(Args);
    const StackArray<IndexType, UDim> indices {{static_cast<IndexType>(args)...}};
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
   * \pre 0 <= idx < shape()[0]
   */
  AXOM_HOST_DEVICE SliceType<1> operator[](const IndexType idx)
  {
    const StackArray<IndexType, 1> slice {{idx}};
    return (*this)[slice];
  }

  /// \overload
  AXOM_HOST_DEVICE ConstSliceType<1> operator[](const IndexType idx) const
  {
    const StackArray<IndexType, 1> slice {{idx}};
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
   * \pre 0 <= args[i] < shape()[i] for i in [0, UDim)
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
   * \note equivalent to *(array.data() + idx * minStride()).
   *
   * \pre 0 <= idx < asDerived().size()
   */
  AXOM_HOST_DEVICE T& flatIndex(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * asDerived().minStride()];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIndex(const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * asDerived().minStride()];
  }
  /// @}

  /// \brief Swaps two ArrayBases
  void swap(ArrayBase& other)
  {
    std::swap(m_shape, other.m_shape);
    std::swap(m_mapping, other.m_mapping);
    std::swap(m_minStride, other.m_minStride);
  }

  /// \brief Returns the dimensions of the Array
  AXOM_HOST_DEVICE const StackArray<IndexType, DIM>& shape() const
  {
    return m_shape;
  }

  /// \brief Returns the multidimensional mapping for the Array
  AXOM_HOST_DEVICE const MDMapping<DIM>& mapping() const { return m_mapping; }

  /*!
   * \brief Returns the memory strides of the Array.
   */
  AXOM_HOST_DEVICE const StackArray<IndexType, DIM>& strides() const
  {
    return m_mapping.strides();
  }

  /*!
   * \brief Returns the minimum stride between adjacent items.
   */
  AXOM_HOST_DEVICE inline IndexType minStride() const { return m_minStride; }

protected:
  /// \brief Set the shape
  AXOM_HOST_DEVICE void setShape(const StackArray<IndexType, DIM>& shape_)
  {
#ifndef NDEBUG
    for(auto s : shape_)
    {
      assert(s >= 0);
    }
#endif

    m_shape = shape_;
    updateStrides();
  }

  /// \brief Set the shape and stride
  AXOM_HOST_DEVICE void setShapeAndStride(const StackArray<IndexType, DIM>& shape,
                                          const StackArray<IndexType, DIM>& stride)
  {
#ifdef AXOM_DEBUG
    validateShapeAndStride(shape, stride);
#endif
    m_shape = shape;
    m_mapping.initializeStrides(stride);
    m_minStride = m_mapping.fastestStrideLength();
  }

  /*!
   * \brief Returns the minimum "chunk size" that should be allocated
   * For example, 2 would be the chunk size of a 2D array whose second dimension is of size 2.
   * This is used when resizing/reallocating; it wouldn't make sense to have a
   * capacity of 3 in the array described above.
   */
  IndexType blockSize() const
  {
    auto slowestDir = m_mapping.slowestDirs()[0];
    return m_mapping.strides()[slowestDir];
  }

  /*!
   * \brief Updates the internal striding information to a row-major format
   * Intended to be called after shape is updated.
   */
  AXOM_HOST_DEVICE void updateStrides(int min_stride = 1)
  {
    // Update m_mapping strides while preserving stride order.
    m_mapping.initializeShape(m_shape, m_mapping.slowestDirs(), min_stride);
    m_minStride = m_mapping.fastestStrideLength();
  }

  /*!
   * \brief Updates the internal dimensions and striding based on the insertion
   *  of a range of elements.
   *  Intended to be called after the insertion of a multidimensional subslice.
   */
  void updateShapeOnInsert(const StackArray<IndexType, DIM>& range_shape)
  {
#ifdef AXOM_DEBUG
    if(!std::equal(m_shape.begin() + 1, m_shape.end(), range_shape.begin() + 1))
    {
      std::cerr << "Cannot append a multidimensional array of incorrect shape.";
      utilities::processAbort();
    }
#endif
    // First update the dimensions - we're adding only to the leading dimension
    m_shape[0] += range_shape[0];
    updateStrides();
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

  //// \brief Memory offset to get to the given multidimensional index.
  AXOM_HOST_DEVICE IndexType offset(const StackArray<IndexType, DIM>& idx) const
  {
    return m_mapping.toFlatIndex(idx);
  }

  /*!
   * \brief Returns the size of the range of memory in which elements are
   *  located. This is equivalent to size() * minStride().
   *
   *  offset() will return a value between [0, memorySize()).
   */
  AXOM_HOST_DEVICE IndexType memorySize() const
  {
    auto slowestDir = m_mapping.slowestDirs()[0];
    return m_mapping.strides()[slowestDir] * m_shape[slowestDir];
  }

  /*!
    \brief Memory offset to a slice at the given lower-dimensional index.

    Allowed only for row-major arrays.

    @pre mapping().getStrideOrder() & ArrayStrideOrder::ROW == true
  */
  template <int UDim>
  AXOM_HOST_DEVICE IndexType offset(const StackArray<IndexType, UDim>& idx) const
  {
    static_assert(UDim <= DIM,
                  "Index dimensions cannot be larger than array dimensions");
    assert(mapping().getStrideOrder() & ArrayStrideOrder::ROW);
    return numerics::dot_product(idx.begin(), m_mapping.strides().begin(), UDim);
  }

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if idx is within bounds */
  AXOM_HOST_DEVICE inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < memorySize();
  }

  /// \brief Checks a shape and stride array for correct bounds.
  AXOM_HOST_DEVICE inline void validateShapeAndStride(
    const StackArray<IndexType, DIM>& shape,
    const StackArray<IndexType, DIM>& stride)
  {
    int sorted_dims[DIM];
    for(int dim = 0; dim < DIM; dim++)
    {
      sorted_dims[dim] = dim;
    }
    // Sort the dimensions by stride.
    axom::utilities::insertionSort(sorted_dims,
                                   DIM,
                                   [&](int dim_a, int dim_b) -> bool {
                                     return stride[dim_a] < stride[dim_b];
                                   });
// Work from the smallest-strided dimension to the largest-strided.
#ifndef NDEBUG
    for(int dim = 0; dim < DIM - 1; dim++)
    {
      const int& minor_dim = sorted_dims[dim];
      const int& major_dim = sorted_dims[dim + 1];
      assert(stride[major_dim] >= stride[minor_dim] * shape[minor_dim]);
      assert(stride[major_dim] % stride[minor_dim] == 0);
    }
#else
    AXOM_UNUSED_VAR(shape);
#endif
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
    const IndexType baseIdx = offset(idx);
    assert(inBounds(baseIdx));
    return asDerived().data()[baseIdx];
  }

  /// \overload
  AXOM_HOST_DEVICE ConstSliceType<DIM> sliceImpl(
    const StackArray<IndexType, DIM>& idx) const
  {
    const IndexType baseIdx = offset(idx);
    assert(inBounds(baseIdx));
    return asDerived().data()[baseIdx];
  }
  /// @}

protected:
  /// \brief The extent in each direction
  StackArray<IndexType, DIM> m_shape;
  /// \brief For converting between multidim indices and offset.
  MDMapping<DIM> m_mapping;
  /*! \brief Cached value for optimization.  @see minStride()

    For some reason, computing min stride in minStride() slows down
    flatIndex() for CUDA and HIP, even though it doesn't seem tricky
    to optimize.  As a work around, we cache the value in m_minStrides
    and update it when m_mapping changes.  BTNG, March 2024.
  */
  IndexType m_minStride;
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

  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, 1>&, int stride = 1)
    : m_stride(stride)
  { }

  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, 1>&,
                             const StackArray<IndexType, 1>& stride)
    : m_stride(stride[0])
  { }

  AXOM_HOST_DEVICE ArrayBase(const StackArray<IndexType, 1>&,
                             const MDMapping<1>& mapping)
    : m_stride(mapping.strides()[0])
  { }

  // Empty implementation because no member data
  template <typename OtherArrayType>
  ArrayBase(const ArrayBase<typename std::remove_const<T>::type, 1, OtherArrayType>&)
  { }

  // Empty implementation because no member data
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

  /// \brief Returns the multidimensional mapping for the Array
  AXOM_HOST_DEVICE MDMapping<1> mapping() const
  {
    return MDMapping<1> {{m_stride}};
  }

  /*!
   * \brief Returns the stride between adjacent items.
   */
  AXOM_HOST_DEVICE IndexType minStride() const { return m_stride; }

  /*!
   * \brief Accessor, returns a reference to the given value.
   * For multidimensional arrays, indexes into the (flat) raw data.
   *
   * \param [in] idx the position of the value to return.
   *
   * \note equivalent to *(array.data() + idx * array.spacing()).
   *
   * \pre 0 <= idx < asDerived().size()
   */
  /// @{
  AXOM_HOST_DEVICE T& operator[](const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * m_stride];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& operator[](const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * m_stride];
  }

  /*!
   * \brief Accessor, returns a reference to the given value.
   * For multidimensional arrays, indexes into the (flat) raw data.
   *
   * \param [in] idx the position of the value to return.
   *
   * \note equivalent to *(array.data() + idx * array.minStride()).
   *
   * \pre 0 <= idx < asDerived().size()
   */
  AXOM_HOST_DEVICE T& flatIndex(const IndexType idx)
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * m_stride];
  }
  /// \overload
  AXOM_HOST_DEVICE RealConstT& flatIndex(const IndexType idx) const
  {
    assert(inBounds(idx));
    return asDerived().data()[idx * m_stride];
  }
  /// @}

  /// \brief Swaps two ArrayBases
  /// No member data, so this is a no-op
  void swap(ArrayBase&) { }

  /// \brief Set the shape
  /// No member data, so this is a no-op
  AXOM_HOST_DEVICE void setShape(const StackArray<IndexType, 1>&) { }

protected:
  /*!
   * \brief Returns the minimum "chunk size" that should be allocated
   */
  IndexType blockSize() const { return m_stride; }

  /*!
   * \brief Updates the internal dimensions and striding based on the insertion
   *  of a range of elements.
   *  No-op, since we don't keep any shape information in this specialization.
   */
  void updateShapeOnInsert(const StackArray<IndexType, 1>&) { }

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

  /*!
   * \brief Returns the size of the range of memory in which elements are
   *  located. This is equivalent to size() * minStride().
   */
  AXOM_HOST_DEVICE IndexType memorySize() const
  {
    return m_stride * shape()[0];
  }

  /*! \brief Test if idx is within bounds */
  AXOM_HOST_DEVICE inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < memorySize();
  }
  /// @}

  int m_stride {1};
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

enum class OperationSpace
{
  Host,
  Device,
  Unified_Device
};

template <typename T, OperationSpace Space>
struct ArrayOpsBase;

template <typename T>
struct ArrayOpsBase<T, OperationSpace::Host>
{
  using DefaultCtorTag = std::is_default_constructible<T>;

  /*!
   * \brief Helper for default-initializing the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
   * \note Specialization for when T is default-constructible.
   */
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
  static void init_impl(T*, IndexType, IndexType, std::false_type) { }

  /*!
   * \brief Default-initializes the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
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
   * \param [in] space the memory space in which values resides
   */
  static void fill_range(T* array,
                         IndexType begin,
                         IndexType nelems,
                         const T* values,
                         MemorySpace space)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(std::is_trivially_copyable<T>::value)
    {
      axom::copy(array + begin, values, sizeof(T) * nelems);
    }
    else
    {
      void* values_buf = nullptr;
      const T* values_host = values;
      if(space == MemorySpace::Device)
      {
        values_buf = ::operator new(sizeof(T) * nelems);
        // "Relocate" the device-side values into host memory, before copying
        // into uninitialized memory
        axom::copy(values_buf, values, sizeof(T) * nelems);
        values_host = static_cast<T*>(values_buf);
      }
      std::uninitialized_copy(values_host, values_host + nelems, array + begin);
      if(values_buf)
      {
        ::operator delete(values_buf);
      }
    }
#else
    AXOM_UNUSED_VAR(space);
    std::uninitialized_copy(values, values + nelems, array + begin);
#endif
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
      auto rbegin = std::reverse_iterator<T*>(array + src_end);
      auto rend = std::reverse_iterator<T*>(array + src_begin);
      auto rdest = std::reverse_iterator<T*>(array + dst_last);
      // Do an "uninitialized-move" in reverse order, to avoid overwriting
      // any existing elements.
      std::uninitialized_copy(std::make_move_iterator(rbegin),
                              std::make_move_iterator(rend),
                              rdest);
    }
    else if(src_begin > dst)
    {
      // This substitutes for std::uninitialized_move(), which is only
      // available in C++17.
      std::uninitialized_copy(std::make_move_iterator(array + src_begin),
                              std::make_move_iterator(array + src_end),
                              array + dst);
    }
  }

  /*!
   * \brief Moves a range of elements to a new allocation.
   *
   * \param [inout] array the array to move the elements to.
   * \param [in] nelems the number of elements to move.
   * \param [in] values the destination index of the range of elements
   */
  static void realloc_move(T* array, IndexType nelems, T* values)
  {
    std::uninitialized_copy(std::make_move_iterator(values),
                            std::make_move_iterator(values + nelems),
                            array);
    destroy(values, 0, nelems);
  }
};

#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
/*!
 * \name Tag types for device initialization
 */
/// @{

/*!
 * \brief Tag type representing that a type can be initialized on the device.
 *
 *  This only applies to types which are trivially device-constructible.
 */
struct InitTypeOnDevice
{ };
/*!
 * \brief Tag type representing that a type can be initialized on the device.
 *
 *  This applies to types which are not trivially default-constructible, but are
 *  trivially-copyable; we can construct a default value on the host and copy-
 *  construct values with it on the device.
 */
struct InitTypeOnDeviceWithCopy
{ };
/*!
 * \brief Tag type representing that a type cannot be initialized on the device.
 */
struct InitTypeOnHost
{ };

/*!
 * \brief Selector type which matches a type to its corresponding initialization
 *  tag type.
 */
template <typename T, typename Enable = void>
struct DeviceInitTag
{
  using Type = InitTypeOnHost;
};

template <typename T>
struct DeviceInitTag<T, std::enable_if_t<std::is_trivially_default_constructible<T>::value>>
{
  using Type = InitTypeOnDevice;
};

template <typename T>
struct DeviceInitTag<T,
                     std::enable_if_t<!std::is_trivially_default_constructible<T>::value &&
                                      std::is_default_constructible<T>::value &&
                                      std::is_trivially_copyable<T>::value>>
{
  using Type = InitTypeOnDeviceWithCopy;
};
/// @}

template <typename T, OperationSpace SPACE>
struct DeviceStagingBuffer;

template <typename T>
struct DeviceStagingBuffer<T, OperationSpace::Device>
{
  /*!
   * \brief Create a staging buffer for a device memory operation.
   *
   * \param [in] data the array to mirror on the CPU
   * \param [in] begin the beginning index of the range to mirror
   * \param [in] nelems the number of elements to mirror
   * \param [in] read_from_data if true, copies existing data range to the
   *  staging buffer during construction
   */
  DeviceStagingBuffer(T* data,
                      IndexType begin,
                      IndexType nelems,
                      bool read_from_data = false)
    : m_data(data)
    , m_begin(begin)
    , m_num_elems(nelems)
  {
    int allocator_id = 0;
  #ifdef AXOM_USE_UMPIRE
    allocator_id = axom::detail::getAllocatorID<axom::MemorySpace::Host>();
  #endif
    m_staging_buf = axom::allocate<T>(nelems, allocator_id);
    if(read_from_data)
    {
      axom::copy(m_staging_buf, m_data + begin, sizeof(T) * nelems);
    }
  }

  DISABLE_COPY_AND_ASSIGNMENT(DeviceStagingBuffer);
  DISABLE_MOVE_AND_ASSIGNMENT(DeviceStagingBuffer);

  ~DeviceStagingBuffer()
  {
    // Copy back staging data to destination buffer.
    axom::copy(m_data + m_begin, m_staging_buf, m_num_elems * sizeof(T));
    axom::deallocate(m_staging_buf);
  }

  T* getStagingBuffer() const { return static_cast<T*>(m_staging_buf); }

  T* m_staging_buf;
  T* m_data;
  IndexType m_begin;
  IndexType m_num_elems;
};

template <typename T>
struct DeviceStagingBuffer<T, OperationSpace::Unified_Device>
{
  DeviceStagingBuffer(T* data,
                      IndexType begin,
                      IndexType nelems,
                      bool read_from_data = false)
    : m_data(data)
    , m_begin(begin)
  {
    AXOM_UNUSED_VAR(nelems);
    AXOM_UNUSED_VAR(read_from_data);
  }

  T* getStagingBuffer() const { return static_cast<T*>(m_data + m_begin); }

  T* m_data;
  IndexType m_begin;
};

template <typename T, OperationSpace SPACE>
struct ArrayOpsBase
{
  #if defined(__CUDACC__)
  using ExecSpace = axom::CUDA_EXEC<256>;
  #else
  using ExecSpace = axom::HIP_EXEC<256>;
  #endif

  static constexpr bool DestroyOnHost = !std::is_trivially_destructible<T>::value;
  static constexpr bool DefaultCtor = std::is_default_constructible<T>::value;

  using HostOp = ArrayOpsBase<T, OperationSpace::Host>;
  using StagingBuffer = DeviceStagingBuffer<T, SPACE>;

  /*!
   * \brief Helper for default-initialization of a range of elements.
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
   * \note Specialization for when T is only initializable on the host.
   */
  static void init_impl(T* data, IndexType begin, IndexType nelems, InitTypeOnHost)
  {
    if(std::is_default_constructible<T>::value)
    {
      // If we instantiated a fill kernel here it would require
      // that T's default ctor is device-annotated which is too
      // strict of a requirement, so we copy a buffer instead.
      StagingBuffer tmp_buf(data, begin, nelems);
      HostOp::init(tmp_buf.getStagingBuffer(), 0, nelems);
    }
  }

  /*!
   * \overload
   * \note Specialization for when T is trivially default-constructible.
   */
  static void init_impl(T* data, IndexType begin, IndexType nelems, InitTypeOnDevice)
  {
    for_all<ExecSpace>(
      begin,
      begin + nelems,
      AXOM_LAMBDA(IndexType i) { new(&data[i]) T(); });
  }

  /*!
   * \overload
   * \note Specialization for when T is trivially copyable.
   */
  static void init_impl(T* data,
                        IndexType begin,
                        IndexType nelems,
                        InitTypeOnDeviceWithCopy)
  {
    T object {};
    for_all<ExecSpace>(
      begin,
      begin + nelems,
      AXOM_LAMBDA(IndexType i) { new(&data[i]) T(object); });
  }

  /*!
   * \brief Default-initializes the "new" segment of an array
   *
   * \param [inout] data The data to initialize
   * \param [in] begin The beginning of the subset of \a data that should be initialized
   * \param [in] nelems the number of elements to initialize
   */
  static void init(T* data, IndexType begin, IndexType nelems)
  {
    init_impl(data, begin, nelems, typename DeviceInitTag<T>::Type {});
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
  static void fill_impl(T* array,
                        IndexType begin,
                        IndexType nelems,
                        const T& value,
                        std::false_type)
  {
    // If we instantiated a fill kernel here it would require
    // that T's copy ctor is device-annotated which is too
    // strict of a requirement, so we copy a buffer instead.
    StagingBuffer tmp_buf(array, begin, nelems);
    HostOp::fill(tmp_buf.getStagingBuffer(), 0, nelems, value);
  }

  /*!
   * \overload
   * \note Specialization for when T is trivially-copyable.
   */
  static void fill_impl(T* array,
                        IndexType begin,
                        IndexType nelems,
                        const T& value,
                        std::true_type)
  {
    for_all<ExecSpace>(
      nelems,
      AXOM_LAMBDA(IndexType i) { new(&array[i + begin]) T(value); });
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
    fill_impl(array, begin, nelems, value, std::is_trivially_copyable<T> {});
  }

  /*!
   * \brief Fills an uninitialized array with a range of objects of type T.
   *
   * \param [inout] array the array to fill
   * \param [in] begin the index at which to begin placing elements
   * \param [in] nelems the number of elements in the range to fill the array with
   * \param [in] values the values to set each array element to
   * \param [in] space the memory space in which values resides
   */
  static void fill_range(T* array,
                         IndexType begin,
                         IndexType nelems,
                         const T* values,
                         MemorySpace space)
  {
    if(std::is_trivially_copyable<T>::value)
    {
      axom::copy(array + begin, values, sizeof(T) * nelems);
    }
    else
    {
      // HostOp::fill_range will handle the copy to our "staging" host buffer,
      // regardless of the source memory space.
      StagingBuffer tmp_buf(array, begin, nelems);
      HostOp::fill_range(tmp_buf.getStagingBuffer(), 0, nelems, values, space);
    }
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
    if(SPACE == OperationSpace::Device)
    {
      // Similar to fill(), except we can allocate stack memory and placement-new
      // the object with a move constructor.
      alignas(T) std::uint8_t host_buf[sizeof(T)];
      T* host_obj = ::new(&host_buf) T(std::forward<Args>(args)...);
      axom::copy(array + i, host_obj, sizeof(T));
    }
    else  // SPACE == OperationSpace::Unified_Device
    {
      // Construct directly in unified/pinned memory.
      ::new(array + i) T(std::forward<Args>(args)...);
    }
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
      StagingBuffer tmp_buf(array, begin, nelems, true);
      HostOp::destroy(tmp_buf.getStagingBuffer(), 0, nelems);
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
    if(!std::is_trivially_copyable<T>::value &&
       SPACE == OperationSpace::Unified_Device)
    {
      // Type might not be trivially-relocatable, move the range on the host.
      // Note that we only do this for objects in unified/pinned memory, since
      // we assume that objects in device-only memory are trivially-relocatable.
      HostOp::move(array, src_begin, src_end, dst);
    }
    else
    {
      // Since this memory is on the device-side, we copy it to a temporary buffer
      // first.
      IndexType nelems = src_end - src_begin;
      T* tmp_buf =
        axom::allocate<T>(nelems,
                          axom::execution_space<ExecSpace>::allocatorID());
      axom::copy(tmp_buf, array + src_begin, nelems * sizeof(T));
      axom::copy(array + dst, tmp_buf, nelems * sizeof(T));
      axom::deallocate(tmp_buf);
    }
  }

  /*!
   * \brief Moves a range of elements to a new allocation.
   *
   * \param [inout] array the array to move the elements to.
   * \param [in] nelems the number of elements to move.
   * \param [in] values the destination index of the range of elements
   */
  static void realloc_move(T* array, IndexType nelems, T* values)
  {
    if(!std::is_trivially_copyable<T>::value &&
       SPACE == OperationSpace::Unified_Device)
    {
      HostOp::realloc_move(array, nelems, values);
    }
    else
    {
      // NOTE: technically this is incorrect for non-trivially relocatable types,
      // but since we only support trivially-relocatable types in device-only
      // memory, a bitcopy will suffice.
      axom::copy(array, values, nelems * sizeof(T));
    }
  }
};
#endif

template <MemorySpace SPACE>
struct MemSpaceTraits
{
  static constexpr OperationSpace Space = OperationSpace::Host;
  // True if memory is accessible by both the host and device. False otherwise.
  static constexpr bool IsUVMAccessible = false;
};

#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
template <>
struct MemSpaceTraits<MemorySpace::Device>
{
  #if defined(AXOM_USE_CUDA)
  // On CUDA platforms, device memory allocated with cudaMalloc can only be
  // touched from a device kernel.
  static constexpr OperationSpace Space = OperationSpace::Device;
  static constexpr bool IsUVMAccessible = false;
  #elif defined(AXOM_USE_HIP)
  // On HIP platforms, device memory allocated with hipMalloc is accessible from
  // the host.
  static constexpr OperationSpace Space = OperationSpace::Unified_Device;
  static constexpr bool IsUVMAccessible = true;
  #endif
};

template <>
struct MemSpaceTraits<MemorySpace::Pinned>
{
  static constexpr OperationSpace Space = OperationSpace::Unified_Device;
  static constexpr bool IsUVMAccessible = true;
};
template <>
struct MemSpaceTraits<MemorySpace::Unified>
{
  static constexpr OperationSpace Space = OperationSpace::Unified_Device;
  static constexpr bool IsUVMAccessible = true;
};

template <>
struct MemSpaceTraits<MemorySpace::Dynamic>
{
  static constexpr bool IsUVMAccessible = true;
};
#endif

template <typename T,
          MemorySpace SPACE,
          bool IsUVMAccessible = MemSpaceTraits<SPACE>::IsUVMAccessible>
struct ArrayOps;

template <typename T, MemorySpace SPACE>
struct ArrayOps<T, SPACE, false>
{
private:
  constexpr static OperationSpace OpSpace = MemSpaceTraits<SPACE>::Space;

  using Base = ArrayOpsBase<T, OpSpace>;

public:
  ArrayOps(int allocId, bool preferDevice)
  {
    AXOM_UNUSED_VAR(allocId);
    AXOM_UNUSED_VAR(preferDevice);
  }

  void init(T* array, IndexType begin, IndexType nelems)
  {
    Base::init(array, begin, nelems);
  }

  void fill(T* array, IndexType begin, IndexType nelems, const T& value)
  {
    Base::fill(array, begin, nelems, value);
  }

  void fill_range(T* array,
                  IndexType begin,
                  IndexType nelems,
                  const T* values,
                  MemorySpace space)
  {
    Base::fill_range(array, begin, nelems, values, space);
  }

  void destroy(T* array, IndexType begin, IndexType nelems)
  {
    if(nelems == 0)
    {
      return;
    }
    Base::destroy(array, begin, nelems);
  }

  void move(T* array, IndexType src_begin, IndexType src_end, IndexType dst)
  {
    if(src_begin >= src_end)
    {
      return;
    }
    Base::move(array, src_begin, src_end, dst);
  }

  void realloc_move(T* array, IndexType nelems, T* values)
  {
    Base::realloc_move(array, nelems, values);
  }

  template <typename... Args>
  void emplace(T* array, IndexType dst, Args&&... args)
  {
    Base::emplace(array, dst, std::forward<Args>(args)...);
  }
};

template <typename T, MemorySpace SPACE>
struct ArrayOps<T, SPACE, true>
{
private:
  using Base = ArrayOpsBase<T, OperationSpace::Host>;
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
  using BaseDevice = ArrayOpsBase<T, OperationSpace::Device>;
  // Works with unified and pinned memory.
  using BaseUM = ArrayOpsBase<T, OperationSpace::Unified_Device>;

  MemorySpace space {MemorySpace::Dynamic};
#endif

public:
  ArrayOps(int allocId, bool preferDevice)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(SPACE == MemorySpace::Dynamic)
    {
      space = getAllocatorSpace(allocId);
    }
    else
    {
      space = SPACE;
    }

    bool isUnifiedSpace = false;
    isUnifiedSpace =
      (space == MemorySpace::Unified || space == MemorySpace::Pinned);
  #if defined(AXOM_USE_HIP)
    isUnifiedSpace = (isUnifiedSpace || space == MemorySpace::Device);
  #endif
    if(!preferDevice && isUnifiedSpace)
    {
      space = MemorySpace::Host;
    }
#else
    AXOM_UNUSED_VAR(allocId);
    AXOM_UNUSED_VAR(preferDevice);
#endif
  }

  void init(T* array, IndexType begin, IndexType nelems)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::init(array, begin, nelems);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::init(array, begin, nelems);
      return;
    }
#endif
    Base::init(array, begin, nelems);
  }

  void fill(T* array, IndexType begin, IndexType nelems, const T& value)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::fill(array, begin, nelems, value);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::fill(array, begin, nelems, value);
      return;
    }
#endif
    Base::fill(array, begin, nelems, value);
  }

  void fill_range(T* array,
                  IndexType begin,
                  IndexType nelems,
                  const T* values,
                  MemorySpace valueSpace)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::fill_range(array, begin, nelems, values, valueSpace);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::fill_range(array, begin, nelems, values, valueSpace);
      return;
    }
#endif
    Base::fill_range(array, begin, nelems, values, valueSpace);
  }

  void destroy(T* array, IndexType begin, IndexType nelems)
  {
    if(nelems == 0)
    {
      return;
    }
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::destroy(array, begin, nelems);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::destroy(array, begin, nelems);
      return;
    }
#endif
    Base::destroy(array, begin, nelems);
  }

  void move(T* array, IndexType src_begin, IndexType src_end, IndexType dst)
  {
    if(src_begin >= src_end)
    {
      return;
    }
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::move(array, src_begin, src_end, dst);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::move(array, src_begin, src_end, dst);
      return;
    }
#endif
    Base::move(array, src_begin, src_end, dst);
  }

  void realloc_move(T* array, IndexType nelems, T* values)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::realloc_move(array, nelems, values);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::realloc_move(array, nelems, values);
      return;
    }
#endif
    Base::realloc_move(array, nelems, values);
  }

  template <typename... Args>
  void emplace(T* array, IndexType dst, Args&&... args)
  {
#if defined(AXOM_USE_GPU) && defined(AXOM_GPUCC) && defined(AXOM_USE_UMPIRE)
    if(space == MemorySpace::Device)
    {
      BaseDevice::emplace(array, dst, std::forward<Args>(args)...);
      return;
    }
    else if(space == MemorySpace::Unified || space == MemorySpace::Pinned)
    {
      BaseUM::emplace(array, dst, std::forward<Args>(args)...);
      return;
    }
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
    : BaseClass(detail::takeLastElems<SliceDim>(array->shape()),
                detail::takeLastElems<SliceDim>(array->strides()))
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

  IndexType spacing() const { return m_array->spacing(); }

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
