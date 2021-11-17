// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAY_HPP_
#define AXOM_ARRAY_HPP_

#include "axom/config.hpp"                    // for compile-time defines
#include "axom/core/Macros.hpp"               // for axom macros
#include "axom/core/utilities/Utilities.hpp"  // for processAbort()
#include "axom/core/Types.hpp"                // for IndexType definition
#include "axom/core/ArrayBase.hpp"
#include "axom/core/ArrayIteratorBase.hpp"

// C/C++ includes
#include <algorithm>  // for std::transform
#include <iostream>   // for std::cerr and std::ostream

namespace axom
{
// Forward declare the templated classes and operator function(s)
template <typename T, int DIM, MemorySpace SPACE>
class Array;

/*!
 * \class Array
 *
 * \brief Provides a generic multidimensional array container.
 *
 *  The Array class provides a generic multidimensional array 
 *  container with dynamic reallocation and insertion.  The dimensionality
 *  of the array must be known at compile time but the extents in each dimension
 *  are dynamic and can be changed at runtime.  Array elements are stored
 *  contiguously.
 *
 *  \note For a multi-component array container, where each element
 *  is a tuple of 1 or more components, Axom provides the MCArray alias, which
 *  corresponds to Array<T, 2>.
 *
 *  The Array class mirrors std::vector, with future support for GPUs
 *  in-development.  The class's multidimensional array functionality roughly
    mirrors the multidimensional array support provided by numpy's ndarray.
 * 
 *  \see https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html
 * 
 *  This class is meant to be a drop-in replacement for std::vector.
 *  However, it differs in its memory management and construction semantics.
 *  Specifically, we do not require axom::Array to initialize/construct
 *  its memory at allocation time and we use axom's memory_management
 *  and allocator ID abstractions rather than std::allocator.
 *
 *  Array always retains exclusive ownership of its data and is responsible for
 *  freeing its memory.
 *
 *  \see ArrayView for non-owning views of one- or multi-dimensional data
 *  Depending on which constructor is used, the Array object can have two
 *  different underlying storage types:
 *
 * \tparam T the type of the values to hold.
 * \tparam DIM The dimension of the array.
 * 
 * \pre T must be CopyAssignable and Erasable
 * \see https://en.cppreference.com/w/cpp/named_req
 *
 */
template <typename T, int DIM = 1, MemorySpace SPACE = MemorySpace::Dynamic>
class Array : public ArrayBase<T, DIM, Array<T, DIM, SPACE>>
{
public:
  static constexpr double DEFAULT_RESIZE_RATIO = 2.0;
  static constexpr IndexType MIN_DEFAULT_CAPACITY = 32;
  using value_type = T;
  static constexpr MemorySpace space = SPACE;
  using ArrayIterator = ArrayIteratorBase<Array<T, DIM, SPACE>>;

public:
  /// \name Native Storage Array Constructors
  /// @{

  /*! 
   * \brief Default constructor. Constructs an Array instance with no elements
   *  and default allocator ID. 
   *
   */
  Array();

  /*!
   * \brief Constructs a 1D Array instance with the given number of elements.
   *
   * \param [in] num_elements the number of elements the Array holds.
   * \param [in] capacity the number of elements to allocate space for.
   * \param [in] allocator_id the ID of the allocator to use (optional)
   *
   * \note If no capacity or capacity less than num_elements is specified
   *  then it will default to at least num_elements * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of elements to store in the
   *  array and does not correspond to the actual bytesize.
   * \note The option to select a capacity is only available for 1-dimensional Arrays
   *
   * \pre num_elements >= 0
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <IndexType SFINAE_DIM = DIM,
            MemorySpace SFINAE_SPACE = SPACE,
            typename std::enable_if<SFINAE_DIM == 1 &&
                                    SFINAE_SPACE == MemorySpace::Dynamic>::type* = nullptr>
  Array(IndexType num_elements,
        IndexType capacity = 0,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /// \overload
  template <IndexType SFINAE_DIM = DIM,
            MemorySpace SFINAE_SPACE = SPACE,
            typename std::enable_if<SFINAE_DIM == 1 &&
                                    SFINAE_SPACE != MemorySpace::Dynamic>::type* = nullptr>
  Array(IndexType num_elements, IndexType capacity = 0);

  /*!
   * \brief Generic constructor for an Array of arbitrary dimension
   *
   * \param [in] args The parameter pack containing the "shape" of the Array
   * \see https://numpy.org/doc/stable/reference/generated/numpy.empty.html#numpy.empty
   *
   * \pre sizeof...(Args) == DIM
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <typename... Args,
            typename std::enable_if<
              detail::all_types_are_integral<Args...>::value>::type* = nullptr>
  Array(Args... args);

  /*! 
   * \brief Copy constructor for an Array instance 
   * 
   * \param [in] allocator_id the ID of the allocator to use (optional)
   */
  Array(const Array& other,
        int allocator_id = axom::detail::getAllocatorID<SPACE>());

  /*! 
   * \brief Move constructor for an Array instance 
   */
  Array(Array&& other);

  /*! 
   * \brief Constructor for transferring between memory spaces
   * 
   * \param [in] other The array in a different memory space to copy from
   */
  template <typename OtherArrayType>
  Array(const ArrayBase<T, DIM, OtherArrayType>& other);

  /// @}

  /// \name Array copy and move operators
  /// @{

  /*! 
   * \brief Copy assignment operator for Array 
   *
   * \note The data will be allocated using the allocator ID of the
   *  copy-assigned Array, not the argument Array.
   * 
   * \pre T must be TriviallyCopyable
   */
  Array& operator=(const Array& other)
  {
    if(this != &other)
    {
      static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE>>&>(*this) = other;
      m_resize_ratio = other.m_resize_ratio;
      initialize(other.size(), other.capacity());
      axom::copy(m_data, other.data(), m_num_elements * sizeof(T));
    }

    return *this;
  }

  /*! 
   * \brief Move assignment operator for Array
   */
  Array& operator=(Array&& other)
  {
    if(this != &other)
    {
      if(m_data != nullptr)
      {
        axom::deallocate(m_data);
      }
      static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE>>&>(*this) =
        std::move(other);

      m_data = other.m_data;
      m_num_elements = other.m_num_elements;
      m_capacity = other.m_capacity;
      m_resize_ratio = other.m_resize_ratio;
      m_allocator_id = other.m_allocator_id;

      other.m_data = nullptr;
      other.m_num_elements = 0;
      other.m_capacity = 0;
      other.m_resize_ratio = DEFAULT_RESIZE_RATIO;
      other.m_allocator_id = INVALID_ALLOCATOR_ID;
    }

    return *this;
  }

  /// @}

  /*!
   * Destructor. Frees the associated buffer.
   */
  virtual ~Array();

  /// \name Array element access operators
  /// @{

  // TODO: Implement View class for the case where sizeof...(Args) < DIM (i.e., where the indexing results in a nonscalar)

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  AXOM_HOST_DEVICE inline T* data() { return m_data; }
  AXOM_HOST_DEVICE inline const T* data() const { return m_data; }

  /// @}

  /// @}

  /// \name Array methods to modify the data.
  /// @{

  /*!
   * \brief Set all the values of the array.
   *
   * \param [in] value the value to set to.
   */
  void fill(const T& value);

  /*!
   * \brief Modify the values of existing elements.
   *
   * \param [in] elements the new elements to write.
   * \param [in] n the number of elements to write.
   * \param [in] pos the position at which to begin writing.
   *
   * \note It's assumed that elements is of length n.
   * \note The size is unchanged by calls to set.
   *
   * \pre pos + n <= m_num_elements.
   */
  void set(const T* elements, IndexType n, IndexType pos);

  /*!
   * \brief Clears the contents of the array
   * 
   * \post size of Array is 0
   * \post capacity is unchanged
   */
  void clear();

  /*!
   * \brief Insert an element into the array at the given position.
   *
   * \param [in] pos the position at which to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   */
  void insert(IndexType pos, const T& value);

  /*!
   * \brief Insert an element into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \return ArrayIterator to inserted value
   */
  ArrayIterator insert(ArrayIterator pos, const T& value);

  /*!
   * \brief Insert elements into the array at the given position.
   *
   * \param [in] pos the position at which to begin the insertion.
   * \param [in] n the number of elements to insert.
   * \param [in] values the element values to insert.
   *
   * \note It's assumed that elements is of length n.
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   *
   * \pre pos <= m_num_elements.
   */
  void insert(IndexType pos, IndexType n, const T* values);

  /*!
   * \brief Insert elements into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] n the number of elements to insert.
   * \param [in] values the element values to insert.
   *
   * \note It's assumed that elements is of length n.
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   *
   * \pre pos <= end()
   *
   * \return ArrayIterator to first element inserted (pos if n == 0)
   */
  ArrayIterator insert(ArrayIterator pos, IndexType n, const T* values);

  /*!
   * \brief Insert n copies of element into the array at the given position.
   *
   * \param [in] pos the position at which to begin the insertion.
   * \param [in] n the number of elements to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   * \note This method is used to create space for elements in the middle of
   *  the array.
   *
   * \pre pos <= m_num_elements.
   */
  void insert(IndexType pos, IndexType n, const T& value);

  /*!
   * \brief Insert n copies of element into the array at the value before pos.
   *
   * \param [in] pos the ArrayIterator before which value will be inserted.
   * \param [in] n the number of elements to insert.
   * \param [in] value the element value to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   * \note This method is used to create space for elements in the middle of
   *  the array.
   *
   * \pre pos <= end()
   *
   * \return ArrayIterator to first element inserted (pos if n == 0)
   */
  ArrayIterator insert(ArrayIterator pos, IndexType n, const T& value);

  // Make the overload "visible"
  using ArrayBase<T, DIM, Array<T, DIM, SPACE>>::insert;

  /*!
   * \brief Appends an Array to the end of the calling object
   *
   * \param [in] other The Array to append
   * \tparam OtherArrayType The underlying type of the other array
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <typename OtherArrayType>
  void append(const ArrayBase<T, DIM, OtherArrayType>& other)
  {
    ArrayBase<T, DIM, Array<T, DIM, SPACE>>::insert(size(), other);
  }

  /*!
   * \brief Erases an element from the Array 
   *
   * \param [in] pos the ArrayIterator to the element in the Array
   *
   * \return An ArrayIterator following the last element removed.
   */
  ArrayIterator erase(ArrayIterator pos);

  /*!
   * \brief Erases elements in the range [first, last) from the Array
   *
   * \param [in] first the ArrayIterator to the beginning of the range.
   * \param [in] last the ArrayIterator to end of range.
   *
   * \return An ArrayIterator following the last element removed. 
   */
  ArrayIterator erase(ArrayIterator first, ArrayIterator last);

  /*!
   * \brief Inserts new element into Array at the given position.
   *
   * \param [in] pos the position to insert element at.
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre T must be MoveAssignable
   */
  template <typename... Args>
  void emplace(IndexType pos, Args&&... args);

  /*!
   * \brief Inserts new element into Array before pos.
   *
   * \param [in] pos the ArrayIterator to insert element before.
   * \param [in] args the arguments to forward to constructor of the element.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre T must be MoveAssignable
   *
   * \return An ArrayIterator to the emplaced element.
   */
  template <typename... Args>
  ArrayIterator emplace(ArrayIterator pos, Args&&... args);

  /// @}

  /// \name Array methods to query and set attributes
  /// @{

  /*!
   * \brief Return the number of elements allocated for the data array.
   */
  IndexType capacity() const { return m_capacity; }

  /*!
   * \brief Increase the capacity. Does nothing if the new capacity is less
   *  than the current capacity.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  void reserve(IndexType capacity)
  {
    if(capacity > m_capacity)
    {
      setCapacity(capacity);
    }
  }

  /*!
   * \brief Returns an ArrayIterator to the first element of the Array
   */
  ArrayIterator begin()
  {
    assert(m_data != nullptr);
    return ArrayIterator(0, this);
  }

  /*!
   * \brief Returns an ArrayIterator to the element following the last
   *  element of the Array.
   */
  ArrayIterator end()
  {
    assert(m_data != nullptr);
    return ArrayIterator(size(), this);
  }

  /*!
   * \brief Shrink the capacity to be equal to the size.
   */
  void shrink() { setCapacity(m_num_elements); }

  /*!
   * \brief Returns true iff the Array stores no elements.
   *
   * \note If the Array is empty the capacity can still be greater than zero.
   */
  bool empty() const { return m_num_elements == 0; }

  /*!
   * \brief Return the number of elements stored in the data array.
   */
  AXOM_HOST_DEVICE inline IndexType size() const { return m_num_elements; }

  /*!
   * \brief Update the number of elements stored in the data array.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  template <typename... Args>
  void resize(Args... args);

  /*!
   * \brief Exchanges the contents of this Array with the other.
   */
  void swap(Array<T, DIM, SPACE>& other);

  /*!
   * \brief Get the ratio by which the capacity increases upon dynamic resize.
   */
  double getResizeRatio() const { return m_resize_ratio; }

  /*!
   * \brief Set the ratio by which the capacity increases upon dynamic resize.
   *
   * \param [in] ratio the new resize ratio.
   */
  void setResizeRatio(double ratio) { m_resize_ratio = ratio; }

  /*!
   * \brief Get the ID for the umpire allocator
   */
  int getAllocatorID() const { return m_allocator_id; }

  /// @}

protected:
  /*!
   * \brief Initialize an Array instance with the given number of elements.
   *
   * \param [in] num_elements the number of elements the Array holds.
   * \param [in] capacity the number of elements to allocate space for.
   *
   * \note If no capacity or capacity less than num_elements is specified
   *  then it will default to at least num_elements * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of elements to store in the
   *  array and does not correspond to the actual bytesize.
   *
   * \pre num_elements >= 0
   *
   * \post capacity() >= size()
   * \post size() == num_elements
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  void initialize(IndexType num_elements, IndexType capacity);

  /*!
   * \brief Make space for a subsequent insertion into the array.
   *
   * \param [in] n the number of elements to insert.
   * \param [in] pos the position at which to begin the insertion.
   *
   * \return a pointer to the beginning of the insertion space.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  T* reserveForInsert(IndexType n, IndexType pos);

  /*!
   * \brief Update the number of elements.
   *
   * \param [in] new_num_elements the new number of elements.
   */
  virtual void updateNumElements(IndexType new_num_elements);

  /*!
   * \brief Set the number of elements allocated for the data array.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  virtual void setCapacity(IndexType new_capacity);

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_elements the number of elements which exceeds the
   *  current capacity.
   */
  virtual void dynamicRealloc(IndexType new_num_elements);

  T* m_data = nullptr;
  /// \brief The full number of elements in the array
  ///  i.e., 3 for a 1D Array of size 3, 9 for a 3x3 2D array, etc
  IndexType m_num_elements = 0;
  IndexType m_capacity = 0;
  double m_resize_ratio = DEFAULT_RESIZE_RATIO;
  int m_allocator_id;
};

/// \brief Helper alias for multi-component arrays
template <typename T>
using MCArray = Array<T, 2>;

//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
Array<T, DIM, SPACE>::Array()
  : m_allocator_id(axom::detail::getAllocatorID<SPACE>())
{ }

template <typename T, int DIM, MemorySpace SPACE>
template <typename... Args,
          typename std::enable_if<detail::all_types_are_integral<Args...>::value>::type*>
Array<T, DIM, SPACE>::Array(Args... args)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE>>(args...)
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
{
  static_assert(sizeof...(Args) == DIM,
                "Array size must match number of dimensions");
  // Intel hits internal compiler error when casting as part of function call
  const IndexType tmp_args[] = {args...};
  assert(detail::allNonNegative(tmp_args));
  initialize(detail::packProduct(tmp_args), 0);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <IndexType SFINAE_DIM,
          MemorySpace SFINAE_SPACE,
          typename std::enable_if<SFINAE_DIM == 1 &&
                                  SFINAE_SPACE == MemorySpace::Dynamic>::type*>
Array<T, DIM, SPACE>::Array(IndexType num_elements,
                            IndexType capacity,
                            int allocator_id)
  : m_allocator_id(allocator_id)
{
  initialize(num_elements, capacity);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <IndexType SFINAE_DIM,
          MemorySpace SFINAE_SPACE,
          typename std::enable_if<SFINAE_DIM == 1 &&
                                  SFINAE_SPACE != MemorySpace::Dynamic>::type*>
Array<T, DIM, SPACE>::Array(IndexType num_elements, IndexType capacity)
  : m_allocator_id(axom::detail::getAllocatorID<SPACE>())
{
  initialize(num_elements, capacity);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
Array<T, DIM, SPACE>::Array(const Array& other, int allocator_id)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE>>(
      static_cast<const ArrayBase<T, DIM, Array<T, DIM, SPACE>>&>(other))
  , m_allocator_id(SPACE == MemorySpace::Dynamic
                     ? allocator_id
                     : axom::detail::getAllocatorID<SPACE>())
{
// We can't template/SFINAE away the allocator_id parameter since this is a copy
// constructor, so we just ignore the allocator ID if the memory space isn't Dynamic.
// We can warn the user that their input is being ignored, though.
#ifdef AXOM_DEBUG
  if(SPACE != MemorySpace::Dynamic &&
     allocator_id != axom::detail::getAllocatorID<SPACE>())
  {
    std::cerr << "Incorrect allocator ID was provided for an Array object with "
                 "explicit memory space\n";
  }
#endif
  initialize(other.size(), other.capacity());
  axom::copy(m_data, other.data(), m_num_elements * sizeof(T));
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
Array<T, DIM, SPACE>::Array(Array&& other)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE>>(
      static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE>>&&>(std::move(other)))
  , m_resize_ratio(0.0)
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
{
  m_data = other.m_data;
  m_num_elements = other.m_num_elements;
  m_capacity = other.m_capacity;
  m_resize_ratio = other.m_resize_ratio;
  m_allocator_id = other.m_allocator_id;

  other.m_data = nullptr;
  other.m_capacity = 0;
  other.m_resize_ratio = DEFAULT_RESIZE_RATIO;
  other.m_allocator_id = INVALID_ALLOCATOR_ID;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename OtherArrayType>
Array<T, DIM, SPACE>::Array(const ArrayBase<T, DIM, OtherArrayType>& other)
  : ArrayBase<T, DIM, Array<T, DIM, SPACE>>(other)
  , m_allocator_id(axom::detail::getAllocatorID<SPACE>())
{
  initialize(static_cast<const OtherArrayType&>(other).size(),
             static_cast<const OtherArrayType&>(other).size());
  // axom::copy is aware of pointers registered in Umpire, so this will handle
  // the transfer between memory spaces
  axom::copy(m_data,
             static_cast<const OtherArrayType&>(other).data(),
             m_num_elements * sizeof(T));
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
Array<T, DIM, SPACE>::~Array()
{
  if(m_data != nullptr)
  {
    axom::deallocate(m_data);
  }

  m_data = nullptr;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::fill(const T& value)
{
  for(IndexType i = 0; i < m_num_elements; i++)
  {
    m_data[i] = value;
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::set(const T* elements, IndexType n, IndexType pos)
{
  assert(elements != nullptr);
  assert(pos >= 0);
  assert(pos + n <= m_num_elements);

  for(IndexType i = 0; i < n; ++i)
  {
    m_data[pos + i] = elements[i];
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::clear()
{
  // This most likely needs to be a call to erase() instead.
  for(IndexType i = 0; i < m_num_elements; ++i)
  {
    m_data[i].~T();
  }

  updateNumElements(0);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::insert(IndexType pos, const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  reserveForInsert(1, pos);
  m_data[pos] = value;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::insert(
  Array<T, DIM, SPACE>::ArrayIterator pos,
  const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), value);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::insert(IndexType pos, IndexType n, const T* values)
{
  assert(values != nullptr);
  reserveForInsert(n, pos);
  for(IndexType i = 0; i < n; ++i)
  {
    m_data[pos + i] = values[i];
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::insert(
  Array<T, DIM, SPACE>::ArrayIterator pos,
  IndexType n,
  const T* values)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), n, values);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::insert(IndexType pos, IndexType n, const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  reserveForInsert(n, pos);
  for(IndexType i = 0; i < n; ++i)
  {
    m_data[pos + i] = value;
  }
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::insert(
  Array<T, DIM, SPACE>::ArrayIterator pos,
  IndexType n,
  const T& value)
{
  static_assert(DIM == 1, "Insertion not supported for multidimensional Arrays");
  assert(pos >= begin() && pos <= end());
  insert(pos - begin(), n, value);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::erase(
  Array<T, DIM, SPACE>::ArrayIterator pos)
{
  assert(pos >= begin() && pos < end());
  int counter = 0;

  while(pos < end() - 1)
  {
    *pos = *(pos + 1);
    pos += 1;
    counter += 1;
  }
  (*pos).~T();

  updateNumElements(m_num_elements - 1);
  return pos - counter;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::erase(
  Array<T, DIM, SPACE>::ArrayIterator first,
  Array<T, DIM, SPACE>::ArrayIterator last)
{
  assert(first >= begin() && first < end());
  assert(last >= first && last <= end());

  // Empty range, return last
  if(first == last)
  {
    return last;
  }

  int count = 0;

  // Erase [first,last) elements
  while(first < last)
  {
    (*first).~T();
    first++;
    count++;
  }

  first -= count;
  int shifted = 0;

  // Shift [last, end) elements over
  while(last < end())
  {
    *first = *last;
    first++;
    last++;
    shifted++;
  }

  updateNumElements(m_num_elements - count);
  return first - shifted;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename... Args>
inline void Array<T, DIM, SPACE>::emplace(IndexType pos, Args&&... args)
{
  reserveForInsert(1, pos);
  m_data[pos] = std::move(T(std::forward<Args>(args)...));
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename... Args>
inline typename Array<T, DIM, SPACE>::ArrayIterator Array<T, DIM, SPACE>::emplace(
  Array<T, DIM, SPACE>::ArrayIterator pos,
  Args&&... args)
{
  assert(pos >= begin() && pos <= end());
  emplace(pos - begin(), args...);
  return pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
template <typename... Args>
inline void Array<T, DIM, SPACE>::resize(Args... args)
{
  static_assert(sizeof...(Args) == DIM,
                "Array size must match number of dimensions");
  // Intel hits internal compiler error when casting as part of function call
  const IndexType tmp_args[] = {args...};
  assert(detail::allNonNegative(tmp_args));
  const auto new_num_elements = detail::packProduct(tmp_args);

  static_cast<ArrayBase<T, DIM, Array<T, DIM, SPACE>>&>(*this) =
    ArrayBase<T, DIM, Array<T, DIM, SPACE>> {static_cast<IndexType>(args)...};

  if(new_num_elements > m_capacity)
  {
    dynamicRealloc(new_num_elements);
  }

  updateNumElements(new_num_elements);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::swap(Array<T, DIM, SPACE>& other)
{
  ArrayBase<T, DIM, Array<T, DIM, SPACE>>::swap(other);
  T* temp_data = m_data;
  IndexType temp_num_elements = m_num_elements;
  IndexType temp_capacity = m_capacity;
  double temp_resize_ratio = m_resize_ratio;

  m_data = other.m_data;
  m_num_elements = other.m_num_elements;
  m_capacity = other.m_capacity;
  m_resize_ratio = other.m_resize_ratio;

  other.m_data = temp_data;
  other.m_num_elements = temp_num_elements;
  other.m_capacity = temp_capacity;
  other.m_resize_ratio = temp_resize_ratio;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::initialize(IndexType num_elements,
                                             IndexType capacity)
{
  assert(num_elements >= 0);

  if(capacity < 0 || num_elements > capacity)
  {
    capacity = 0;
  }

  if(capacity == 0)
  {
    capacity = (num_elements > MIN_DEFAULT_CAPACITY) ? num_elements
                                                     : MIN_DEFAULT_CAPACITY;
  }
  setCapacity(capacity);
  updateNumElements(num_elements);

  // quick checks
  assert(m_data != nullptr);
  assert(m_num_elements >= 0);
  assert(m_capacity >= m_num_elements);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline T* Array<T, DIM, SPACE>::reserveForInsert(IndexType n, IndexType pos)
{
  assert(n >= 0);
  assert(pos >= 0);
  assert(pos <= m_num_elements);

  if(n == 0)
  {
    return m_data + pos;
  }

  IndexType new_size = m_num_elements + n;
  if(new_size > m_capacity)
  {
    dynamicRealloc(new_size);
  }

  T* const insert_pos = m_data + pos;
  T* cur_pos = m_data + m_num_elements - 1;
  for(; cur_pos >= insert_pos; --cur_pos)
  {
    *(cur_pos + n) = *cur_pos;
  }

  updateNumElements(new_size);
  return insert_pos;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::updateNumElements(IndexType new_num_elements)
{
  assert(new_num_elements >= 0);
  assert(new_num_elements <= m_capacity);
  m_num_elements = new_num_elements;
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::setCapacity(IndexType new_capacity)
{
  assert(new_capacity >= 0);

  if(new_capacity < m_num_elements)
  {
    updateNumElements(new_capacity);
  }

  m_data = axom::reallocate<T>(m_data, new_capacity, m_allocator_id);
  m_capacity = new_capacity;

  assert(m_data != nullptr || m_capacity <= 0);
}

//------------------------------------------------------------------------------
template <typename T, int DIM, MemorySpace SPACE>
inline void Array<T, DIM, SPACE>::dynamicRealloc(IndexType new_num_elements)
{
  assert(m_resize_ratio >= 1.0);
  IndexType new_capacity = new_num_elements * m_resize_ratio + 0.5;
  const IndexType block_size = this->blockSize();
  const IndexType remainder = new_capacity % block_size;
  if(remainder != 0)
  {
    new_capacity += block_size - remainder;
  }

  if(m_resize_ratio < 1.0)
  {
    std::cerr << "ERROR: resize ratio must be greater than 1.0.\n";
    std::cerr << "Set a valid resize ratio via calling setResizeRatio() with "
              << "an appropriate value.\n";

    utilities::processAbort();
  }

  m_data = axom::reallocate<T>(m_data, new_capacity, m_allocator_id);
  m_capacity = new_capacity;

  assert(m_data != nullptr || m_capacity <= 0);
}

} /* namespace axom */

#endif /* AXOM_ARRAY_HPP_ */
