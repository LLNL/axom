// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ARRAY_HPP_
#define AXOM_ARRAY_HPP_

#include "axom/config.hpp"                    // for compile-time defines
#include "axom/core/Macros.hpp"               // for axom macros
#include "axom/core/memory_management.hpp"    // for memory allocation functions
#include "axom/core/utilities/Utilities.hpp"  // for processAbort()
#include "axom/core/Types.hpp"                // for IndexType definition

// C/C++ includes
#include <cstring>   // for std::memcpy
#include <iostream>  // for std::cerr

namespace axom
{
/* Provided so that 0 doesn't convert to nullptr and lead to ambiguous
 * constructor calls. */
namespace internal
{
constexpr IndexType ZERO = 0;
}

/*!
 * \class Array
 *
 * \brief Provides a generic multi-component array container.
 *
 *  The Array class provides a generic multi-component array container with
 *  dynamic re-allocation and insertion. Each element in the array is a tuple
 *  consisting of 1 or more components, which are stored contiguously.
 *
 *  Depending on which constructor is used, the Array object can have two
 *  different underlying storage types:
 *
 *  * <b> Native Storage </b> <br />
 *
 *     When using native storage, the Array object manages all memory.
 *     Typically, the Array object will allocate extra space to facilitate
 *     the insertion of new elements and minimize the number of reallocations.
 *     The actual capacity of the array (i.e., total number of tuples that the
 *     Array can hold) can be queried by calling the capacity() function.
 *     When allocated memory is used up, inserting a new element triggers a
 *     re-allocation.  At each re-allocation, extra space is allocated
 *     according to the <em> resize_ratio </em> parameter, which is set to 2.0
 *     by default. To return all extra memory, an application can call
 *     `shrink()`.
 *
 *     \note The Array destructor deallocates and returns all memory associated
 *      with it to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    An Array object may be constructed from an external, user-supplied buffer
 *    consisting of the given number of tuples and specified number of
 *    components per tuple.  In this case, the Array object does not own the
 *    memory.  Instead, the Array object makes a shallow copy of the pointer.
 *
 *    \warning An Array object that points to an external buffer has a fixed
 *     size and cannot be dynamically resized.
 *
 *    \note The Array destructor does not deallocate a user-supplied buffer,
 *     since it does not manage that memory.
 *
 * \warning Reallocations tend to be costly operations in terms of performance.
 *  Use `reserve()` when the number of nodes is known a priori, or opt to
 *  use a constructor that takes an actual size and capacity when possible.
 *
 * \tparam T the type of the values to hold.
 *
 */
template <typename T>
class Array
{
public:
  static constexpr double DEFAULT_RESIZE_RATIO = 2.0;
  static constexpr IndexType MIN_DEFAULT_CAPACITY = 32;

public:
  /// \name Native Storage Array Constructors
  /// @{

  /*!
   * \brief Constructs an Array instance with the given number of tuples.
   *
   * \param [in] num_tuples the number of tuples the Array holds.
   * \param [in] num_components the number of values per tuple. If not
   *  specified defaults to 1.
   * \param [in] capacity the number of tuples to allocate space for.
   *
   * \note If no capacity or capacity less than num_tuples is specified
   *  then it will default to at least num_tuples * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of tuples to store in the
   *  array and does not correspond to the actual bytesize.
   *
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   *
   * \post capacity() >= size()
   * \post size() == num_tuples
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  Array(IndexType num_tuples, IndexType num_components = 1, IndexType capacity = 0);

  /// @}

  /// \name External Storage Array Constructors
  /// @{

  /*!
   * \brief Constructs an Array instance with the given number of tuples from
   *  an external data buffer.
   *
   * \param [in] data the external data this Array will wrap.
   * \param [in] num_tuples the number of tuples in the Array.
   * \param [in] num_components the number of values per tuple. If not
   *  specified defaults to 1.
   * \param [in] capacity the capacity of the external buffer.
   *
   * \pre data != nullptr
   * \pre num_tuples > 0
   * \pre num_components >= 1
   *
   * \post numComponents() == num_components
   * \post getResizeRatio == 0.0
   *
   * \note a capacity is specified for the number of tuples to store in the
   *  array and does not correspond to the actual bytesize.
   * \note If no capacity or capacity less than num_tuples is specified then
   *  it will default to the number of tuples.
   *
   * \note This constructor wraps the supplied buffer and does not own the data.
   *  Consequently, the Array instance cannot be reallocated.
   */
  Array(T* data,
        IndexType num_tuples,
        IndexType num_components = 1,
        IndexType capacity = 0);

  /// @}

  /*!
   * Destructor. Frees the associated buffer unless the memory is external.
   */
  virtual ~Array();

  /// \name Array tuple access operators
  /// @{

  /*!
   * \brief Accessor, returns a reference to the given component of the
   *  specified tuple.
   *
   * \param [in] pos the tuple to query.
   * \param [in] component the component to return.
   *
   * \pre 0 <= pos < size()
   * \pre 0 <= component < numComponents()
   */
  /// @{

  inline T& operator()(IndexType pos, IndexType component = 0)
  {
    assert(inBounds(pos, component));

    return m_data[pos * m_num_components + component];
  }

  inline const T& operator()(IndexType pos, IndexType component = 0) const
  {
    assert(inBounds(pos, component));

    return m_data[pos * m_num_components + component];
  }

  /// @}

  /*!
   * \brief Accessor, returns a reference to the given value.
   *
   * \param [in] idx the position of the value to return.
   *
   * \note equivalent to *(array.getData() + idx).
   *
   * \pre 0 <= idx < m_num_tuples * m_num_components
   */
  /// @{

  T& operator[](IndexType idx)
  {
    assert(inBounds(idx));

    return m_data[idx];
  }

  const T& operator[](IndexType idx) const
  {
    assert(inBounds(idx));

    return m_data[idx];
  }

  /// @}

  /*!
   * \brief Return a pointer to the array of data.
   */
  /// @{

  T* getData() { return m_data; }
  const T* getData() const { return m_data; }

  /// @}

  /// @}

  /// \name Array methods to modify the data.
  /// @{

  /*!
   * \brief Set all the values of the array.
   *
   * \param [in] value the value to set to.
   */
  void fill(const T& value)
  {
    std::fill_n(m_data, m_num_tuples * m_num_components, value);
  }

  /*!
   * \brief Append a value to the end of the array.
   *
   * \param [in] value the value to append.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   *
   * \pre m_num_components == 1.
   */
  void append(const T& value);

  /*!
   * \brief Append tuples to the end of the array.
   *
   * \param [in] tuples the tuples to append.
   * \param [in] n the number of tuples to append.
   *
   * \note It's assumed that tuples is of length n * m_num_components.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void append(const T* tuples, IndexType n);

  /*!
   * \brief Modify the values of existing tuples.
   *
   * \param [in] tuples the new tuples to write.
   * \param [in] n the number of tuples to write.
   * \param [in] pos the position at which to begin writing.
   *
   * \note It's assumed that tuples is of length n * m_num_components.
   * \note The size is unchanged by calls to set.
   *
   * \pre pos + n <= m_num_tuples.
   */
  void set(const T* tuples, IndexType n, IndexType pos);

  /*!
   * \brief Insert a tuple into the array at the given position.
   *
   * \param [in] value the value to insert.
   * \param [in] pos the position at which to insert.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by 1.
   *
   * \pre numComponents() == 1.
   */
  void insert(const T& value, IndexType pos);

  /*!
   * \brief Insert tuples into the array at the given position.
   *
   * \param [in] tuples the tuples to insert.
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   *
   * \note It's assumed that tuples is of length n * m_num_components.
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   *
   * \pre pos <= m_num_tuples.
   */
  void insert(const T* tuples, IndexType n, IndexType pos);

  /*!
   * \brief Insert multiple copies of the same value at the given position.
   *
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position to insert at.
   * \param [in] value the value for each component of the new tuples. If not
   *  specified defaults to the default value of T (zero for most numeric
   *  types).
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   * \note The size increases by n.
   * \note This method is used to create space for tuples in the middle of the
   *  array.
   *
   * \pre pos <= m_num_tuples.
   */
  void emplace(IndexType n, IndexType pos, const T& value = T());

  /// @}

  /// \name Array methods to query and set attributes
  /// @{

  /*!
   * \brief Return the number of tuples allocated for the data array.
   */
  IndexType capacity() const { return m_capacity; }

  /*!
   * \brief Increase the capacity. Does nothing if the new capacity is less
   *  than the current capacity.
   *
   * \param [in] capacity the new number of tuples to allocate.
   */
  void reserve(IndexType capacity)
  {
    if(capacity > m_capacity)
    {
      setCapacity(capacity);
    }
  }

  /*!
   * \brief Shrink the capacity to be equal to the size.
   */
  void shrink() { setCapacity(m_num_tuples); }

  /*!
   * \brief Returns true iff the Array stores no elements.
   *
   * \note If the Array is empty the capacity can still be greater than zero.
   */
  bool empty() const { return m_num_tuples == 0; }

  /*!
   * \brief Return the number of tuples stored in the data array.
   */
  IndexType size() const { return m_num_tuples; }

  /*!
   * \brief Update the number of tuples stored in the data array.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  void resize(IndexType new_num_tuples);

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
   * \brief Return the number of components per tuple.
   */
  IndexType numComponents() const { return m_num_components; }

  /*!
   * \brief Return true iff the external buffer constructor was called.
   */
  bool isExternal() const { return m_is_external; }

  /*!
   * \brief Return true iff a sidre constructor was called.
   */
  virtual bool isInSidre() const { return false; }

  /// @}

protected:
  /*! \brief Default constructor supports infrastructure in subclasses. */
  Array();

  /*!
   * \brief Initialize an Array instance with the given number of tuples.
   *
   * \param [in] num_tuples the number of tuples the Array holds.
   * \param [in] num_components the number of values per tuple. If not
   *  specified defaults to 1.
   * \param [in] capacity the number of tuples to allocate space for.
   *
   * \note If no capacity or capacity less than num_tuples is specified
   *  then it will default to at least num_tuples * DEFAULT_RESIZE_RATIO.
   * \note a capacity is specified for the number of tuples to store in the
   *  array and does not correspond to the actual bytesize.
   *
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   *
   * \post capacity() >= size()
   * \post size() == num_tuples
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  void initialize(IndexType num_tuples,
                  IndexType num_components,
                  IndexType capacity);

  /*!
   * \brief Make space for a subsequent insertion into the array.
   *
   * \param [in] n the number of tuples to insert.
   * \param [in] pos the position at which to begin the insertion.
   *
   * \return a pointer to the beginning of the insertion space.
   *
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  T* reserveForInsert(IndexType n, IndexType pos);

  /*!
   * \brief Update the number of tuples.
   *
   * \param [in] new_num_tuples the new number of tuples.
   */
  virtual void updateNumTuples(IndexType new_num_tuples);

  /*!
   * \brief Set the number of tuples allocated for the data array.
   *
   * \param [in] capacity the new number of tuples to allocate.
   */
  virtual void setCapacity(IndexType new_capacity);

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_tuples the number of tuples which exceeds the current
   *  capacity.
   */
  virtual void dynamicRealloc(IndexType new_num_tuples);

  /// \name Internal bounds-checking routines
  /// @{

  /*! \brief Test if pos and component are within bounds */
  inline bool inBounds(IndexType pos, IndexType component) const
  {
    return (pos >= 0 && pos < m_num_tuples) &&
      (component >= 0 && component < m_num_components);
  }

  /*! \brief Test if idx is within bounds */
  inline bool inBounds(IndexType idx) const
  {
    return idx >= 0 && idx < m_num_tuples * m_num_components;
  }
  /// @}

  T* m_data;
  IndexType m_num_tuples;
  IndexType m_capacity;
  IndexType m_num_components;
  double m_resize_ratio;
  bool const m_is_external;

  DISABLE_COPY_AND_ASSIGNMENT(Array);
  DISABLE_MOVE_AND_ASSIGNMENT(Array);
};

//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T>
Array<T>::Array()
  : m_data(nullptr)
  , m_num_tuples(0)
  , m_capacity(0)
  , m_num_components(1)
  , m_resize_ratio(DEFAULT_RESIZE_RATIO)
  , m_is_external(false)
{ }

//------------------------------------------------------------------------------
template <typename T>
Array<T>::Array(IndexType num_tuples, IndexType num_components, IndexType capacity)
  : m_data(nullptr)
  , m_num_tuples(0)
  , m_capacity(0)
  , m_num_components(0)
  , m_resize_ratio(DEFAULT_RESIZE_RATIO)
  , m_is_external(false)
{
  initialize(num_tuples, num_components, capacity);
}

//------------------------------------------------------------------------------
template <typename T>
Array<T>::Array(T* data,
                IndexType num_tuples,
                IndexType num_components,
                IndexType capacity)
  : m_data(data)
  , m_num_tuples(num_tuples)
  , m_capacity(0)
  , m_num_components(num_components)
  , m_resize_ratio(0.0)
  , m_is_external(true)
{
  m_capacity = (capacity < num_tuples) ? num_tuples : capacity;

  assert(m_num_tuples >= 0);
  assert(m_num_components >= 1);
  assert(m_num_tuples <= m_capacity);
  assert(m_data != nullptr || m_capacity <= 0);
}

//------------------------------------------------------------------------------
template <typename T>
Array<T>::~Array()
{
  if(m_data != nullptr && !m_is_external)
  {
    axom::deallocate(m_data);
  }

  m_data = nullptr;
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::append(const T& value)
{
  assert(m_num_components == 1);

  IndexType new_size = m_num_tuples + 1;
  if(new_size > m_capacity)
  {
    dynamicRealloc(new_size);
  }

  m_data[m_num_tuples] = value;
  updateNumTuples(new_size);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::append(const T* tuples, IndexType n)
{
  IndexType new_size = m_num_tuples + n;
  if(new_size > m_capacity)
  {
    dynamicRealloc(new_size);
  }

  T* cur_end = m_data + m_num_tuples * m_num_components;
  std::memcpy(cur_end, tuples, n * m_num_components * sizeof(T));
  updateNumTuples(new_size);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::set(const T* tuples, IndexType n, IndexType pos)
{
  assert(tuples != nullptr);
  assert(pos >= 0);
  assert(pos + n <= m_num_tuples);

  T* set_position = &m_data[pos * m_num_components];
  IndexType byte_size = n * m_num_components * sizeof(T);
  std::memcpy(set_position, tuples, byte_size);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::insert(const T& value, IndexType pos)
{
  assert(m_num_components == 1);
  reserveForInsert(1, pos);
  m_data[pos] = value;
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::insert(const T* tuples, IndexType n, IndexType pos)
{
  assert(tuples != nullptr);
  T* insert_pos = reserveForInsert(n, pos);
  std::memcpy(insert_pos, tuples, n * m_num_components * sizeof(T));
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::emplace(IndexType n, IndexType pos, const T& value)
{
  T* insert_pos = reserveForInsert(n, pos);
  std::fill_n(insert_pos, n * numComponents(), value);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::resize(IndexType new_num_tuples)
{
  assert(new_num_tuples >= 0);

  if(new_num_tuples > m_capacity)
  {
    dynamicRealloc(new_num_tuples);
  }

  updateNumTuples(new_num_tuples);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::initialize(IndexType num_tuples,
                                 IndexType num_components,
                                 IndexType capacity)
{
  assert(num_tuples >= 0);
  assert(num_components > 0);

  m_num_tuples = num_tuples;
  m_num_components = num_components;

  if(capacity < 0 || num_tuples > capacity)
  {
    capacity = 0;
  }

  if(capacity == 0)
  {
    capacity =
      (num_tuples > MIN_DEFAULT_CAPACITY) ? num_tuples : MIN_DEFAULT_CAPACITY;
  }
  setCapacity(capacity);

  // sanity checks
  assert(m_data != nullptr);
  assert(m_num_tuples >= 0);
  assert(m_capacity >= m_num_tuples);
  assert(m_num_components >= 1);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* Array<T>::reserveForInsert(IndexType n, IndexType pos)
{
  assert(n >= 0);
  assert(pos >= 0);
  assert(pos <= m_num_tuples);

  if(n == 0)
  {
    return m_data + pos * m_num_components;
  }

  IndexType new_size = m_num_tuples + n;
  if(new_size > m_capacity)
  {
    dynamicRealloc(new_size);
  }

  T* const insert_pos = m_data + pos * m_num_components;
  T* cur_pos = m_data + (m_num_tuples * m_num_components) - 1;
  for(; cur_pos >= insert_pos; --cur_pos)
  {
    *(cur_pos + n * m_num_components) = *cur_pos;
  }

  updateNumTuples(new_size);
  return insert_pos;
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::updateNumTuples(IndexType new_num_tuples)
{
  assert(new_num_tuples >= 0);
  assert(new_num_tuples <= m_capacity);
  m_num_tuples = new_num_tuples;
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::setCapacity(IndexType new_capacity)
{
  assert(new_capacity >= 0);

  if(m_is_external && new_capacity <= m_capacity)
  {
    return;
  }

  if(m_is_external)
  {
    std::cerr << "Cannot reallocate an externally provided buffer.";
    utilities::processAbort();
  }

  if(new_capacity < m_num_tuples)
  {
    updateNumTuples(new_capacity);
  }

  m_data = axom::reallocate(m_data, new_capacity * m_num_components);
  m_capacity = new_capacity;

  assert(m_data != nullptr || m_capacity <= 0);
}

//------------------------------------------------------------------------------
template <typename T>
inline void Array<T>::dynamicRealloc(IndexType new_num_tuples)
{
  if(m_is_external)
  {
    std::cerr << "Cannot reallocate an externally provided buffer.";
    utilities::processAbort();
  }

  assert(m_resize_ratio >= 1.0);
  const IndexType new_capacity =
    static_cast<IndexType>(new_num_tuples * m_resize_ratio + 0.5);

  if(m_resize_ratio < 1.0)
  {
    std::cerr << "ERROR: resize ratio must be greater than 1.0.\n";
    std::cerr << "Set a valid resize ratio via calling setResizeRatio() with "
              << "an appropriate value.\n";

    utilities::processAbort();
  }

  m_data = axom::reallocate(m_data, new_capacity * m_num_components);
  m_capacity = new_capacity;

  assert(m_data != nullptr || m_capacity <= 0);
}

} /* namespace axom */

#endif /* AXOM_ARRAY_HPP_ */
