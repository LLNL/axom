// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Array.hpp"               /* for axom::Array */
#include "axom/core/utilities/Utilities.hpp" /* for utilities::max */

#include "axom/slic/core/SimpleLogger.hpp" /* for SimpleLogger */
#include "axom/slic/interface/slic.hpp"    /* for slic macros */

#include "gtest/gtest.h" /* for TEST and EXPECT_* macros */

#include "axom/sidre/core/sidre.hpp"

// C/C++ includes
#include <algorithm> /* for std::fill_n */

namespace axom
{
namespace sidre
{
const char IGNORE_OUTPUT[] = ".*";

namespace internal
{
/*!
 * \brief Calculate the new capacity for and Array given an increase in the
 *  size.
 * \param [in] v, the Array in question.
 * \param [in] increase, the amount the size will increase by
 * \return the new capacity.
 */
template <typename T>
axom::IndexType calc_new_capacity(Array<T>& v, axom::IndexType increase)
{
  axom::IndexType new_num_tuples = v.size() + increase;
  if(new_num_tuples > v.capacity())
  {
    return new_num_tuples * v.getResizeRatio() + 0.5;
  }

  return v.capacity();
}

/*!
 * \brief Check if two Arrays are equivalent. Does not check the resize ratio.
 * \param [in] lhs, the first Array to compare.
 * \param [in] rhs, the second Array to compare.
 * \return the new capacity.
 */
template <typename T>
void check_copy(const Array<T>& lhs, const Array<T>& rhs)
{
  EXPECT_EQ(lhs.size(), rhs.size());
  EXPECT_EQ(lhs.numComponents(), rhs.numComponents());
  EXPECT_EQ(lhs.capacity(), rhs.capacity());

  const T* lhs_data = lhs.getData();
  const T* rhs_data = rhs.getData();
  EXPECT_EQ(lhs_data, rhs_data);
}

/*!
 * \brief Check that the storage of an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_storage(Array<T>& v)
{
  EXPECT_TRUE(v.empty());
  EXPECT_EQ(v.size(), 0);

  axom::IndexType capacity = v.capacity();
  axom::IndexType num_components = v.numComponents();
  const T* data_ptr = v.getData();

  /* Append up to half the capacity. */
  if(num_components == 1)
  {
    for(T i = 0; i < capacity / 2; ++i)
    {
      v.append(i);
    }
  }
  else
  {
    T* tuple = new T[num_components];
    for(axom::IndexType i = 0; i < capacity / 2; ++i)
    {
      for(axom::IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(tuple, 1);
    }
    delete[] tuple;
    tuple = nullptr;
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity / 2);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);

  /* Append up to the full capacity. */
  if(num_components == 1)
  {
    for(T i = capacity / 2; i < capacity; ++i)
    {
      v.append(i);
    }
  }
  else
  {
    T* tuple = new T[num_components];
    for(axom::IndexType i = capacity / 2; i < capacity; ++i)
    {
      for(axom::IndexType j = 0; j < num_components; ++j)
      {
        tuple[j] = i * num_components + j;
      }
      v.append(tuple, 1);
    }
    delete[] tuple;
    tuple = nullptr;
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);

  /* Check the array data using the () and [] operators and the raw pointer. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
      EXPECT_EQ(v[i * num_components + j], i * num_components + j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * num_components + j);
    }
  }

  /* Set the array data to new values using the () operator. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Check the array data using the () and [] operators and the raw pointer. */
  for(axom::IndexType i = 0; i < capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
      EXPECT_EQ(v[i * num_components + j], i * j - 5 * i + 7 * j);
      EXPECT_EQ(data_ptr[i * num_components + j], i * j - 5 * i + 7 * j);
    }
  }

  /* Check the array metadata. */
  EXPECT_TRUE(!v.empty());
  EXPECT_EQ(v.size(), capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_ptr);
}

/*!
 * \brief Check that the fill method is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_fill(Array<T>& v)
{
  constexpr T MAGIC_NUM_0 = 55;
  constexpr T MAGIC_NUM_1 = 6834;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType size = v.size();
  const axom::IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  /* Fill the Array with MAGIC_NUM_0. */
  v.fill(MAGIC_NUM_0);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the entries are all MAGIC_NUM_0. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM_0);
    }
  }

  /* Fill the Array with MAGIC_NUM_1. */
  v.fill(MAGIC_NUM_1);

  /* Check the meta data. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the entries are all MAGIC_NUM_1. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM_1);
    }
  }
}

/*!
 * \brief Check that the set method is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_set(Array<T>& v)
{
  constexpr T ZERO = 0;
  const axom::IndexType capacity = v.capacity();
  const axom::IndexType size = v.size();
  const axom::IndexType num_components = v.numComponents();
  const double ratio = v.getResizeRatio();
  const T* const data_ptr = v.getData();

  /* Allocate a buffer half the size of the array. Fill it up with sequential
   * values. */
  const axom::IndexType buffer_size = size / 2;
  T* buffer = new T[buffer_size * num_components];
  for(axom::IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i;
  }

  /* Set all the values in the array to zero. */
  v.fill(ZERO);

  /* Set the first half of the tuples in the array to the sequential values in
   * buffer. */
  v.set(buffer, buffer_size, 0);

  /* Check the array metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that the first half of the tuples in the array are equivalent to
   * those in buffer. */
  for(axom::IndexType i = 0; i < buffer_size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), buffer[i * num_components + j]);
    }
  }

  /* Check that the second half of the tuples in the array are all zero. */
  for(axom::IndexType i = buffer_size; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), ZERO);
    }
  }

  /* Reset the values in buffer to the next sequential values. */
  for(axom::IndexType i = 0; i < buffer_size * num_components; ++i)
  {
    buffer[i] = i + buffer_size * num_components;
  }

  /* Set the second half of the tuples in the array to the new sequential values
   * in buffer. */
  v.set(buffer, buffer_size, buffer_size);

  /* Check the array metadata. */
  EXPECT_EQ(capacity, v.capacity());
  EXPECT_EQ(size, v.size());
  EXPECT_EQ(num_components, v.numComponents());
  EXPECT_EQ(ratio, v.getResizeRatio());
  EXPECT_EQ(data_ptr, v.getData());

  /* Check that all the tuples in the array now hold sequential values. */
  for(axom::IndexType i = 0; i < 2 * buffer_size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }
}

/*!
 * \brief Check that the resizing of an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_resize(Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;
  axom::IndexType num_components = v.numComponents();

  /* Check that the size equals the capacity. */
  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  T* tuple = new T[num_components];
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.append(tuple, 1);
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Prepare 1000 tuples to be appended. */
  const axom::IndexType n_tuples = 1000;
  T* values = new T[n_tuples * num_components];
  for(axom::IndexType i = 0; i < n_tuples; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      axom::IndexType i_real = i + size;
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  /* Append the new tuples. */
  capacity = calc_new_capacity(v, n_tuples);
  v.append(values, n_tuples);
  size += n_tuples;

  /* Check that size and capacity are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is still intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reduce the size down to 500 tuples */
  T* data_address = v.getData();
  size = 500;
  v.resize(size);

  /* Check the metadata. */
  EXPECT_EQ(v.size(), size);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.getData(), data_address);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Shrink the vector */
  capacity = size;
  v.shrink();

  /* Check that the capacity and size are as expected. */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append a new tuple, should resize. */
  old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.append(tuple, 1);
  size++;

  /* Check the new size and capacity. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check that the data is intact. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Reset the data */
  T* data_ptr = v.getData();
  for(axom::IndexType i = 0; i < size * num_components; ++i)
  {
    data_ptr[i] = i;
  }

  /* Append a bunch of tuples to fill in up to the capacity. Resize should
   * not occur. */
  old_capacity = capacity;
  for(axom::IndexType i = size; i < old_capacity; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }

    v.append(tuple, 1);
    size++;
    EXPECT_EQ(v.capacity(), old_capacity);
    EXPECT_EQ(v.size(), size);
    EXPECT_EQ(v.getData(), data_ptr);
  }

  EXPECT_EQ(v.size(), old_capacity);

  /* Append a final tuple that should trigger a resize. */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * num_components + j;
  }

  capacity = calc_new_capacity(v, old_capacity - size + 1);
  v.append(tuple, 1);
  size++;

  /* Check the new capacity and size. */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  delete[] tuple;
  tuple = nullptr;

  delete[] values;
  values = nullptr;
}

/*!
 * \brief Check that the insertion into an Array is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_insert(Array<T>& v)
{
  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;
  axom::IndexType num_components = v.numComponents();

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * j - 5 * i + 7 * j;
    }
  }

  /* Append a new tuple, should resize. */
  axom::IndexType old_capacity = capacity;
  capacity = calc_new_capacity(v, 1);
  T* tuple = new T[num_components];
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    tuple[j] = size * j - 5 * size + 7 * j;
  }
  v.insert(tuple, 1, v.size());
  size++;

  /* Check that it resized properly */
  EXPECT_GT(capacity, old_capacity);
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  /* Append 1000 tuples */
  const axom::IndexType n_tuples = 1000;
  T* values = new T[n_tuples * num_components];
  for(axom::IndexType i = 0; i < n_tuples; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      axom::IndexType i_real = i + size;
      values[i * num_components + j] = i_real * j - 5 * i_real + 7 * j;
    }
  }

  capacity = calc_new_capacity(v, n_tuples);
  v.insert(values, n_tuples, size);
  size += n_tuples;

  /* Check that it resizes properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - 5 * i + 7 * j);
    }
  }

  capacity = size;
  v.shrink();
  axom::IndexType n_insert_front = 100;

  /* Reset the data */
  T* data_ptr = v.getData();
  for(axom::IndexType i = 0; i < size * num_components; ++i)
  {
    data_ptr[i] = i + num_components * n_insert_front;
  }

  /* Insert into the front of the Array. */
  for(axom::IndexType i = n_insert_front - 1; i >= 0; i--)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      tuple[j] = i * num_components + j;
    }
    capacity = calc_new_capacity(v, 1);
    v.insert(tuple, 1, 0);
    size++;
  }

  /* Check that the insertion worked as expected */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  delete[] tuple;
  tuple = nullptr;

  delete[] values;
  values = nullptr;
}

/*!
 * \brief Check that the emplace method is working properly.
 * \param [in] v the Array to check.
 */
template <typename T>
void check_emplace(Array<T>& v)
{
  constexpr T MAGIC_NUM = 52706;

  /* Resize the array up to the capacity */
  axom::IndexType capacity = v.capacity();
  v.resize(capacity);
  axom::IndexType size = capacity;
  const axom::IndexType num_components = v.numComponents();

  EXPECT_EQ(v.size(), v.capacity());

  /* Set the existing data in v */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * num_components + j;
    }
  }

  /* Emplace 1 tuple at the end of the array with the default value,
   * should resize. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, v.size());

  /* Emplace 9 tuples at the end of the array with MAGIC_NUM. */
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, v.size(), MAGIC_NUM);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the data */
  for(axom::IndexType i = 0; i < size - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * num_components + j);
    }
  }

  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(axom::IndexType i = size - 9; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 10 tuples at the beginning of the array with MAGIC_NUM. */
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, 0, MAGIC_NUM);

  /* Emplace 1 tuple at the beginning of the array with the default value. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, 0);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(axom::IndexType i = 1; i < 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the middle */
  for(axom::IndexType i = 10; i < size - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the end */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(axom::IndexType i = size - 9; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Emplace 9 tuples in the middle of the array with MAGIC_NUM. */
  axom::IndexType middle = size / 2;
  capacity = calc_new_capacity(v, 9);
  v.emplace(9, middle, MAGIC_NUM);

  /* Emplace 1 tuple in the middle of the array with the default value. */
  capacity = calc_new_capacity(v, 1);
  v.emplace(1, middle);
  size += 10;

  /* Check that it resized properly */
  EXPECT_EQ(v.capacity(), capacity);
  EXPECT_EQ(v.size(), size);

  /* Check the beginning */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(0, j), T());
  }

  for(axom::IndexType i = 1; i < 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the first section */
  for(axom::IndexType i = 10; i < middle; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 10) * num_components + j);
    }
  }

  /* Check the middle */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(middle, j), T());
  }

  for(axom::IndexType i = middle + 1; i < middle + 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }

  /* Check the second section */
  for(axom::IndexType i = middle + 10; i < size - 10; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), (i - 20) * num_components + j);
    }
  }

  /* Check the end */
  for(axom::IndexType j = 0; j < num_components; ++j)
  {
    EXPECT_EQ(v(size - 10, j), T());
  }

  for(axom::IndexType i = size - 9; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), MAGIC_NUM);
    }
  }
}

/*!
 * \brief Make a copy of an Array through sidre and check it for defects.
 * \param [in] v the Array to copy.
 */
template <typename T>
void check_sidre(Array<T>& v)
{
  ASSERT_TRUE(v.isInSidre());
  /* Create a copy. */
  Array<T> cpy(const_cast<sidre::View*>(v.getView()));
  cpy.setResizeRatio(v.getResizeRatio());

  /* Check that the copy holds the same data. */
  check_copy(v, cpy);

  /* Resize the copy and check that it functions correctly. */
  cpy.resize(0);
  check_storage(cpy);
  check_insert(cpy);
}

/*!
 * \brief Check an external array for defects.
 * \param [in] v the external array to check.
 */
template <typename T>
void check_external(Array<T>& v)
{
  ASSERT_TRUE(v.isExternal());

  /* Check that the array is full. */
  ASSERT_EQ(v.size(), v.capacity());

  const axom::IndexType size = v.size();
  const axom::IndexType num_components = v.numComponents();
  const axom::IndexType num_values = size * num_components;
  T* const data_ptr = v.getData();

  /* Set the tuples in the array. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      v(i, j) = i * num_components + j;
    }
  }

  /* Check the tuples using the raw pointer. */
  for(axom::IndexType i = 0; i < num_values; ++i)
  {
    EXPECT_EQ(data_ptr[i], i);
  }

  /* Set the tuples using the raw pointer. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      data_ptr[i * num_components + j] = i * j - i - j;
    }
  }

  /* Check the tuples using the () operator. */
  for(axom::IndexType i = 0; i < size; ++i)
  {
    for(axom::IndexType j = 0; j < num_components; ++j)
    {
      EXPECT_EQ(v(i, j), i * j - i - j);
    }
  }

  EXPECT_EQ(size, v.size());
  EXPECT_EQ(size, v.capacity());
  EXPECT_EQ(data_ptr, v.getData());

  /* Since the array is full all of the following calls should require a
   * reallocation and cause a fatal error. */
  T* tuple = new T[num_components];
  EXPECT_DEATH_IF_SUPPORTED(v.append(tuple, 1), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(v.insert(tuple, 1, 0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(v.reserve(size + 1), IGNORE_OUTPUT);

  delete[] tuple;
  tuple = nullptr;
}

} /* end namespace internal */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkStorage)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      Array<int> v_int_sidre(root->createView("int"), ZERO, n_components, capacity);
      internal::check_storage(v_int_sidre);

      Array<double> v_double_sidre(root->createView("double"),
                                   ZERO,
                                   n_components,
                                   capacity);
      internal::check_storage(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkFill)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      Array<int> v_int_sidre(root->createView("int"), size, n_components, capacity);
      internal::check_fill(v_int_sidre);

      Array<double> v_double_sidre(root->createView("double"),
                                   size,
                                   n_components,
                                   capacity);
      internal::check_fill(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkSet)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  for(axom::IndexType capacity = 2; capacity < 512; capacity *= 2)
  {
    axom::IndexType size = capacity / 2;
    for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
    {
      Array<int> v_int_sidre(root->createView("int"), size, n_components, capacity);
      internal::check_set(v_int_sidre);

      Array<double> v_double_sidre(root->createView("double"),
                                   size,
                                   n_components,
                                   capacity);
      internal::check_set(v_double_sidre);

      root->destroyViewsAndData();
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkResize)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 4; n_components++)
      {
        Array<int> v_int_sidre(root->createView("int"),
                               ZERO,
                               n_components,
                               capacity);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_resize(v_int_sidre);

        Array<double> v_double_sidre(root->createView("double"),
                                     ZERO,
                                     n_components,
                                     capacity);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_resize(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array_DeathTest, checkResize)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* Resizing isn't allowed with a ratio less than 1.0. */
  Array<int> v_int(root->createView("int"), internal::ZERO, 1, 100);
  v_int.setResizeRatio(0.99);
  EXPECT_DEATH_IF_SUPPORTED(internal::check_resize(v_int), IGNORE_OUTPUT);
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkInsert)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        Array<int> v_int_sidre(root->createView("int"),
                               ZERO,
                               n_components,
                               capacity);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_insert(v_int_sidre);

        Array<double> v_double_sidre(root->createView("double"),
                                     ZERO,
                                     n_components,
                                     capacity);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_insert(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkEmplace)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        Array<int> v_int_sidre(root->createView("int"),
                               ZERO,
                               n_components,
                               capacity);
        v_int_sidre.setResizeRatio(ratio);
        internal::check_emplace(v_int_sidre);

        Array<double> v_double_sidre(root->createView("double"),
                                     ZERO,
                                     n_components,
                                     capacity);
        v_double_sidre.setResizeRatio(ratio);
        internal::check_emplace(v_double_sidre);

        root->destroyViewsAndData();
      }
    }
  }
}

/* Sidre specific tests */

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkSidre)
{
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        Array<int> v_int(root->createView("int"), ZERO, n_components, capacity);
        v_int.setResizeRatio(ratio);
        internal::check_storage(v_int);
        internal::check_sidre(v_int);

        Array<double> v_double(root->createView("double"),
                               ZERO,
                               n_components,
                               capacity);
        v_double.setResizeRatio(ratio);
        internal::check_storage(v_double);
        internal::check_sidre(v_double);

        root->destroyViewsAndData();
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(sidre_core_array, checkSidrePermanence)
{
  constexpr double MAGIC_NUM = 5683578.8;

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr axom::IndexType ZERO = 0;

  for(double ratio = 1.0; ratio <= 2.0; ratio += 0.5)
  {
    for(axom::IndexType capacity = 2; capacity <= 512; capacity *= 2)
    {
      for(axom::IndexType n_components = 1; n_components <= 3; n_components++)
      {
        const double* array_data_ptr;
        axom::IndexType num_values;
        { /* Begin scope */
          Array<double> v(root->createView("double"), ZERO, n_components, capacity);
          array_data_ptr = v.getData();
          num_values = v.size() * v.numComponents();
          v.setResizeRatio(ratio);
          internal::check_storage(v);

          /* Set v's data to MAGIC_NUM */
          for(axom::IndexType i = 0; i < v.size(); ++i)
          {
            for(axom::IndexType j = 0; j < v.numComponents(); ++j)
            {
              v(i, j) = MAGIC_NUM;
            }
          }
        } /* End scope, v has been deallocated. */

        /* Check that the data still exists in sidre */
        sidre::View* view = root->getView("double");
        const double* view_data_ptr =
          static_cast<const double*>(view->getVoidPtr());
        EXPECT_EQ(view_data_ptr, array_data_ptr);
        EXPECT_EQ(view->getNumDimensions(), 2);

        IndexType dims[2];
        view->getShape(2, dims);
        EXPECT_EQ(dims[0], capacity);
        EXPECT_EQ(dims[1], n_components);

        for(axom::IndexType i = 0; i < num_values; ++i)
        {
          EXPECT_EQ(view_data_ptr[i], MAGIC_NUM);
        }

        root->destroyViewsAndData();
      }
    }
  }
}

} /* end namespace sidre */
} /* end namespace axom */

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
