// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_UTILITIES_SORTING_HPP
#define AXOM_CORE_UTILITIES_SORTING_HPP
#include <axom/core/utilities/Utilities.hpp>
#include <axom/core/NumericLimits.hpp>

namespace axom
{
namespace utilities
{
/**
 * \brief This is a template suitable for sorting small arrays on device.
 *
 * \tparam T The data type being sorted.
 * \tparam N The biggest array size that will be used.
 *
 * \note For very short arrays, a simpler sort is faster. As the array size
 *       increases, the algorithm switches to qsort. Also, this is designed as
 *       a template class so it can be specialized.
 */
template <typename T, int N = axom::numeric_limits<int>::max()>
struct Sorting
{
  /**
   * \brief Sort an array of values in place.
   *
   * \param[inout] values The array of values to sort.
   * \param n The number of values in the array. 
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values, int n)
  {
    if(n < 11)
      insertionSort(values, n);
    else
      qsort(values, n);
  }

  /**
   * \brief Computes stack size for qsort.
   * \return A number of elements for an array-based stack.
   */
  AXOM_HOST_DEVICE
  constexpr static int stack_size()
  {
    int v = N;
    int i = 0;
    while(v > 0)
    {
      i++;
      v /= 2;
    }
    return i * 2;
  }

  /**
   * \brief Sort the input array using qsort.
   *
   * \param[inout] values The array to be sorted.
   * \param n The number of values in the array.
   */
  AXOM_HOST_DEVICE
  static void qsort(T *values, int n)
  {
    if(n <= 1) return;
    int stack[stack_size()][2];
    int stack_count = 1;
    stack[0][0] = 0;
    stack[0][1] = n - 1;

    while(stack_count > 0)
    {
      stack_count--;
      int low = stack[stack_count][0];
      int high = stack[stack_count][1];

      if(low < high)
      {
        int pivot_index = partition(values, low, high);

        stack[stack_count][0] = low;
        stack[stack_count][1] = pivot_index - 1;
        stack_count++;

        stack[stack_count][0] = pivot_index + 1;
        stack[stack_count][1] = high;
        stack_count++;
      }
    }
  }

  /**
   * \brief Create a partition for qsort.
   *
   * \param[inout] values The array to be sorted.
   * \param low The lower bound of the input partition.
   * \param high The upper bound of the input partition.
   *
   * \return A new pivot.
   */
  AXOM_HOST_DEVICE
  static int partition(T *values, int low, int high)
  {
    const T pivot = values[high];
    int i = low - 1;
    for(int j = low; j < high; j++)
    {
      if(values[j] <= pivot)
      {
        i++;
        axom::utilities::swap(values[i], values[j]);
      }
    }
    axom::utilities::swap(values[i + 1], values[high]);
    return i + 1;
  }

  /**
   * \brief Sort the input array using insertion sort.
   *
   * \param[inout] values The array to be sorted.
   * \param n The number of values in the array.
   */
  AXOM_HOST_DEVICE
  inline static void insertionSort(T* values, int n)
  {
    for(int i = 1; i < n; i++)
    {
      int j = i;
      // Keep swapping elements until we're not out-of-order.
      while(j > 0 && (values[j] < values[j - 1]))
      {
        axom::utilities::swap(values[j], values[j - 1]);
        j--;
      }
    }
  }
};

namespace detail
{

/**
 * \brief Swap 2 elements if b < a.
 * param a The first value.
 * param b The second value.
 */
template <typename T>
AXOM_HOST_DEVICE inline void ifswap(T &a, T &b)
{
  if(b < a)
  {
    T tmp = a;
    a = b;
    b = tmp;
  }
}

} // end namespace detail

/**
 * \brief Template specialization for sorting arrays with 3 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <typename T>
struct Sorting<T, 3>
{
  /**
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values)
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
  }
};

/**
 * \brief Template specialization for sorting arrays with 4 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <typename T>
struct Sorting<T, 4>
{
  /**
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values)
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
  }
};

}  // end namespace utilities
}  // end namespace axom

#endif
