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

/*!
 * \brief Swap one set of values, move to next ones.
 *
 * \param idx1 The first array index.
 * \param idx2 The second array index.
 * \param values The current array values being swapped.
 * \param args A parameter pack of additional arrays that need to be swapped.
 */
template <typename T, typename... Args>
AXOM_HOST_DEVICE
inline void internal_swap(int idx1, int idx2, T values, Args... args)
{
  axom::utilities::swap(values[idx1], values[idx2]);
  internal_swap(idx1, idx2, args...);
}

/*!
 * \brief This method is called when there are no more values to swap.
 */
AXOM_HOST_DEVICE
inline void internal_swap(int, int)
{
}

template <typename T, typename... Args>
AXOM_HOST_DEVICE
inline void ifswap_multiple(int idx1, int idx2, T values, Args... args)
{
  if(values[idx2] < values[idx1])
  {
    internal_swap(idx1, idx2, values, args...);
  }
}

}  // end namespace detail

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
template <int N = axom::numeric_limits<int>::max()>
class Sorting
{
  struct Delimiter{};
public:
  static constexpr int MAX_SIZE = N;
  static constexpr int SORT_SIZE_CUTOFF = 11;

  /**
   * \brief Sort an array of values in place.
   *
   * \param[inout] values The array of values to sort.
   * \param n The number of values in the array. 
   */
  template <typename T>
  AXOM_HOST_DEVICE
  inline static void sort(T values, int n)
  {
    if(n < SORT_SIZE_CUTOFF)
      insertionSort(values, n);
    else
      qsort(values, n);
  }

  /*!
   * \brief Sort multiple arrays, using the first array as the key for sorting.
   *        All other arrays are sorted the same.
   *
   * \tparam T The type of the first array.
   * \tparam Args Parameter pack containing multiple array types and an int.
   *
   * \param first The first data array.
   * \param args A parameter pack of all other args.
   */
  AXOM_HOST_DEVICE
  template <typename T, typename... Args>
  inline static void sort_multiple(T first, Args... args)
  {
    sort_multiple_internal(first, Delimiter{}, args...);
  }

  //----------------------------------------------------------------------------
  // Internal methods (sort)
  //----------------------------------------------------------------------------

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
  template <typename T>
  AXOM_HOST_DEVICE
  static void qsort(T values, int n)
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
  template <typename T>
  AXOM_HOST_DEVICE
  static int partition(T values, int low, int high)
  {
    const auto pivot = values[high];
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
  template <typename T>
  AXOM_HOST_DEVICE
  inline static void insertionSort(T values, int n)
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

  //----------------------------------------------------------------------------
  // Internal methods (sort_multiple)
  //----------------------------------------------------------------------------

  /*!
   * \brief Shift args to end until the array length is at the desired position.
   */
  AXOM_HOST_DEVICE
  template <typename T, typename S, typename... Args>
  inline static void sort_multiple_internal(T first, Delimiter d, S second, Args... args)
  {
    sort_multiple_internal(first, d, args..., second);
  }

  /*!
   * \brief Terminal case for argument reordering where we know all of the
   *        array arguments as well as the array length argument.
   *
   * \param first The first data array used for sorting.
   * \param d Delimiter
   * \param n The length of the arrays.
   * \param args A parameter pack containing all other arrays to be sorted.
   */
  template <typename T, typename... Args>
  AXOM_HOST_DEVICE
  inline static void sort_multiple_internal(T first, Delimiter, int n, Args... args)
  {
    if(n < SORT_SIZE_CUTOFF)
    {
      insertionSortMultiple(n, first, args...);
    }
    else
    {
      qsortMultiple(n, first, args...);
    }
  }

  /**
   * \brief Sort the input array using qsort.
   *
   * \param[inout] values The array to be sorted.
   * \param n The number of values in the array.
   */
  template <typename T, typename... Args>
  AXOM_HOST_DEVICE
  static void qsortMultiple(int n, T values, Args... args)
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
        int pivot_index = partitionMultiple(low, high, values, args...);

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
  template <typename T, typename... Args>
  AXOM_HOST_DEVICE
  static int partitionMultiple(int low, int high, T values, Args... args)
  {
    const T pivot = values[high];
    int i = low - 1;
    for(int j = low; j < high; j++)
    {
      if(values[j] <= pivot)
      {
        i++;
        detail::internal_swap(i, j, values, args...);
      }
    }
    detail::internal_swap(i + 1, high, values, args...);
    return i + 1;
  }

  /*!
   * \brief Perform actual sort on multiple data arrays.
   *
   * \param n The length of the data arrays.
   * \param values The array being sorted.
   * \param args The rest of the arrays to be sorted the same way.
   */
  template <typename T, typename... Args>
  AXOM_HOST_DEVICE
  static void insertionSortMultiple(int n, T values, Args... args)
  {
    for(int i = 1; i < n; i++)
    {
      int j = i;
      // Keep swapping elements until we're not out-of-order.
      while(j > 0 && (values[j] < values[j - 1]))
      {
        detail::internal_swap(j, j - 1, values, args...);
        j--;
      }
    }
  }
};

/**
 * \brief Template specialization for sorting arrays with 3 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <>
struct Sorting<3>
{
  /**
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  template <typename T>
  inline static void sort(T values)
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
  }

  /**
   * \brief Sort the input arrays.
   *
   * \param[inout] values The first array to be sorted.
   * \param[inout] args A parameter pack of other arrays to sort.
   */
  AXOM_HOST_DEVICE
  template <typename T, typename... Args>
  inline static void sort_multiple(T values, Args... args)
  {
    detail::ifswap_multiple(0, 1, values, args...);
    detail::ifswap_multiple(1, 2, values, args...);
    detail::ifswap_multiple(0, 1, values, args...);
  }
};

/**
 * \brief Template specialization for sorting arrays with 4 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <>
struct Sorting<4>
{
  /**
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  template <typename T>
  inline static void sort(T values)
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
  }

  /**
   * \brief Sort the input arrays.
   *
   * \param[inout] values The first array to be sorted.
   * \param[inout] args A parameter pack of other arrays to sort.
   */
  AXOM_HOST_DEVICE
  template <typename T, typename... Args>
  inline static void sort_multiple(T values, Args... args)
  {
    detail::ifswap_multiple(0, 1, values, args...);
    detail::ifswap_multiple(2, 3, values, args...);
    detail::ifswap_multiple(1, 2, values, args...);
    detail::ifswap_multiple(0, 1, values, args...);
    detail::ifswap_multiple(2, 3, values, args...);
    detail::ifswap_multiple(1, 2, values, args...);
  }
};

}  // end namespace utilities
}  // end namespace axom

#endif
