// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file StridePolicies.hpp
 *
 * \brief Stride policies for SLAM
 *
 * Stride policies are meant to represent the fixed distance between consecutive
 * elements of an OrderedSet
 * A valid stride policy must support the following interface:
 *   * [required]
 *   * DEFAULT_VALUE is a public static const IntType
 *   * stride() : IntType  -- returns the stride
 *   * isValid() : bool -- indicates whether the Stride policy of the set is
 *     valid
 *   * [optional]
 *   * operator(): IntType -- alternate accessor for the stride value
 *
 * \note All non-zero stride values are valid.
 */

#ifndef SLAM_POLICIES_STRIDE_H_
#define SLAM_POLICIES_STRIDE_H_

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slam
{
namespace policies
{
/**
 * \name OrderedSet_Stride_Policies
 * \brief A few default policies for the stride of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for the stride in a set.
 * When using this class, the stride can be set at runtime.
 */
template <typename IntType>
struct RuntimeStride
{
public:
  static const IntType DEFAULT_VALUE = IntType(1);

  RuntimeStride(IntType stride = DEFAULT_VALUE) : m_stride(stride) { }

  inline IntType stride() const { return m_stride; }
  inline IntType& stride() { return m_stride; }

  void setStride(IntType str) { m_stride = str; }

  inline IntType operator()() const { return stride(); }
  inline IntType& operator()() { return stride(); }

  /** All non-zero strides are valid     */
  inline bool isValid(bool) const { return (m_stride != 0); }

  //inline bool hasStride() const       { return m_stride != IntType(); }
private:
  IntType m_stride;
};

/**
 * \brief A policy class for a compile-time known stride
 */
template <typename IntType, IntType INT_VAL>
struct CompileTimeStride
{
  static const IntType DEFAULT_VALUE = INT_VAL;

  CompileTimeStride(IntType val = DEFAULT_VALUE) { setStride(val); }

  inline IntType stride() const { return INT_VAL; }
  inline IntType operator()() const { return stride(); }

  void setStride(IntType AXOM_DEBUG_PARAM(val))
  {
    SLIC_ASSERT_MSG(
      val == INT_VAL,
      "slam::CompileTimeStride -- tried to set a compile time stride"
        << " with value (" << val << " ) that differs from the template"
        << " parameter of " << INT_VAL << ".");
  }

  /** All non-zero strides are valid     */
  inline bool isValid(bool) const { return (INT_VAL != 0); }
};

/**
 * \brief A policy class for a set with stride one (i.e. the default stride)
 */
template <typename IntType>
struct StrideOne
{
  static const IntType DEFAULT_VALUE;

  /**
   * This constructor only exists to allow the derived class to not have
   * to specialize for when the stride is known at compile time
   */
  StrideOne(IntType val = DEFAULT_VALUE) { setStride(val); }

  inline const IntType stride() const { return DEFAULT_VALUE; }
  inline const IntType operator()() const { return stride(); }

  void setStride(IntType AXOM_DEBUG_PARAM(val))
  {
    SLIC_ASSERT_MSG(
      val == DEFAULT_VALUE,
      "slam::StrideOne policy -- tried to set a stride-one StridePolicy"
        << " with value (" << val << "), but should always be 1.");
  }
  inline bool isValid(bool) const { return true; }
};

template <typename IntType>
const IntType StrideOne<IntType>::DEFAULT_VALUE = IntType {1};

/// \}

}  // end namespace policies
}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_POLICIES_STRIDE_H_
