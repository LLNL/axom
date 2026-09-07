// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file StridePolicies.hpp
 *
 * \brief Stride policies for SLAM
 *
 * In an OrderedSet, stride() is the distance between consecutive indices before
 * indirection. Scalar policies accept any nonzero stride, including a negative
 * stride for reverse traversal.
 *
 * In a Map, stride() is the number of components per set element and must be
 * positive. MultiDimStride describes a multidimensional component shape with
 * positive dimensions and a product representable by its index type.
 *
 * RuntimeStride and CompileTimeStride use the storage and validity checks in
 * ValuePolicies.hpp. Concepts.hpp defines the operations required by each owner.
 */

#include "axom/core/Macros.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/slam/policies/ValuePolicies.hpp"
#include "axom/slam/detail/SizeChecks.hpp"

namespace axom::slam::policies
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
struct RuntimeStride : RuntimeValue<StrideTag<IntType>>
{
private:
  using BaseType = RuntimeValue<StrideTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE;
  static const bool IS_COMPILE_TIME = false;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return StrideTag<IntType>::defaultValue(); }

  using BaseType::BaseType;

  /// \brief Returns the stride between consecutive elements.
  AXOM_HOST_DEVICE constexpr IntType stride() const { return this->value(); }
  AXOM_HOST_DEVICE constexpr IntType& stride() { return this->value(); }

  /*!
   * \brief Returns the shape of the inner data for a given stride.
   *  This only has meaning when used with Map-based types.
   */
  AXOM_HOST_DEVICE constexpr IntType shape() const { return this->value(); }

  void setStride(IntType str) { this->m_value = str; }
};

template <typename IntType>
const IntType RuntimeStride<IntType>::DEFAULT_VALUE = StrideTag<IntType>::defaultValue();

/// \brief A policy class for a compile-time known stride
template <typename IntType, IntType INT_VAL>
struct CompileTimeStride : CompileTimeValue<INT_VAL, StrideTag<IntType>>
{
private:
  using BaseType = CompileTimeValue<INT_VAL, StrideTag<IntType>>;

public:
  static const IntType DEFAULT_VALUE = INT_VAL;
  static const bool IS_COMPILE_TIME = true;
  constexpr static int NumDims = 1;

  using IndexType = IntType;
  using ShapeType = IntType;

  static constexpr IntType DefaultSize() { return DEFAULT_VALUE; }

  using BaseType::BaseType;

  AXOM_HOST_DEVICE constexpr IntType stride() const { return INT_VAL; }
  AXOM_HOST_DEVICE constexpr IntType shape() const { return INT_VAL; }

  AXOM_HOST_DEVICE void setStride(IntType AXOM_DEBUG_PARAM(val))
  {
    SLIC_ASSERT_MSG(val == INT_VAL,
                    "slam::CompileTimeStride -- tried to set a compile time stride"
                      << " with value (" << val << " ) that differs from the template"
                      << " parameter of " << INT_VAL << ".");
  }
};

/// \brief A policy with stride one.
template <typename IntType>
using StrideOne = CompileTimeStride<IntType, 1>;

/// \brief A row-major multidimensional component shape for maps.
template <typename IntType, int Dims>
struct MultiDimStride
{
  static_assert(Dims > 0, "MultiDimStride requires at least one dimension");
  using IndexType = IntType;
  using ShapeType = StackArray<IntType, Dims>;
  constexpr static int NumDims = Dims;

  static ShapeType DefaultSize()
  {
    ShapeType array;
    for(int i = 0; i < Dims; i++)
    {
      array[i] = 1;
    }
    return array;
  }

  AXOM_HOST_DEVICE MultiDimStride(StackArray<IntType, Dims> shape) : m_shape(shape)
  {
    IntType product = 1;
    for(int i = Dims - 1; i >= 0; --i)
    {
      if(m_shape[i] <= 0 || !::axom::slam::detail::nonnegativeProductFits(product, m_shape[i]))
      {
#ifndef AXOM_DEVICE_CODE
        SLIC_ERROR(
          "MultiDimStride requires positive dimensions and a representable component count.");
#else
        SLIC_ASSERT_MSG(false, "Invalid MultiDimStride shape.");
#endif
        // Avoid invalid arithmetic even when error logging is configured not to abort.
        m_strides = {};
        return;
      }
      m_strides[i] = product;
      product *= m_shape[i];
    }
  }

  /// \brief Returns the total number of components, the product of the shape dimensions.
  AXOM_HOST_DEVICE inline IntType stride() const { return m_shape[0] * m_strides[0]; }

  inline IntType operator()() const { return stride(); }
  inline IntType operator()() { return stride(); }

  AXOM_HOST_DEVICE bool isValid(bool = false) const
  {
    IntType product = 1;
    for(int i = Dims - 1; i >= 0; --i)
    {
      if(m_shape[i] <= 0 || !::axom::slam::detail::nonnegativeProductFits(product, m_shape[i]) ||
         m_strides[i] != product)
      {
        return false;
      }
      product *= m_shape[i];
    }
    return true;
  }

  /// \brief Returns the strides for each indexing dimension.
  AXOM_HOST_DEVICE inline ShapeType strides() const { return m_strides; }
  /*!
   * \brief Returns the multi-dimensional shape of the inner data.
   *  This only has meaning when used with Map-based types.
   */
  AXOM_HOST_DEVICE inline ShapeType shape() const { return m_shape; }

private:
  ShapeType m_shape;
  ShapeType m_strides;
};

/// \}

}  // end namespace axom::slam::policies
