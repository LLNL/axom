/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file OffsetPolicies.hpp
 *
 * \brief Offset policies for SLAM
 *
 * Offset policies are meant to represent the offset to the first element of
 * SLAM ordered set.
 *
 * A valid offset policy must support the following interface:
 *   * [required]
 *     * DEFAULT_VALUE is a public static constant of type IntType
 *     * offset() : IntType  -- returns the offset
 *     * isValid() : bool -- indicates whether the Offset policy of the set is
 *       valid  [optional]
 *     * operator(): IntType -- alternate accessor for the offset value
 */

#ifndef SLAM_POLICIES_OFFSET_H_
#define SLAM_POLICIES_OFFSET_H_

#include "axom/core/Macros.hpp"


namespace axom
{
namespace slam
{
namespace policies
{

/**
 * \name OrderedSet_Offset_Policies
 * \brief A few default policies for the offset of an OrderedSet
 */

/// \{

/**
 * \brief A policy class for the offset in a set.  The offset can be set at
 * runtime.
 */
template<typename IntType>
struct RuntimeOffset
{
public:
  static const IntType DEFAULT_VALUE = IntType();

  RuntimeOffset(IntType off = DEFAULT_VALUE) : m_off(off) {}

  inline IntType          offset() const { return m_off; }
  inline IntType&         offset() { return m_off; }

  inline IntType operator ()() const { return offset(); }
  inline IntType& operator()() { return offset(); }

  inline bool             isValid(bool) const { return true; }
private:
  IntType m_off;
};


/**
 * \brief A policy class for a compile-time known set offset
 */
template<typename IntType, IntType INT_VAL>
struct CompileTimeOffset
{
  static const IntType DEFAULT_VALUE = INT_VAL;

  CompileTimeOffset(IntType val = DEFAULT_VALUE) {
    AXOM_DEBUG_VAR(val);
    SLIC_ASSERT_MSG(
      val == INT_VAL,
      "slam::CompileTimeOffset -- tried to initialize a compile time "
      << "offset with value (" << val << " ) that differs from "
      << "the template parameter of " << INT_VAL << ".");
  }

  inline IntType          offset() const { return INT_VAL; }

  inline IntType operator ()() const { return offset(); }

  inline bool             isValid(bool) const { return true; }
};

/**
 * \brief A policy class for when we have no offset
 */
template<typename IntType>
struct ZeroOffset
{
  static const IntType DEFAULT_VALUE = IntType();

  ZeroOffset(IntType val = DEFAULT_VALUE)
  {
    AXOM_DEBUG_VAR(val);
    SLIC_ASSERT_MSG(
      val == DEFAULT_VALUE,
      "slam::ZeroOffset policy -- tried to initialize a NoOffset policy"
      << " with (" << val << ", but should always be 0");
  }

  inline IntType          offset() const { return DEFAULT_VALUE; }

  inline IntType operator ()() const { return offset(); }

  inline bool             isValid(bool) const { return true; }
};

/// \}

} // end namespace policies
} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICIES_OFFSET_H_
