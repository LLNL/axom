// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * \brief A few utility functions used by the SLAM component.
 */
#ifndef SLAM_UTILITIES_H_
#define SLAM_UTILITIES_H_

#include "axom/core/Types.hpp"

#include <string>

namespace axom
{
namespace slam
{
using DefaultPositionType = axom::IndexType;
using DefaultElementType = axom::IndexType;

class NotImplementedException
{ };

namespace util
{
/** \brief A helper class to print the name of a few types */
template <typename T>
struct TypeToString
{
  static std::string to_string() { return "<unspecialized>"; }
};

/** \brief A helper class to print the name of integers as 'int' */
template <>
struct TypeToString<int>
{
  static std::string to_string() { return "int"; }
};

/** \brief A helper class to print the name of doubles as 'double' */
template <>
struct TypeToString<double>
{
  static std::string to_string() { return "double"; }
};

}  // end namespace util
}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_UTILITIES_H_
