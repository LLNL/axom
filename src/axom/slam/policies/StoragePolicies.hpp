// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_POLICIES_STORAGE_HPP
#define SLAM_POLICIES_STORAGE_HPP

#include "axom/core/Array.hpp"
#include "axom/core/memory_management.hpp"

#include <vector>

namespace axom
{
namespace slam
{
namespace policies
{
/*!
 * \brief A policy class for storage backed by an std::vector<T>.
 */
template <typename T>
struct STLVectorStorage
{
  using ElementType = T;
  using StorageType = std::vector<T>;

  static StorageType create(int nelems,
                            T defaultValue = {},
                            int allocatorID = axom::getDefaultAllocatorID())
  {
    AXOM_UNUSED_VAR(allocatorID);
    return StorageType(nelems, defaultValue);
  }
};

/*!
 * \brief A policy class for storage backed by an axom::Array<T>.
 */
template <typename T>
struct ArrayStorage
{
  using ElementType = T;
  using StorageType = axom::Array<T>;

  static StorageType create(int nelems,
                            T defaultValue = {},
                            int allocatorID = axom::getDefaultAllocatorID())
  {
    StorageType buf(nelems, nelems, allocatorID);
    buf.fill(defaultValue);
    return buf;
  }
};

}  // namespace policies
}  // namespace slam
}  // namespace axom

#endif  // SLAM_POLICIES_STORAGE_HPP
