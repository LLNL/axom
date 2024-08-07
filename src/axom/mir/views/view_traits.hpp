// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEW_TRAITS_HPP_
#define AXOM_MIR_VIEW_TRAITS_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/// General traits for topology views.
template <typename TopologyView>
struct view_traits
{
  static constexpr bool supports_strided_structured() { return false; }
};

/// If StructuredTopologyView was instantiated with StridedStructuredIndexing
/// (of varying dimensions) then say that strided structured is supported.
template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 3>>>
{
  static constexpr bool supports_strided_structured() { return true; }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 2>>>
{
  static constexpr bool supports_strided_structured() { return true; }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 1>>>
{
  static constexpr bool supports_strided_structured() { return true; }
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
