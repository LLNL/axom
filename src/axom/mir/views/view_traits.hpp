// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEW_TRAITS_HPP_
#define AXOM_MIR_VIEW_TRAITS_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/Shapes.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/// Return which shapes to expect for a given dimension.
static constexpr int shapes_for_dimension(int dimension)
{
  int shapes = 0;
  switch(dimension)
  {
  case 2:
    axom::utilities::setBitOn(shapes, Tri_ShapeID);
    axom::utilities::setBitOn(shapes, Quad_ShapeID);
    axom::utilities::setBitOn(shapes, Mixed_ShapeID);
    break;
  case 3:
    axom::utilities::setBitOn(shapes, Tet_ShapeID);
    axom::utilities::setBitOn(shapes, Pyramid_ShapeID);
    axom::utilities::setBitOn(shapes, Wedge_ShapeID);
    axom::utilities::setBitOn(shapes, Hex_ShapeID);
    axom::utilities::setBitOn(shapes, Mixed_ShapeID);
    break;
  }
  return shapes;
}

/// General traits for topology views.
template <typename TopologyView>
struct view_traits
{
  static constexpr bool supports_strided_structured() { return false; }
  static constexpr int selected_shapes()
  {
    return shapes_for_dimension(TopologyView::dimension());
  }
};

/// If StructuredTopologyView was instantiated with StridedStructuredIndexing
/// (of varying dimensions) then say that strided structured is supported.
template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 3>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(3); }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 2>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(2); }
};

template <typename IndexT>
struct view_traits<StructuredTopologyView<StridedStructuredIndexing<IndexT, 1>>>
{
  static constexpr bool supports_strided_structured() { return true; }
  static constexpr int selected_shapes() { return shapes_for_dimension(1); }
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
