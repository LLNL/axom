// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_COORDSET_HPP_
#define AXOM_MIR_DISPATCH_COORDSET_HPP_

#include "axom/mir/views/UnstructuredTopologySingleShapeView.hpp"
#include "axom/mir/views/UnstructuredTopologyPolyhedralView.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

namespace axom
{
namespace mir
{
namespace views
{

constexpr int AnyShape = -1;

/**
 * \brief This function dispatches a Conduit topology to the right view type
 *        and passes that view to the supplied function/lambda.
 *
 * \tparam ShapeTypes Allows us to limit which shape types get compiled in.
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \
 */
template <int ShapeTypes = AnyShape, typename FuncType>
void dispatch_unstructured_topology(const conduit::Node &topo, FuncType &&func)
{
  const std::string type = topo["type"].as_string();
  if(type == "unstructured")
  {
    const conduit::Node *coordset_ptr = conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    const conduit::Node &coordset = *coordset_ptr;
    dispatch_coordset(coordset, [&](auto coordsetView)
    {
      const std::string shape = topo["elements/shape"].as_string();
      bool eligible = true;

      // Conditionally add polyhedron support.
      if constexpr (ShapeTypes & UnstructuredTopologyPolyhedralView::PolyhedralShape::id())
      {
        if(shape == "polyhedral")
        {
          IndexNode_to_ArrayView_same(
            topo["subelements/connectivity"], topo["subelements/sizes"], topo["subelements/offsets"],
            topo["elements/connectivity"], topo["elements/sizes"], topo["elements/offsets"],
            [&](auto seConnView, auto seSizesView, auto seOffsetsView, auto connView, auto sizesView, auto offsetsView)
            {
              using IndexType = typename decltype(seConnView)::value_type;
              UnstructuredTopologyPolyhedralView<IndexType> ugView(seConnView, seSizesView, seOffsetsView, connView, sizesView, offsetsView);
              func(shape, ugView, coordsetView);
            });
          eligible = false;
        }
      }

      if constexpr (ShapeTypes & PolygonShape::id())
      {
        if(eligible && shape == "polygon")
        {
          IndexNode_to_ArrayView_same(topo["elements/connectivity"],
                                      topo["elements/sizes"],
                                      topo["elements/offsets"], [&](auto connView, auto sizesView, auto offsetsView)
          {
            using IndexType = typename decltype(connView)::value_type;
            UnstructuredTopologySingleShapeView<IndexType, PolygonShape<IndexType> > ugView(connView, sizesView, offsetsView);
            func(shape, ugView, coordsetView);
          });
          eligible = false;
        }
      }

      if constexpr (ShapeTypes & AnyShape)
      {
        if(eligible && shape == "mixed")
        {
          // TODO: handle mixed zone types.
          eligible = false;
        }
      }

      IndexNode_to_ArrayView(topo["elements/connectivity"], [&](auto connView)
      {
        using IndexType = typename decltype(connView)::value_type;
        // TODO: points, lines
        if constexpr (ShapeTypes & TriShape::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<IndexType, TriShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & QuadShape::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<IndexType, QuadShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & TetShape::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<IndexType, TetShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & PyramidShape::id())
        {
          if(eligible && shape == "pyramid")
          {
            UnstructuredTopologySingleShapeView<IndexType, PyramidShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & WedgeShape::id())
        {
          if(eligible && shape == "wedge")
          {
            UnstructuredTopologySingleShapeView<IndexType, WedgeShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & HexShape::id())
        {
          if(eligible && shape == "hex")
          {
            UnstructuredTopologySingleShapeView<IndexType, HexShape<IndexType> > ugView(connView);
            func(shape, ugView, coordsetView);
            eligible = false;
          }
        }
      });
    });
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
