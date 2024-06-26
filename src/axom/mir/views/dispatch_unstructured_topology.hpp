// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_

#include "axom/mir/views/UnstructuredTopologySingleShapeView.hpp"
#include "axom/mir/views/UnstructuredTopologyPolyhedralView.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
#include "axom/mir/views/Shapes.hpp"

#include <conduit/conduit_blueprint.hpp>

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
      const std::string shape = topo["elements/shape"].as_string();
      bool eligible = true;

      // Conditionally add polyhedron support.
      if constexpr (ShapeTypes & UnstructuredTopologyPolyhedralView<IndexType>::PolyhedronShape::id())
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
              func(shape, ugView);
            });
          eligible = false;
        }
      }
#if 0
// TODO: Can't use polygon with single shape view because its sizes are not known at compile time.
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
            func(shape, ugView);
          });
          eligible = false;
        }
      }
#endif
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
        if constexpr (ShapeTypes & TriShape<IndexType>::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<TriShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & QuadShape<IndexType>::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<QuadShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & TetShape<IndexType>::id())
        {
          if(eligible && shape == "tet")
          {
            UnstructuredTopologySingleShapeView<TetShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & PyramidShape<IndexType>::id())
        {
          if(eligible && shape == "pyramid")
          {
            UnstructuredTopologySingleShapeView<PyramidShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & WedgeShape<IndexType>::id())
        {
          if(eligible && shape == "wedge")
          {
            UnstructuredTopologySingleShapeView<WedgeShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }
        if constexpr (ShapeTypes & HexShape<IndexType>::id())
        {
          if(eligible && shape == "hex")
          {
            UnstructuredTopologySingleShapeView<HexShape<IndexType>> ugView(connView);
            func(shape, ugView);
            eligible = false;
          }
        }

        // TODO: handle mixed shapes.

      });
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
