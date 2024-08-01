// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_UNSTRUCTURED_TOPOLOGY_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/UnstructuredTopologySingleShapeView.hpp"
#include "axom/mir/views/UnstructuredTopologyPolyhedralView.hpp"
#include "axom/mir/views/UnstructuredTopologyMixedShapeView.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
#include "axom/mir/views/Shapes.hpp"

#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace mir
{
namespace views
{
// Turn on all bits so all shapes will be enabled.
constexpr int AnyShape = -1;

/**
 * \brief This function dispatches a Conduit polyhedral unstructured topology.
 *
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
 */
template <typename FuncType>
void dispatch_unstructured_polyhedral_topology(const conduit::Node &topo,
                                               FuncType &&func)
{
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "polyhedral")
  {
    IndexNode_to_ArrayView_same(
      topo["subelements/connectivity"],
      topo["subelements/sizes"],
      topo["subelements/offsets"],
      topo["elements/connectivity"],
      topo["elements/sizes"],
      topo["elements/offsets"],
      [&](auto seConnView,
          auto seSizesView,
          auto seOffsetsView,
          auto connView,
          auto sizesView,
          auto offsetsView) {
        using ConnType = typename decltype(seConnView)::value_type;
        UnstructuredTopologyPolyhedralView<ConnType> ugView(seConnView,
                                                            seSizesView,
                                                            seOffsetsView,
                                                            connView,
                                                            sizesView,
                                                            offsetsView);
        func(shape, ugView);
      });
  }
}

/**
 * \brief This function dispatches a Conduit mixed unstructured topology.
 *
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
 *
 * \note When this function makes the view, the view keeps a reference to
 *       the shape_map within the topology so we can build our own shape map
 *       later in the for_all_zones method.
 */
template <typename FuncType>
void dispatch_unstructured_mixed_topology(const conduit::Node &topo,
                                          FuncType &&func)
{
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "mixed")
  {
    IndexNode_to_ArrayView_same(
      topo["elements/connectivity"],
      topo["elements/shapes"],
      topo["elements/sizes"],
      topo["elements/offsets"],
      [&](auto connView, auto shapesView, auto sizesView, auto offsetsView) {
        using ConnType = typename decltype(connView)::value_type;

        UnstructuredTopologyMixedShapeView<ConnType> ugView(topo,
                                                            connView,
                                                            shapesView,
                                                            sizesView,
                                                            offsetsView);
        func(shape, ugView);
      });
  }
}

/**
 * \brief This function dispatches a Conduit topology to the right view type
 *        and passes that view to the supplied function/lambda.
 *
 * \tparam ShapeTypes Allows us to limit which shape types get compiled in.
 * \tparam FuncType The function/lambda type that will be invoked on the view.
 *
 * \param topo The node that contains the topology.
 * \param func The function/lambda to call with the topology view.
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
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Polyhedron_ShapeID))
    {
      if(shape == "polyhedral")
      {
        dispatch_unstructured_polyhedral_topology(topo, func);
        eligible = false;
      }
    }

#if 0
// TODO: Can't use polygon with single shape view because its sizes are not known at compile time.
      if constexpr (axom::utilities::bitIsSet(ShapeTypes, Polygon_ShapeID))
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
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Mixed_ShapeID))
    {
      if(eligible && shape == "mixed")
      {
        dispatch_unstructured_mixed_topology(topo, func);
        eligible = false;
      }
    }

    IndexNode_to_ArrayView(topo["elements/connectivity"], [&](auto connView) {
      using ConnType = typename decltype(connView)::value_type;
      // TODO: points, lines
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Tri_ShapeID))
      {
        if(eligible && shape == "tet")
        {
          UnstructuredTopologySingleShapeView<TriShape<ConnType>> ugView(connView);
          func(shape, ugView);
          eligible = false;
        }
      }
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Quad_ShapeID))
      {
        if(eligible && shape == "tet")
        {
          UnstructuredTopologySingleShapeView<QuadShape<ConnType>> ugView(
            connView);
          func(shape, ugView);
          eligible = false;
        }
      }
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Tet_ShapeID))
      {
        if(eligible && shape == "tet")
        {
          UnstructuredTopologySingleShapeView<TetShape<ConnType>> ugView(connView);
          func(shape, ugView);
          eligible = false;
        }
      }
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Pyramid_ShapeID))
      {
        if(eligible && shape == "pyramid")
        {
          UnstructuredTopologySingleShapeView<PyramidShape<ConnType>> ugView(
            connView);
          func(shape, ugView);
          eligible = false;
        }
      }
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Wedge_ShapeID))
      {
        if(eligible && shape == "wedge")
        {
          UnstructuredTopologySingleShapeView<WedgeShape<ConnType>> ugView(
            connView);
          func(shape, ugView);
          eligible = false;
        }
      }
      if constexpr(axom::utilities::bitIsSet(ShapeTypes, Hex_ShapeID))
      {
        if(eligible && shape == "hex")
        {
          UnstructuredTopologySingleShapeView<HexShape<ConnType>> ugView(connView);
          func(shape, ugView);
          eligible = false;
        }
      }
    });
  }
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
