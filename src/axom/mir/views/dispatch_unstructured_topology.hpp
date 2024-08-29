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
#include "axom/mir/blueprint_utilities.hpp"

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

template <typename ConnType, typename FuncType>
void typed_dispatch_unstructured_polyhedral_topology(const conduit::Node &topo,
                                                     FuncType &&func)
{
  namespace bputils = axom::mir::utilities::blueprint;
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "polyhedral")
  {
    auto seConnView =
      bputils::make_array_view<ConnType>(topo["subelements/connectivity"]);
    auto seSizesView =
      bputils::make_array_view<ConnType>(topo["subelements/sizes"]);
    auto seOffsetsView =
      bputils::make_array_view<ConnType>(topo["subelements/offsets"]);
    auto connView =
      bputils::make_array_view<ConnType>(topo["elements/connectivity"]);
    auto sizesView = bputils::make_array_view<ConnType>(topo["elements/sizes"]);
    auto offsetsView =
      bputils::make_array_view<ConnType>(topo["elements/offsets"]);

    UnstructuredTopologyPolyhedralView<ConnType> ugView(seConnView,
                                                        seSizesView,
                                                        seOffsetsView,
                                                        connView,
                                                        sizesView,
                                                        offsetsView);
    func(shape, ugView);
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

template <typename ConnType, typename FuncType>
void typed_dispatch_unstructured_mixed_topology(const conduit::Node &topo,
                                                FuncType &&func)
{
  namespace bputils = axom::mir::utilities::blueprint;
  const std::string shape = topo["elements/shape"].as_string();
  if(shape == "mixed")
  {
    auto connView =
      bputils::make_array_view<ConnType>(topo["elements/connectivity"]);
    auto shapesView =
      bputils::make_array_view<ConnType>(topo["elements/shapes"]);
    auto sizesView = bputils::make_array_view<ConnType>(topo["elements/sizes"]);
    auto offsetsView =
      bputils::make_array_view<ConnType>(topo["elements/offsets"]);

    UnstructuredTopologyMixedShapeView<ConnType> ugView(topo,
                                                        connView,
                                                        shapesView,
                                                        sizesView,
                                                        offsetsView);
    func(shape, ugView);
  }
}

template <typename... Args>
constexpr int encode_shapes(Args... args)
{
  return (... | args);
}

template <typename... Args>
constexpr int select_shapes(Args... args)
{
  return encode_types((1 << args)...);
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
template <typename ConnType, int ShapeTypes = AnyShape, typename FuncType>
void typed_dispatch_unstructured_topology(const conduit::Node &topo,
                                          FuncType &&func)
{
  namespace bputils = axom::mir::utilities::blueprint;
  const std::string type = topo["type"].as_string();
  if(type == "unstructured")
  {
    const std::string shape = topo["elements/shape"].as_string();
    const auto connView =
      bputils::make_array_view<ConnType>(topo["elements/connectivity"]);
    bool eligible = true;

    // Conditionally add polyhedron support.
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Polyhedron_ShapeID))
    {
      if(shape == "polyhedral")
      {
        typed_dispatch_unstructured_polyhedral_topology<ConnType>(topo, func);
        eligible = false;
      }
    }

    // TODO: add polygon

    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Mixed_ShapeID))
    {
      if(eligible && shape == "mixed")
      {
        typed_dispatch_unstructured_mixed_topology<ConnType>(topo, func);
        eligible = false;
      }
    }

    // TODO: points, lines

    // Make sizes / offsets views if the values are present.
    axom::ArrayView<ConnType> sizesView, offsetsView;
    if(topo.has_path("elements/sizes"))
      sizesView = bputils::make_array_view<ConnType>(
        topo.fetch_existing("elements/sizes"));
    if(topo.has_path("elements/offsets"))
      offsetsView = bputils::make_array_view<ConnType>(
        topo.fetch_existing("elements/offsets"));

    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Tri_ShapeID))
    {
      if(eligible && shape == "tri")
      {
        UnstructuredTopologySingleShapeView<TriShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Quad_ShapeID))
    {
      if(eligible && shape == "quad")
      {
        UnstructuredTopologySingleShapeView<QuadShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Tet_ShapeID))
    {
      if(eligible && shape == "tet")
      {
        UnstructuredTopologySingleShapeView<TetShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Pyramid_ShapeID))
    {
      if(eligible && shape == "pyramid")
      {
        UnstructuredTopologySingleShapeView<PyramidShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Wedge_ShapeID))
    {
      if(eligible && shape == "wedge")
      {
        UnstructuredTopologySingleShapeView<WedgeShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
    if constexpr(axom::utilities::bitIsSet(ShapeTypes, Hex_ShapeID))
    {
      if(eligible && shape == "hex")
      {
        UnstructuredTopologySingleShapeView<HexShape<ConnType>> ugView(
          connView,
          sizesView,
          offsetsView);
        func(shape, ugView);
        eligible = false;
      }
    }
  }
}

/// Dispatch in a way that does not care about the connectivity type.
template <int ShapeTypes = AnyShape, typename FuncType>
void dispatch_unstructured_topology(const conduit::Node &topo, FuncType &&func)
{
  IndexNode_to_ArrayView(topo["elements/connectivity"], [&](auto connView) {
    using ConnType = typename decltype(connView)::value_type;
    typed_dispatch_unstructured_topology<ConnType, ShapeTypes>(topo, func);
  });
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
