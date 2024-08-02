// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_STRUCTURED_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"
#include "axom/mir/views/StridedStructuredIndexing.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include "axom/mir/views/dispatch_uniform_topology.hpp"
#include "axom/mir/views/dispatch_rectilinear_topology.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
namespace views
{
/**
 * \brief Fill an array from a Conduit node, filling the destination array if the values do not exist.
 *
 * \tparam ArrayType The array type to use.
 *
 * \param n The conduit node that contains the named array, if it exists.
 * \param key The name of the node.
 * \param[out] The array to be filled.
 * \param fillValue the value to use if the array is not found.
 */
template <typename ArrayType>
bool fillFromNode(const conduit::Node &n,
                  const std::string &key,
                  ArrayType &arr,
                  int fillValue)
{
  bool found = false;
  if((found = n.has_path(key)) == true)
  {
    const auto acc = n.fetch_existing(key).as_int_accessor();
    for(int i = 0; i < arr.size(); i++) arr[i] = acc[i];
  }
  else
  {
    for(int i = 0; i < arr.size(); i++) arr[i] = fillValue;
  }
  return found;
}

/**
 * \brief Base template for strided structured topology creation
 */
template <int NDIMS>
struct make_strided_structured {};

/**
 * \brief Create a 3D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<3>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
  const std::string offsetsKey("elements/offsets");
  const std::string stridesKey("elements/strides");

      LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();
      zoneDims[1] = topo.fetch_existing("elements/dims/j").as_int();
      zoneDims[2] = topo.fetch_existing("elements/dims/k").as_int();

      LogicalIndex offsets, strides;
      fillFromNode(topo, offsetsKey, offsets, 0);
      if(!fillFromNode(topo, stridesKey, strides, 1))
      {
        strides[1] = zoneDims[0];
        strides[2] = zoneDims[0] * zoneDims[1];
      }

      Indexing zoneIndexing(zoneDims, offsets, strides);
    return TopoView(zoneIndexing);
  }
};

/**
 * \brief Create a 2D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<2>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
  const std::string offsetsKey("elements/offsets");
  const std::string stridesKey("elements/strides");
    LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();
      zoneDims[1] = topo.fetch_existing("elements/dims/j").as_int();

        LogicalIndex offsets, strides;
        fillFromNode(topo, offsetsKey, offsets, 0);
        if(!fillFromNode(topo, stridesKey, strides, 1))
        {
          strides[1] = zoneDims[0];
        }

        Indexing zoneIndexing(
          zoneDims,
          offsets,
          strides);
    return TopoView(zoneIndexing);

  }
};

/**
 * \brief Create a 1D structured topology view with strided structured indexing.
 */
template <>
struct make_strided_structured<1>
{
  using Indexing = views::StridedStructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
  const std::string offsetsKey("elements/offsets");
  const std::string stridesKey("elements/strides");

      LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();

        LogicalIndex offsets, strides;
        fillFromNode(topo, offsetsKey, offsets, 0);
        fillFromNode(topo, stridesKey, strides, 1);

        Indexing zoneIndexing(
          zoneDims,
          offsets,
          strides);
        return TopoView(zoneIndexing);
  }
};

/**
 * \brief Base template for structured topology creation
 */
template <int NDIMS>
struct make_structured {};

/**
 * \brief Create a 3D structured topology view with normal structured indexing.
 */
template <>
struct make_structured<3>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
      LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();
      zoneDims[1] = topo.fetch_existing("elements/dims/j").as_int();
      zoneDims[2] = topo.fetch_existing("elements/dims/k").as_int();

      Indexing zoneIndexing(zoneDims);
    return TopoView(zoneIndexing);
  }
};

/**
 * \brief Create a 2D structured topology view with normal structured indexing.
 */
template <>
struct make_structured<2>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
      LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();
      zoneDims[1] = topo.fetch_existing("elements/dims/j").as_int();

      Indexing zoneIndexing(zoneDims);
    return TopoView(zoneIndexing);
  }
};

/**
 * \brief Create a 1D structured topology view with normal structured indexing.
 */
template <>
struct make_structured<1>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
      LogicalIndex zoneDims;
      zoneDims[0] = topo.fetch_existing("elements/dims/i").as_int();

      Indexing zoneIndexing(zoneDims);
    return TopoView(zoneIndexing);
  }
};

/**
 * \brief Creates a topology view compatible with structured topologies and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_structured_topology(const conduit::Node &topo, FuncType &&func)
{
  int ndims = 1;
  ndims += topo.has_path("elements/dims/j") ? 1 : 0;
  ndims += topo.has_path("elements/dims/k") ? 1 : 0;
  const std::string offsetsKey("elements/offsets");
  const std::string stridesKey("elements/strides");

  switch(ndims)
  {
  case 3:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      const std::string shape("hex");
      if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
      {
        auto topoView = make_strided_structured<3>::view(topo);
        func(shape, topoView);
      }
      else
      {
        auto topoView = make_structured<3>::view(topo);
        func(shape, topoView);
      }
    }
    break;
  case 2:
    if constexpr(dimension_selected(SelectedDimensions, 2))
    {
      const std::string shape("quad");
      if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
      {
        auto topoView = make_strided_structured<2>::view(topo);
        func(shape, topoView);
      }
      else
      {
        auto topoView = make_structured<2>::view(topo);
        func(shape, topoView);
      }
    }
    break;
  case 1:
    if constexpr(dimension_selected(SelectedDimensions, 1))
    {
      const std::string shape("line");
      if(topo.has_path(offsetsKey) || topo.has_path(stridesKey))
      {
        auto topoView = make_strided_structured<1>::view(topo);
        func(shape, topoView);
      }
      else
      {
        auto topoView = make_structured<1>::view(topo);
        func(shape, topoView);
      }
    }
  }
}

/**
 * \brief Creates a topology view compatible with various logically "structured" topologies (uniform, rectilinear, structured) and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set. dimension
 *
 * \param topo     The node that contains the topology.
 * \param coordset The coordset node that contains the topology dimensions.
 * \param func     The function to invoke using the view. It should accept a string with the shape name and an auto parameter for the view.

 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_structured_topologies(const conduit::Node &topo,
                                    const conduit::Node &coordset,
                                    FuncType &&func)
{
  const std::string type = topo["type"].as_string();
  if(type == "uniform")
    dispatch_uniform_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "rectilinear")
    dispatch_rectilinear_topology<SelectedDimensions>(topo, coordset, func);
  else if(type == "structured")
    dispatch_structured_topology<SelectedDimensions>(topo, func);
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
