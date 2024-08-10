// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_
#define AXOM_MIR_DISPATCH_RECTILINEAR_TOPOLOGY_HPP_

#include "axom/mir/views/StructuredTopologyView.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"
#include "axom/mir/views/dispatch_utilities.hpp"
#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint_mesh_utils.hpp>

namespace axom
{
namespace mir
{
namespace views
{
template <int NDIMS>
struct make_rectilinear
{ };

/**
 * \brief Create a 3D structured topology view with normal structured indexing.
 */
template <>
struct make_rectilinear<3>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 3>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topo.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(*coordset);
    const conduit::Node &values = coordset->fetch_existing("values");
    LogicalIndex zoneDims;
    zoneDims[0] =
      values.fetch_existing(axes[0]).dtype().number_of_elements() - 1;
    zoneDims[1] =
      values.fetch_existing(axes[1]).dtype().number_of_elements() - 1;
    zoneDims[2] =
      values.fetch_existing(axes[2]).dtype().number_of_elements() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topo.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Create a 2D structured topology view with normal structured indexing.
 */
template <>
struct make_rectilinear<2>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 2>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(*coordset);
    const conduit::Node &values = coordset->fetch_existing("values");
    LogicalIndex zoneDims;
    zoneDims[0] =
      values.fetch_existing(axes[0]).dtype().number_of_elements() - 1;
    zoneDims[1] =
      values.fetch_existing(axes[1]).dtype().number_of_elements() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Create a 1D structured topology view with normal structured indexing.
 */
template <>
struct make_rectilinear<1>
{
  using Indexing = views::StructuredIndexing<axom::IndexType, 1>;
  using LogicalIndex = typename Indexing::LogicalIndex;
  using TopoView = views::StructuredTopologyView<Indexing>;

  /**
   * \brief Create the indexing and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The indexing.
   */
  static Indexing indexing(const conduit::Node &topo)
  {
    const conduit::Node *coordset =
      conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    SLIC_ASSERT(coordset != nullptr);
    const auto axes = conduit::blueprint::mesh::utils::coordset::axes(*coordset);
    const conduit::Node &values = coordset->fetch_existing("values");
    LogicalIndex zoneDims;
    zoneDims[0] =
      values.fetch_existing(axes[0]).dtype().number_of_elements() - 1;
    return Indexing(zoneDims);
  }

  /**
   * \brief Create the topology view and initialize it from the topology.
   * \param topo The node containing the topology.
   * \return The topology view.
   */
  static TopoView view(const conduit::Node &topo)
  {
    return TopoView(indexing(topo));
  }
};

/**
 * \brief Creates a topology view compatible with rectilinear topologies and passes that view to the supplied function.
 *
 * \tparam FuncType The function/lambda type to invoke on the view.
 * \tparam SelectedDimensions  An integer whose bits indicate which dimensions are set.
 *
 * \param topo     The node that contains the rectilinear topology.
 * \param func     The function to invoke using the view.
 */
template <int SelectedDimensions = select_dimensions(1, 2, 3), typename FuncType>
void dispatch_rectilinear_topology(const conduit::Node &topo, FuncType &&func)
{
  const conduit::Node *coordset =
    conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
  SLIC_ASSERT(coordset != nullptr);
  const auto axes = conduit::blueprint::mesh::utils::coordset::axes(*coordset);
  switch(axes.size())
  {
  case 3:
    if constexpr(dimension_selected(SelectedDimensions, 3))
    {
      auto topoView = make_rectilinear<3>::view(topo);
      const std::string shape("hex");
      func(shape, topoView);
    }
    break;
  case 2:
    if constexpr(dimension_selected(SelectedDimensions, 2))
    {
      auto topoView = make_rectilinear<2>::view(topo);
      const std::string shape("quad");
      func(shape, topoView);
    }
    break;
  case 1:
    if constexpr(dimension_selected(SelectedDimensions, 1))
    {
      auto topoView = make_rectilinear<1>::view(topo);
      const std::string shape("line");
      func(shape, topoView);
    }
    break;
  default:
    break;
  }
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
