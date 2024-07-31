// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_
#define AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_

#include "axom/mir/views/Shapes.hpp"

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief This class provides a view for Conduit/Blueprint single shape unstructured grids.
 *
 * \tparam IndexT The index type that will be used for connectivity, etc.
 * \tparam ShapeT The shape type.
 */
template <typename ShapeT>
class UnstructuredTopologySingleShapeView
{
public:
  using ShapeType = ShapeT;
  using ConnectivityType = typename ShapeType::ConnectivityType;
  using ConnectivityView = typename ShapeType::ConnectivityView;

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn) : m_connectivity(conn), m_sizes(), m_offsets()
  {
  }

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn,
                                      const ConnectivityView &sizes,
                                      const ConnectivityView &offsets) : m_connectivity(conn), m_sizes(sizes), m_offsets(offsets)
  {
    SLIC_ASSERT(m_sizes.size() != 0);
    SLIC_ASSERT(m_offsets.size() != 0);
    SLIC_ASSERT(m_offsets.size() == m_sizes.size());
  }

  /**
   * \brief Return the dimension of the shape.
   *
   * \return The dimension of the shape.
   */
  static constexpr int dimension() { return ShapeT::dimension(); }

  /**
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  IndexType numberOfZones() const
  {
    return (m_sizes.size() != 0) ? m_sizes.size() : (m_connectivity.size() / ShapeType::numberOfNodes());
  }

  /**
   * \brief Execute a function for each zone in the mesh.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename FuncType>
  void for_all_zones(FuncType &&func) const
  {
    const auto nzones = numberOfZones();

    ConnectivityView connectivityView(m_connectivity);
    if constexpr (ShapeType::is_variable_size())
    {
      ConnectivityView sizesView(m_sizes);
      ConnectivityView offsetsView(m_offsets);
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        const ConnectivityView shapeDataView(connectivityView.data() + offsetsView[zoneIndex], sizesView[zoneIndex]);
        const ShapeType shape(shapeDataView);
        func(zoneIndex, shape);
      });
    }
    else
    {
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        const ConnectivityView shapeDataView(connectivityView.data() + ShapeType::zoneOffset(zoneIndex), ShapeType::numberOfNodes());
        const ShapeType shape(shapeDataView);
        func(zoneIndex, shape);
      });
    }
  }

  /**
   * \brief Execute a function for each zone in the mesh.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param selectedIdsView A view containing selected zone ids.
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename ViewType, typename FuncType>
  void for_selected_zones(const ViewType &selectedIdsView, FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();

    ConnectivityView connectivityView(m_connectivity);
    const ViewType localSelectedIdsView(selectedIdsView);
    if constexpr (ShapeType::is_variable_size())
    {
      ConnectivityView sizesView(m_sizes);
      ConnectivityView offsetsView(m_offsets);
      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(auto selectIndex)
      {
        const auto zoneIndex = localSelectedIdsView[selectIndex];
        const ConnectivityView shapeDataView(connectivityView.data() + offsetsView[zoneIndex], sizesView[zoneIndex]);
        const ShapeType shape(shapeDataView);
        func(zoneIndex, shape);
      });
    }
    else
    {
      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(auto selectIndex)
      {
        const auto zoneIndex = localSelectedIdsView[selectIndex];
        const ConnectivityView shapeData(connectivityView.data() + ShapeType::zoneOffset(zoneIndex), ShapeType::numberOfNodes());
        const ShapeType shape(shapeData);
        func(zoneIndex, shape);
      });
    }
  }

private:
  ConnectivityView m_connectivity;
  ConnectivityView m_sizes;
  ConnectivityView m_offsets;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
