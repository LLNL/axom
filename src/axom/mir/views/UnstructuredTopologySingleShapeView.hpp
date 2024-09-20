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
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn)
    : m_connectivityView(conn)
    , m_sizesView()
    , m_offsetsView()
  { }

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn,
                                      const ConnectivityView &sizes,
                                      const ConnectivityView &offsets)
    : m_connectivityView(conn)
    , m_sizesView(sizes)
    , m_offsetsView(offsets)
  {
    SLIC_ASSERT(m_offsetsView.size() == m_sizesView.size());
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
    return (m_sizesView.size() != 0)
      ? m_sizesView.size()
      : (m_connectivityView.size() / ShapeType::numberOfNodes());
  }

  /**
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  IndexType connectivitySize() const { return m_connectivityView.size(); }

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

    ConnectivityView connectivityView(m_connectivityView);
    if constexpr(ShapeType::is_variable_size())
    {
      ConnectivityView sizesView(m_sizesView);
      ConnectivityView offsetsView(m_offsetsView);
      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(auto zoneIndex) {
          const ConnectivityView shapeIdsView(
            connectivityView.data() + offsetsView[zoneIndex],
            sizesView[zoneIndex]);
          const ShapeType shape(shapeIdsView);
          func(zoneIndex, shape);
        });
    }
    else
    {
      ConnectivityView sizesView(m_sizesView);
      ConnectivityView offsetsView(m_offsetsView);
      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(auto zoneIndex) {
          ConnectivityView shapeIdsView {};
          if(sizesView.empty())
          {
            shapeIdsView = ConnectivityView(
              connectivityView.data() + ShapeType::zoneOffset(zoneIndex),
              ShapeType::numberOfNodes());
          }
          else
          {
            shapeIdsView =
              ConnectivityView(connectivityView.data() + offsetsView[zoneIndex],
                               sizesView[zoneIndex]);
          }
          const ShapeType shape(shapeIdsView);
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

    ConnectivityView connectivityView(m_connectivityView);
    const ViewType localSelectedIdsView(selectedIdsView);
    if constexpr(ShapeType::is_variable_size())
    {
      ConnectivityView sizesView(m_sizesView);
      ConnectivityView offsetsView(m_offsetsView);
      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(auto selectIndex) {
          const auto zoneIndex = localSelectedIdsView[selectIndex];
          const ConnectivityView shapeIdsView(
            connectivityView.data() + offsetsView[zoneIndex],
            sizesView[zoneIndex]);
          const ShapeType shape(shapeIdsView);
          func(selectIndex, zoneIndex, shape);
        });
    }
    else
    {
      ConnectivityView sizesView(m_sizesView);
      ConnectivityView offsetsView(m_offsetsView);
      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(auto selectIndex) {
          const auto zoneIndex = localSelectedIdsView[selectIndex];
          ConnectivityView shapeIdsView {};
          if(sizesView.empty())
          {
            shapeIdsView = ConnectivityView(
              connectivityView.data() + ShapeType::zoneOffset(zoneIndex),
              ShapeType::numberOfNodes());
          }
          else
          {
            shapeIdsView =
              ConnectivityView(connectivityView.data() + offsetsView[zoneIndex],
                               sizesView[zoneIndex]);
          }
          const ShapeType shape(shapeIdsView);
          func(selectIndex, zoneIndex, shape);
        });
    }
  }

private:
  ConnectivityView m_connectivityView;
  ConnectivityView m_sizesView;
  ConnectivityView m_offsetsView;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
