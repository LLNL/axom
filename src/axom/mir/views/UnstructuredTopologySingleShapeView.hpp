// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_
#define AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir/views/Shapes.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/*!
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
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /*!
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  AXOM_HOST_DEVICE
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn)
    : m_connectivityView(conn)
    , m_sizesView()
    , m_offsetsView()
  { }

  /*!
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  AXOM_HOST_DEVICE
  UnstructuredTopologySingleShapeView(const ConnectivityView &conn,
                                      const ConnectivityView &sizes,
                                      const ConnectivityView &offsets)
    : m_connectivityView(conn)
    , m_sizesView(sizes)
    , m_offsetsView(offsets)
  {
#if defined(AXOM_DEBUG)
  #if defined(AXOM_DEVICE_CODE)
    assert(m_offsetsView.size() == m_sizesView.size());
  #else
    SLIC_ASSERT(m_offsetsView.size() == m_sizesView.size());
  #endif
#endif
  }

  /*!
   * \brief Return the dimension of the shape.
   *
   * \return The dimension of the shape.
   */
  AXOM_HOST_DEVICE static constexpr int dimension()
  {
    return ShapeT::dimension();
  }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE IndexType numberOfZones() const
  {
    return (m_sizesView.size() != 0)
      ? m_sizesView.size()
      : (m_connectivityView.size() / ShapeType::numberOfNodes());
  }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE IndexType connectivitySize() const
  {
    return m_connectivityView.size();
  }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The requested zone.
   *
   * \return The requested zone.
   */
  /// @{
  template <bool _variable_size = ShapeType::is_variable_size()>
  AXOM_HOST_DEVICE typename std::enable_if<_variable_size, ShapeType>::type zone(
    axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG)
  #if defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
  #else
    SLIC_ASSERT(zoneIndex < numberOfZones());
  #endif
#endif

    return ShapeType(
      ConnectivityView(m_connectivityView.data() + m_offsetsView[zoneIndex],
                       m_sizesView[zoneIndex]));
  }

  template <bool _variable_size = ShapeType::is_variable_size()>
  AXOM_HOST_DEVICE typename std::enable_if<!_variable_size, ShapeType>::type zone(
    axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG)
  #if defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
  #else
    SLIC_ASSERT(zoneIndex < numberOfZones());
  #endif
#endif

    ConnectivityView shapeIdsView {};
    if(m_sizesView.empty())
    {
      shapeIdsView = ConnectivityView(
        m_connectivityView.data() + ShapeType::zoneOffset(zoneIndex),
        ShapeType::numberOfNodes());
    }
    else
    {
      shapeIdsView =
        ConnectivityView(m_connectivityView.data() + m_offsetsView[zoneIndex],
                         m_sizesView[zoneIndex]);
    }
    return ShapeType(shapeIdsView);
  }
  /// @}

private:
  ConnectivityView m_connectivityView;
  ConnectivityView m_sizesView;
  ConnectivityView m_offsetsView;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
