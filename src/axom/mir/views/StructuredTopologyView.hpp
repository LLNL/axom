// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_
#define AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/Shapes.hpp"

#include <type_traits>

namespace axom
{
namespace mir
{
namespace views
{
/*!
 * \brief This class provides a view for Conduit/Blueprint structured grid types.
 *
 * \tparam IndexPolicy The policy for making/using indices.
 */
template <typename IndexPolicy>
class StructuredTopologyView
{
public:
  using IndexingPolicy = IndexPolicy;
  using IndexType = typename IndexingPolicy::IndexType;
  using LogicalIndex = typename IndexingPolicy::LogicalIndex;
  using ConnectivityType = IndexType;
  using Shape1D = LineShape<axom::StackArray<ConnectivityType, 2>>;
  using Shape2D = QuadShape<axom::StackArray<ConnectivityType, 4>>;
  using Shape3D = HexShape<axom::StackArray<ConnectivityType, 8>>;
  using ShapeType = typename std::conditional<IndexingPolicy::dimension() == 3, Shape3D, typename std::conditional<IndexingPolicy::dimension() == 2, Shape2D, Shape1D>::type>::type;

  /*!
   * \brief Return the number of dimensions.
   *
   * \return The number of dimensions.
   */
  AXOM_HOST_DEVICE constexpr static int dimension() { return IndexingPolicy::dimension(); }

  /*!
   * \brief Constructor
   */
  AXOM_HOST_DEVICE StructuredTopologyView() : m_zoneIndexing(), m_nodeIndexing() { }

  /*!
   * \brief Constructor
   *
   * \param indexing The indexing policy for the topology (num zones in each dimension).
   */
  AXOM_HOST_DEVICE StructuredTopologyView(const IndexingPolicy &indexing) : m_zoneIndexing(indexing), m_nodeIndexing(indexing.expand())
  {
  }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE IndexType size() const { return m_zoneIndexing.size(); }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE IndexType numberOfZones() const { return size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE IndexType connectivitySize() const
  {
    IndexType nodesPerElem = 1;
    for(int d = 0; d < dimension(); d++)
    {
      nodesPerElem *= 2;
    }
    return numberOfZones() * nodesPerElem;
  }

  /*!
   * \brief Return the mesh logical dimensions.
   *
   * \return The mesh logical dimensions.
   */
  AXOM_HOST_DEVICE const LogicalIndex &logicalDimensions() const
  {
    return m_zoneIndexing.logicalDimensions();
  }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  AXOM_HOST_DEVICE IndexingPolicy &indexing() { return m_zoneIndexing; }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  AXOM_HOST_DEVICE const IndexingPolicy &indexing() const { return m_zoneIndexing; }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 3D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 3, Shape3D>::type
  zone(axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
#endif
    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);
    const auto jp = m_nodeIndexing.jStride();
    const auto kp = m_nodeIndexing.kStride();

    Shape3D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;
    data[2] = data[1] + jp;
    data[3] = data[2] - 1;
    data[4] = data[0] + kp;
    data[5] = data[1] + kp;
    data[6] = data[2] + kp;
    data[7] = data[3] + kp;

    return shape;
  }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 2D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 2, Shape2D>::type
  zone(axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
#endif
    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);
    const auto jp = m_nodeIndexing.jStride();

    Shape2D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;
    data[2] = data[1] + jp;
    data[3] = data[2] - 1;

    return shape;
  }

  /*!
   * \brief Return a zone.
   *
   * \param index The index of the zone to return.
   *
   * \return The requested zone.
   *
   * \note 1D implementation.
   */
  template <int _ndims = IndexingPolicy::dimension()>
  AXOM_HOST_DEVICE typename std::enable_if<_ndims == 1, Shape1D>::type
  zone(axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
#endif
    const auto localLogical = m_zoneIndexing.IndexToLogicalIndex(zoneIndex);

    Shape1D shape;
    auto &data = shape.getIdsStorage();
    data[0] = m_nodeIndexing.GlobalToGlobal(m_nodeIndexing.LocalToGlobal(localLogical));
    data[1] = data[0] + 1;

    return shape;
  }

  /*!
   * \brief Execute a function for each zone in the mesh using axom::for_all.
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

    // Q: Should we make a for_all() that iterates over multiple ranges?
    // Q: Should the logical index be passed to the lambda?

    if constexpr(IndexingPolicy::dimension() == 3)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          const auto kp = nodeIndexing.kStride();
          Shape3D shape;
          auto &data = shape.getIdsStorage();
          data[0] = nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;
          data[4] = data[0] + kp;
          data[5] = data[1] + kp;
          data[6] = data[2] + kp;
          data[7] = data[3] + kp;

          func(zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 2)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          Shape2D shape;
          auto &data = shape.getIdsStorage();
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;

          func(zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 1)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          Shape1D shape;
          auto &data = shape.getIdsStorage();
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;

          func(zoneIndex, shape);
        });
    }
  }

  /*!
   * \brief Execute a function for each zone in the mesh using axom::for_all.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam ViewType  A concrete ArrayView that contains selected zone ids.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename ViewType, typename FuncType>
  void for_selected_zones(const ViewType &selectedIdsView,
                          const FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();
    ViewType idsView(selectedIdsView);

    if constexpr(IndexingPolicy::dimension() == 3)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          const auto kp = nodeIndexing.kStride();
          Shape3D shape;
          auto &data = shape.getIdsStorage();
          data[0] = nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;
          data[4] = data[0] + kp;
          data[5] = data[1] + kp;
          data[6] = data[2] + kp;
          data[7] = data[3] + kp;

          func(selectIndex, zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 2)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          Shape2D shape;
          auto &data = shape.getIdsStorage();

          data[0] = nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;

          func(selectIndex, zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 1)
    {
      const IndexingPolicy zoneIndexing = m_zoneIndexing;
      const IndexingPolicy nodeIndexing = m_zoneIndexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          Shape1D shape;
          auto &data = shape.getIdsStorage();

          data[0] = nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;

          func(selectIndex, zoneIndex, shape);
        });
    }
  }

private:
  IndexingPolicy m_zoneIndexing;
  IndexingPolicy m_nodeIndexing;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
