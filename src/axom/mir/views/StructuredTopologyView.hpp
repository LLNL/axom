// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_
#define AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_

#include "axom/mir/views/Shapes.hpp"

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
  using ShapeType = typename std::conditional<IndexingPolicy::dimension() == 3, HexShape<ConnectivityType>, typename std::conditional<IndexingPolicy::dimension() == 2, QuadShape<ConnectivityType>, LineShape<ConnectivityType>>::type>::type;

  /*!
   * \brief Return the number of dimensions.
   *
   * \return The number of dimensions.
   */
  AXOM_HOST_DEVICE constexpr static int dimension() { return IndexingPolicy::dimension(); }

  /*!
   * \brief Constructor
   */
  StructuredTopologyView() : m_indexing() { }

  /*!
   * \brief Constructor
   *
   * \param indexing The indexing policy for the topology (num zones in each dimension).
   */
  StructuredTopologyView(const IndexingPolicy &indexing) : m_indexing(indexing)
  { }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  IndexType size() const { return m_indexing.size(); }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  IndexType numberOfZones() const { return size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  IndexType connectivitySize() const
  {
    IndexType nodesPerElem = 1;
    for(int d = 0; d < dimension(); d++) nodesPerElem *= 2;
    return numberOfZones() * nodesPerElem;
  }

  /*!
   * \brief Return the mesh logical dimensions.
   *
   * \return The mesh logical dimensions.
   */
  const LogicalIndex &logicalDimensions() const
  {
    return m_indexing.logicalDimensions();
  }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  IndexingPolicy &indexing() { return m_indexing; }

  /*!
   * \brief Return indexing object.
   *
   * \return The indexing object.
   */
  const IndexingPolicy &indexing() const { return m_indexing; }

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
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {
          //using ShapeType = HexShape<IndexType>;

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          const auto kp = nodeIndexing.kStride();
          ConnectivityType data[8];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;
          data[4] = data[0] + kp;
          data[5] = data[1] + kp;
          data[6] = data[2] + kp;
          data[7] = data[3] + kp;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 8));
          func(zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 2)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(auto zoneIndex) {
          //using ShapeType = QuadShape<IndexType>;

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          ConnectivityType data[4];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 4));
          func(zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 1)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nzones,
        AXOM_LAMBDA(axom::IndexType zoneIndex) {
          //using ShapeType = LineShape<IndexType>;

          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          ConnectivityType data[2];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 2));
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

    // Q: Should we make a for_all() that iterates over multiple ranges?
    // Q: Should the logical index be passed to the lambda?

    if constexpr(IndexingPolicy::dimension() == 3)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {
          //using ShapeType = HexShape<IndexType>;

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          const auto kp = nodeIndexing.kStride();
          ConnectivityType data[8];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;
          data[4] = data[0] + kp;
          data[5] = data[1] + kp;
          data[6] = data[2] + kp;
          data[7] = data[3] + kp;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 8));
          func(selectIndex, zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 2)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {
          //using ShapeType = QuadShape<IndexType>;

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          const auto jp = nodeIndexing.jStride();
          ConnectivityType data[4];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;
          data[2] = data[1] + jp;
          data[3] = data[2] - 1;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 4));
          func(selectIndex, zoneIndex, shape);
        });
    }
    else if constexpr(IndexingPolicy::dimension() == 1)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(
        0,
        nSelectedZones,
        AXOM_LAMBDA(axom::IndexType selectIndex) {
          //using ShapeType = LineShape<IndexType>;

          const auto zoneIndex = idsView[selectIndex];
          const auto localLogical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
          ConnectivityType data[2];
          data[0] =
            nodeIndexing.GlobalToGlobal(nodeIndexing.LocalToGlobal(localLogical));
          data[1] = data[0] + 1;

          const ShapeType shape(axom::ArrayView<ConnectivityType>(data, 2));
          func(selectIndex, zoneIndex, shape);
        });
    }
  }

private:
  IndexingPolicy m_indexing;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
