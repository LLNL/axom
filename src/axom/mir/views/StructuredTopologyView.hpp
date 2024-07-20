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

/**
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

  /**
   * \brief Return the number of dimensions.
   *
   * \return The number of dimensions.
   */
  AXOM_HOST_DEVICE
  constexpr static int dimension() { return IndexingPolicy::dimensions(); }

  /**
   * \brief Constructor
   *
   * \param indexing The indexing policy for the topology (num zones in each dimension).
   */
  AXOM_HOST_DEVICE
  StructuredTopologyView(const IndexingPolicy &indexing) : m_indexing(indexing)
  {
  }

  /**
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE
  IndexType size() const
  {
     return m_indexing.size();
  }

  /**
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE IndexType numberOfZones() const { return size(); }

  /**
   * \brief Return the mesh logical dimensions.
   *
   * \return The mesh logical dimensions.
   */
  const LogicalIndex &logicalDimensions() const { return m_indexing.logicalDimensions(); }

  /**
   * \brief Execute a function for each zone in the mesh using axom::for_all.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename FuncType>
  AXOM_HOST
  void for_all_zones(FuncType &&func) const
  {
    const auto nzones = numberOfZones();

    // Q: Should we make a for_all() that iterates over multiple ranges?
    // Q: Should the logical index be passed to the lambda?

    if constexpr (IndexingPolicy::dimensions() == 3)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        using ShapeType = HexShape<IndexType>;

        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        const auto jp = nodeIndexing.jStride();
        const auto kp = nodeIndexing.kStride();
        IndexType data[8];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;
        data[4] = data[0] + kp;
        data[5] = data[1] + kp;
        data[6] = data[2] + kp;
        data[7] = data[3] + kp;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 8));
        func(zoneIndex, shape);
      });
    }
    else if constexpr (IndexingPolicy::dimensions() == 2)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        using ShapeType = QuadShape<IndexType>;

        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        const auto jp = nodeIndexing.jStride();
        IndexType data[4];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 4));
        func(zoneIndex, shape);
      });
    }
    else if constexpr (IndexingPolicy::dimensions() == 1)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
      {
        using ShapeType = LineShape<IndexType>;

        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        IndexType data[2];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 2));
        func(zoneIndex, shape);
      });
    }
  }

  /**
   * \brief Execute a function for each zone in the mesh using axom::for_all.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam ViewType  A concrete ArrayView that contains selected zone ids.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename ViewType, typename FuncType>
  void for_selected_zones(const ViewType &selectedIdsView, const FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();
    ViewType idsView(selectedIdsView);

    // Q: Should we make a for_all() that iterates over multiple ranges?
    // Q: Should the logical index be passed to the lambda?

    if constexpr (IndexingPolicy::dimensions() == 3)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(auto selectIndex)
      {
        using ShapeType = HexShape<IndexType>;

        const auto zoneIndex = idsView[selectIndex];
        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        const auto jp = nodeIndexing.jStride();
        const auto kp = nodeIndexing.kStride();
        IndexType data[8];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;
        data[4] = data[0] + kp;
        data[5] = data[1] + kp;
        data[6] = data[2] + kp;
        data[7] = data[3] + kp;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 8));
        func(zoneIndex, shape);
      });
    }
    else if constexpr (IndexingPolicy::dimensions() == 2)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(auto selectIndex)
      {
        using ShapeType = QuadShape<IndexType>;

        const auto zoneIndex = idsView[selectIndex];
        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        const auto jp = nodeIndexing.jStride();
        IndexType data[4];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 4));
        func(zoneIndex, shape);
      });
    }
    else if constexpr (IndexingPolicy::dimensions() == 1)
    {
      const IndexingPolicy zoneIndexing = m_indexing;
      const IndexingPolicy nodeIndexing = m_indexing.expand();

      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(auto selectIndex)
      {
        using ShapeType = LineShape<IndexType>;

        const auto zoneIndex = idsView[selectIndex];
        const auto logical = zoneIndexing.IndexToLogicalIndex(zoneIndex);
        IndexType data[2];
        data[0] = nodeIndexing.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 2));
        func(zoneIndex, shape);
      });
    }
  }

private:
  IndexingPolicy m_indexing;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
