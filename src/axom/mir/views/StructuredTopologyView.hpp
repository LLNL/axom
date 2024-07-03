// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_
#define AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_VIEW_HPP_

#include "axom/mir/views/Shapes.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"

namespace axom
{
namespace mir
{
namespace views
{

// NOTE: we could subclass this one and make a strided structured view that lets one define zones where the zones do not span all of the nodes.

/**
 * \brief This class provides a view for Conduit/Blueprint single shape unstructured grids.
 *
 * \tparam IndexT The index type that will be used for connectivity, etc.
 * \tparam ShapeT The shape type.
 */
template <typename IndexT, int NDIMS>
class StructuredTopologyView
{
public:
  using IndexType = IndexT;
  using LogicalIndexType = axom::StackArray<IndexT, NDIMS>; //typename StructuredIndexing<IndexType, NDIMS>::LogicalIndex;

  constexpr static int dimension() { return NDIMS; }

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  AXOM_HOST_DEVICE
  StructuredTopologyView(const LogicalIndexType &dims) : m_shape(dims)
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
     return m_shape.size();
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
  const LogicalIndexType &logicalDimensions() const { return m_shape; }

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

    if constexpr (NDIMS == 3)
    {
      const StructuredIndexing<IndexType, NDIMS> zoneShape{m_shape},
                                                 nodeShape{{m_shape.m_dimensions[0] + 1,
                                                            m_shape.m_dimensions[1] + 1,
                                                            m_shape.m_dimensions[2] + 1}};
      const auto jp = nodeShape.jStride();
      const auto kp = nodeShape.kStride();
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        using ShapeType = HexShape<IndexType>;

        const auto logical = zoneShape.IndexToLogicalIndex(zoneIndex);
        IndexType data[8];
        data[0] = nodeShape.LogicalIndexToIndex(logical);
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
    else if constexpr (NDIMS == 2)
    {
      // Q: Should we make a for_all() that iterates over multiple ranges?
      // Q: Should the logical index be passed to the lambda?

      const StructuredIndexing<IndexType, NDIMS> zoneShape{m_shape},
                                                 nodeShape{{m_shape.m_dimensions[0] + 1,
                                                            m_shape.m_dimensions[1] + 1}};
      const auto jp = nodeShape.jStride();

      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        using ShapeType = QuadShape<IndexType>;

        const auto logical = zoneShape.IndexToLogicalIndex(zoneIndex);
        IndexType data[4];
        data[0] = nodeShape.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 4));
        func(zoneIndex, shape);
      });
    }
    // TODO: NDIMS == 1
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

    if constexpr (NDIMS == 3)
    {
      const StructuredIndexing<IndexType, NDIMS> zoneShape{m_shape},
                                                 nodeShape{{m_shape.m_dimensions[0] + 1,
                                                            m_shape.m_dimensions[1] + 1,
                                                            m_shape.m_dimensions[2] + 1}};
      const auto jp = nodeShape.jStride();
      const auto kp = nodeShape.kStride();
      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(int selectIndex)
      {
        using ShapeType = HexShape<IndexType>;

        const auto zoneIndex = idsView[selectIndex];
        const auto logical = zoneShape.IndexToLogicalIndex(zoneIndex);
        IndexType data[8];
        data[0] = nodeShape.LogicalIndexToIndex(logical);
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
    else if constexpr (NDIMS == 2)
    {
      // Q: Should we make a for_all() that iterates over multiple ranges?
      // Q: Should the logical index be passed to the lambda?

      const StructuredIndexing<IndexType, NDIMS> zoneShape{m_shape},
                                                 nodeShape{{m_shape.m_dimensions[0] + 1,
                                                            m_shape.m_dimensions[1] + 1}};
      const auto jp = nodeShape.jStride();

      axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(int selectIndex)
      {
        using ShapeType = QuadShape<IndexType>;

        const auto zoneIndex = idsView[selectIndex];
        const auto logical = zoneShape.IndexToLogicalIndex(zoneIndex);
        IndexType data[4];
        data[0] = nodeShape.LogicalIndexToIndex(logical);
        data[1] = data[0] + 1;
        data[2] = data[1] + jp;
        data[3] = data[2] - 1;

        const ShapeType shape(axom::ArrayView<IndexType>(data, 4));
        func(zoneIndex, shape);
      });
    }
    // TODO: NDIMS == 1
  }

private:
  StructuredIndexing<IndexType, NDIMS> m_shape;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
