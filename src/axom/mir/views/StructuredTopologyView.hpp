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
  using LogicalIndexType = StackArray<IndexType, NDIMS>;
  constexpr static int Dimensions = NDIMS;

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  AXOM_HOST_DEVICE
  StructuredTopologyView(const LogicalIndexType &dims) : m_dimensions(dims)
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
     IndexType sz = 1;
     for(int i = 0; i < NDIMS; i++)
       sz *= m_dimensions[i];
     return sz;
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
  const LogicalIndexType &dimensions() const { return m_dimensions; }

  /**
   * \brief Execute a function for each zone in the mesh.
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

    if constexpr (NDIMS == 2)
    {
      // Q: Should we make a for_all() that iterates over multiple ranges?
      // Q: Should the logical index be passed to the lambda?

      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        using ShapeType = QuadShape<IndexType>;
        IndexType data[4];
        const ShapeType shape(axom::ArrayView<IndexType>(data, 4));
        func(zoneIndex, shape);
      });
    }
    if constexpr (NDIMS == 3)
    {
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        using ShapeType = HexShape<IndexType>;
        IndexType data[8];
        const ShapeType shape(axom::ArrayView<IndexType>(data, 8));
        func(zoneIndex, shape);
      });
    }
  }

private:
  StructuredTopologyView() = delete;

  LogicalIndexType m_dimensions;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
