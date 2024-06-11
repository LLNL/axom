// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_
#define AXOM_MIR_VIEWS_STRUCTURED_TOPOLOGY_SINGLE_SHAPE_VIEW_HPP_

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
template <typename IndexT, int NDIMS>
class StructuredTopologyView
{
public:
  using IndexType = IndexT;
  using LogicalIndexType = StackArray<IndexType, NDIMS>;

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  AXOM_HOST_VIEW
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

    axom::ArrayView<IndexType> connectivity(m_connectivity);
    axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
    {
      const ShapeType shape(axom::ArrayView<IndexType>(connectivity.data() + ShapeType::zoneOffset(zoneIndex), ShapeType::numberOfNodes()));
      func(zoneIndex, shape);
    });
  }

private:
  StructuredTopologyView() = delete;

  LogicalIndexType m_dimensions;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
