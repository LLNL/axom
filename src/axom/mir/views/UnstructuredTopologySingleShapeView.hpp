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
template <typename IndexT, typename ShapeT>
class UnstructuredTopologySingleShapeView
{
public:
  using IndexType = IndexT;
  using ShapeType = ShapeT;

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   */
  UnstructuredTopologySingleShapeView(const axom::ArrayView<IndexType> &conn) : m_connectivity(conn), m_sizes(), m_offsets()
  {
  }

  /**
   * \brief Constructor
   *
   * \param conn The mesh connectivity.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  UnstructuredTopologySingleShapeView(const axom::ArrayView<IndexType> &conn,
                                      const axom::ArrayView<IndexType> &sizes,
                                      const axom::ArrayView<IndexType> &offsets) : m_connectivity(conn), m_sizes(sizes), m_offsets(offsets)
  {
  }

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

    axom::ArrayView<IndexType> connectivity(m_connectivity);
    if constexpr (ShapeType::is_variable_size())
    {
      assert(m_sizes.size() != 0);
      assert(m_offsets.size() != 0);
      assert(m_offsets.size() == m_sizes.size());

      axom::ArrayView<IndexType> sizes(m_sizes);
      axom::ArrayView<IndexType> offsets(m_offsets);
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        const ShapeType shape(axom::ArrayView<IndexType>(connectivity.data() + offsets[zoneIndex], sizes[zoneIndex]));
        func(zoneIndex, shape);
      });
    }
    else
    {
      axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
      {
        const ShapeType shape(axom::ArrayView<IndexType>(connectivity.data() + ShapeType::zoneOffset(zoneIndex), ShapeType::numberOfNodes()));
        func(zoneIndex, shape);
      });
    }
  }

private:
  axom::ArrayView<IndexType> m_connectivity;
  axom::ArrayView<IndexType> m_sizes;
  axom::ArrayView<IndexType> m_offsets;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
