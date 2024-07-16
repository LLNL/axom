// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UNIFORM_COORDSET_VIEW_HPP_
#define AXOM_MIR_UNIFORM_COORDSET_VIEW_HPP_

#include "axom/core/StackArray.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \class This class provides a view for Conduit/Blueprint uniform coordsets.
 *
 * \tparam DataType The underlying type used for the coordinates.
 * \tparam NDIMS The number of dimensions in each point.
 *
 */
template <typename DataType, int NDIMS = 2>
class UniformCoordsetView
{
public:
  using LogicalIndex = axom::StackArray<axom::IndexType, NDIMS>;
  using ExtentsType = axom::StackArray<double, NDIMS>;;
  using IndexType = axom::IndexType;
  using value_type = DataType;
  using PointType = axom::primal::Point<DataType, NDIMS>;

  /**
   * \brief Constructor
   *
   * \param dims    The logical dimensions of the coordset.
   * \param origin  The origin of the coordset's coordinate system.
   * \param spacing The spacing inbetween points.
   */
  AXOM_HOST_DEVICE
  UniformCoordsetView(const LogicalIndex &dims,
                      const ExtentsType &origin,
                      const ExtentsType &spacing) : m_indexing(dims), m_origin(origin), m_spacing(spacing)
  {
  }

  /**
   * \brief Return the number of points in the coordset.
   *
   * \return The number of points in the coordset.
   */
  /// @{
  AXOM_HOST_DEVICE
  IndexType size() const
  {
    return m_indexing.size();
  }

  AXOM_HOST_DEVICE
  IndexType numberOfNodes() const
  {
    return m_indexing.size();
  }
  /// @}

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The logical index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(const LogicalIndex &vertex_index) const
  {
    PointType pt;
    for(int i = 0; i < NDIMS; i++)
      pt[i] = m_origin[i] + vertex_index[i] * m_spacing[i];
    return pt;
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The logical index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType
  operator[](const LogicalIndex &vertex_index) const
  {
    return getPoint(vertex_index);
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType
  operator[](IndexType vertex_index) const
  {
    return getPoint(m_indexing.IndexToLogicalIndex(vertex_index));
  }

  StructuredIndexing<IndexType, NDIMS> m_indexing;
  ExtentsType                          m_origin;
  ExtentsType                          m_spacing;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
