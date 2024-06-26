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
  using LogicalIndexType = typename StructuredIndexing<axom::IndexType, NDIMS>::LogicalIndex;
  using ExtentsType = axom::StackArray<DataType, NDIMS>;
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
  UniformCoordsetView(const LogicalIndexType dims,
                      const ExtentsType origin,
                      const ExtentsType spacing) : m_shape{dims}, m_origin{origin}, m_spacing{spacing}
  {
  }

  /**
   * \brief Return the number of points in the coordset.
   *
   * \return The number of points in the coordset.
   */
  AXOM_HOST_DEVICE
  IndexType size() const
  {
    return m_shape.size();
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The logical index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(const LogicalIndexType &vertex_index) const
  {
    PointType pt;
    for(int i = 0; i < NDIMS; i++)
      pt.x[i] = m_origin[i] + vertex_index[i] * m_spacing[i];
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
  operator[](const LogicalIndexType &vertex_index) const
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
    return getPoint(m_shape.IndexToLogicalIndex(vertex_index));
  }

  StructuredIndexing<IndexType, NDIMS> m_shape;
  ExtentsType                          m_origin;
  ExtentsType                          m_spacing;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
