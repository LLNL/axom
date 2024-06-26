// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UNIFORM_COORDSET_VIEW_HPP_
#define AXOM_MIR_UNIFORM_COORDSET_VIEW_HPP_

#include "axom/core/StackArray.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/primal/geometry/Point.hpp"

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
  using LogicalIndexType = axom::StackArray<axom::IndexType, NDIMS>;
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
                      const ExtentsType spacing) : m_dimensions{dims}, m_origin{origin}, m_spacing{spacing}
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
     IndexType sz = 1;
     for(int i = 0; i < NDIMS; i++)
       sz *= m_dimensions[i];
     return sz;
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The logical index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(LogicalIndexType vertex_index) const
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
  operator[](LogicalIndexType vertex_index) const
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
    return getPoint(IndexToLogicalIndex(vertex_index));
  }

  
private:
  /**
   * \brief Turn an index into a logical IJ index.
   *
   * \param index The index to convert.
   *
   * \return The logical index that corresponds to the \a index.
   */
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 2, LogicalIndexType>::type
  IndexToLogicalIndex(IndexType index) const
  {
    LogicalIndexType logical;
    const auto nx = m_dimensions[0];
    assert(index >= 0);
    logical[0] = index % nx;
    logical[1] = index / nx;
    return logical;
  }

  /**
   * \brief Turn an index into a logical IJK index.
   *
   * \param index The index to convert.
   *
   * \return The logical index that corresponds to the \a index.
   */
  template <size_t _ndims = NDIMS>
  AXOM_HOST_DEVICE
  typename std::enable_if<_ndims == 3, LogicalIndexType>::type
  IndexToLogicalIndex(IndexType index) const
  {
    LogicalIndexType logical;
    const auto nx = m_dimensions[0];
    const auto nxy = nx * m_dimensions[1];
    logical[0] = index % nx;
    logical[1] = (index % nxy) / nx;
    logical[2] = index / nxy;
    return logical;
  }

  LogicalIndexType m_dimensions;
  ExtentsType      m_origin;
  ExtentsType      m_spacing;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
