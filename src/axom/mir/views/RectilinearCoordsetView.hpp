// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_UNIFORm_coordinates_VIEW_HPP_
#define AXOM_MIR_UNIFORm_coordinates_VIEW_HPP_

#include "axom/core/StackArray.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/mir/views/StructuredIndexing.hpp"

#include <initializer_list>

namespace axom
{
namespace mir
{
namespace views
{

/// NOTE: The rectilinear coordset views could be combined into a single RectilinearCoordset
///       view that is templated on NDIMS but the resulting SFINAE would just overcomplicate it.

/**
 * \class This class provides a view for Conduit/Blueprint 2D rectilinear coordsets.
 */
template <typename DataType>
class RectilinearCoordsetView2
{
public:
  using LogicalIndex = axom::StackArray<axom::IndexType, 2>;
  using IndexType = axom::IndexType;
  using value_type = DataType;
  using PointType = axom::primal::Point<DataType, 2>;

  /**
   * \brief Constructor
   *
   * \param x The first coordinate component.
   * \param y The second coordinate component.
   */
  AXOM_HOST_DEVICE
  RectilinearCoordsetView2(const axom::ArrayView<DataType> &x,
                           const axom::ArrayView<DataType> &y) : m_coordinates{x, y}, m_indexing(LogicalIndex{{x.size(), y.size()}})
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
    return m_indexing.size();
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The logical index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(LogicalIndex vertex_index) const
  {
    return PointType(std::initializer_list<value_type>{m_coordinates[0][vertex_index[0]],
                                                       m_coordinates[1][vertex_index[1]]});
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(IndexType vertex_index) const
  {
    return getPoint(m_indexing.IndexToLogicalIndex(vertex_index));
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
  operator[](LogicalIndex vertex_index) const
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

private:
  axom::ArrayView<DataType>        m_coordinates[2];
  StructuredIndexing<IndexType, 2> m_indexing;
};

/**
 * \class This class provides a view for Conduit/Blueprint 3D rectilinear coordsets.
 */
template <typename DataType>
class RectilinearCoordsetView3
{
public:
  using LogicalIndex = axom::StackArray<axom::IndexType, 3>;
  using IndexType = axom::IndexType;
  using value_type = DataType;
  using PointType = axom::primal::Point<DataType, 3>;

  /**
   * \brief Constructor
   *
   * \param x The first coordinate component.
   * \param y The second coordinate component.
   * \param z The third coordinate component.
   */
  AXOM_HOST_DEVICE
  RectilinearCoordsetView3(const axom::ArrayView<DataType> &x,
                           const axom::ArrayView<DataType> &y,
                           const axom::ArrayView<DataType> &z) : m_coordinates{x, y, z},  m_indexing(LogicalIndex{{x.size(), y.size(), z.size()}})
  {
  }

  /**
   * \brief Return the number of nodes in the coordset.
   *
   * \return The number of nodes in the coordset.
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
  PointType getPoint(LogicalIndex vertex_index) const
  {
    return PointType(std::initializer_list<value_type>{m_coordinates[0][vertex_index[0]],
                                                       m_coordinates[1][vertex_index[1]],
                                                       m_coordinates[2][vertex_index[2]]});
  }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType getPoint(IndexType vertex_index) const
  {
    return getPoint(m_indexing.IndexToLogicalIndex(vertex_index));
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
  operator[](LogicalIndex vertex_index) const
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

private:
  axom::ArrayView<DataType>        m_coordinates[3];
  StructuredIndexing<IndexType, 3> m_indexing;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
