// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_EXPLICIT_COORDSET_VIEW_HPP_
#define AXOM_MIR_EXPLICIT_COORDSET_VIEW_HPP_

#include "axom/core/ArrayView.hpp"
#include "axom/primal/geometry/Point.hpp"

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \class This class provides a view for Conduit/Blueprint explicit coordsets.
 */
template <typename DataType>
class ExplicitCoordsetView2
{
public:
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
  ExplicitCoordsetView2(const axom::ArrayView<DataType> &x,
                        const axom::ArrayView<DataType> &y) : m_coordinates{x,y}
  {
    assert(x.size() == y.size());
  }

  /**
   * \brief Return the number of points in the coordset.
   *
   * \return The number of points in the coordset.
   */
  AXOM_HOST_DEVICE
  IndexType size() const { return m_coordinates[0].size(); }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType
  getPoint(IndexType vertex_index) const
  {
    assert(vertex_index < size());
    return PointType(m_coordinates[0][vertex_index],
                     m_coordinates[1][vertex_index]);
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
    return getPoint(vertex_index);
  }

private:
  axom::ArrayView<DataType> m_coordinates[2];
};

/**
 * \class This class provides a view for Conduit/Blueprint explicit coordsets.
 */
template <typename DataType>
class ExplicitCoordsetView3
{
public:
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
  ExplicitCoordsetView3(const axom::ArrayView<DataType> &x,
                        const axom::ArrayView<DataType> &y,
                        const axom::ArrayView<DataType> &z) : m_coordinates{x,y,z}
  {
    assert(x.size() == y.size() && x.size() == z.size());
  }

  /**
   * \brief Return the number of points in the coordset.
   *
   * \return The number of points in the coordset.
   */
  AXOM_HOST_DEVICE
  IndexType size() const { return m_coordinates[0].size(); }

  /**
   * \brief Return the requested point from the coordset.
   *
   * \param vertex_index The index of the point to return.
   *
   * \return A point that corresponds to \a vertex_index.
   */
  AXOM_HOST_DEVICE
  PointType
  getPoint(IndexType vertex_index) const
  {
    return PointType(m_coordinates[0][vertex_index],
                     m_coordinates[1][vertex_index],
                     m_coordinates[2][vertex_index]);
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
    return getPoint(vertex_index);
  }

private:
  axom::ArrayView<DataType> m_coordinates[3];
};


} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
