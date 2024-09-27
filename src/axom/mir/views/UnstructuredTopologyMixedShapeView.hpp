// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_MIXED_SHAPE_VIEW_HPP_
#define AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_MIXED_SHAPE_VIEW_HPP_

#include "axom/mir/views/Shapes.hpp"
#include "axom/mir/utilities.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/*!
 * \brief Given a shape value, we can get the Shape::id() that is used internally.
 *
 * \note If the view was to renumber the shapes array to use the Shape::id() values
 *       then operator[] could return its input value and skip bsearch.
 */
class ShapeMap
{
public:
  /*!
   * \brief Constructor
   */
  AXOM_HOST_DEVICE ShapeMap() : m_shape_values(), m_shape_ids() { }

  /*!
   * \brief Constructor
   *
   * \param shape_values A view of sorted values used in the Conduit data.
   * \param shape_ids A view of shape ids that correspond to the values from Conduit.
   */
  AXOM_HOST_DEVICE ShapeMap(const axom::ArrayView<IndexType> &shape_values,
                            const axom::ArrayView<IndexType> &shape_ids)
    : m_shape_values(shape_values)
    , m_shape_ids(shape_ids)
  { }

  /*!
   * \brief Return the size of the shape map.
   * \return The number of entries in the shape map.
   */
  AXOM_HOST_DEVICE IndexType size() const { return m_shape_values.size(); }

  /*!
   * \brief Return whether the shape map is empty.
   * \return True if the map is empty; False otherwise.
   */
  AXOM_HOST_DEVICE bool empty() const { return m_shape_values.empty(); }

  /*!
   * \brief Given a shape value (as in the Conduit shapes array), return the shape id.
   *
   * \param value A value from the shapes array that we want to map to a shape id.
   *
   * \return A shape id.
   */
  AXOM_HOST_DEVICE IndexType operator[](IndexType value) const
  {
    const auto index = axom::mir::utilities::bsearch(value, m_shape_values);
    return (index >= 0) ? m_shape_ids[index] : 0;
  }

private:
  axom::ArrayView<IndexType> m_shape_values;
  axom::ArrayView<IndexType> m_shape_ids;
};

/*!
 * \brief Populate the shape map values/ids arrays using data in the topology's shape_map.
 *
 * \param n_topo The topology that contains the shape map.
 * \param[out] values The sorted values used for shapes in the topology.
 * \param[out] ids The Shape ids that correspond to the shape values.
 * \param allocatorID The allocator to use when creating the arrays.
 */
ShapeMap buildShapeMap(const conduit::Node &n_topo,
                       axom::Array<IndexType> &values,
                       axom::Array<IndexType> &ids,
                       int allocatorID);
/*!
 * \brief This class provides a view for Conduit/Blueprint mixed shape unstructured grids.
 *
 * \tparam IndexT The index type that will be used for connectivity, etc.
 * \tparam ShapeT The shape type.
 *
 * \note This view does not support topologies that also contain polyhedral elements.
 */
template <typename ConnT>
class UnstructuredTopologyMixedShapeView
{
public:
  using ConnectivityType = ConnT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using ShapeType = VariableShape<ConnectivityType>;

  /*!
   * \brief Constructor
   *
   * \param topo A reference to the topology.
   * \param conn The mesh connectivity.
   * \param shapes The shape in each zone.
   * \param sizes The number of nodes in each zone.
   * \param offsets The offset to each zone in the connectivity.
   */
  UnstructuredTopologyMixedShapeView(const conduit::Node &topo,
                                     const ConnectivityView &conn,
                                     const ConnectivityView &shapes,
                                     const ConnectivityView &sizes,
                                     const ConnectivityView &offsets,
                                     const ShapeMap &shapemap)
    : m_topo(topo)
    , m_connectivity(conn)
    , m_shapes(shapes)
    , m_sizes(sizes)
    , m_offsets(offsets)
    , m_shapeMap(shapemap)
  {
#if defined(AXOM_DEBUG)
#if defined(AXOM_DEVICE_CODE)
    assert(m_shapes.size() != 0);
    assert(m_sizes.size() != 0);
    assert(m_offsets.size() != 0);
    assert(m_offsets.size() == m_sizes.size() &&
           m_offsets.size() == m_shapes.size());
#else
    SLIC_ASSERT(m_shapes.size() != 0);
    SLIC_ASSERT(m_sizes.size() != 0);
    SLIC_ASSERT(m_offsets.size() != 0);
    SLIC_ASSERT(m_offsets.size() == m_sizes.size() &&
                m_offsets.size() == m_shapes.size());
#endif
#endif
  }

  /*!
   * \brief Return the dimension of the shape.
   *
   * \return -1 for unknown dimension. We'd have to look at the shapes.
   */
  AXOM_HOST_DEVICE static constexpr int dimension() { return -1; }

  /*!
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  AXOM_HOST_DEVICE IndexType numberOfZones() const { return m_sizes.size(); }

  /*!
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  AXOM_HOST_DEVICE IndexType connectivitySize() const { return m_connectivity.size(); }

  /*!
   * \brief Return a zone.
   *
   * \param zoneIndex The index of the zone to return.
   *
   * \return The requested zone.
   */
  ShapeType zone(axom::IndexType zoneIndex) const
  {
#if defined(AXOM_DEBUG)
#if defined(AXOM_DEVICE_CODE)
    assert(zoneIndex < numberOfZones());
#else
    SLIC_ASSERT(zoneIndex < numberOfZones());
#endif
#endif
    const ConnectivityView shapeData(
          m_connectivity.data() + m_offsets[zoneIndex],
          m_sizes[zoneIndex]);
    const auto shapeID = m_shapeMap[m_shapes[zoneIndex]];
#if defined(AXOM_DEBUG)
#if defined(AXOM_DEVICE_CODE)
    assert(shapeID > 0);
#else
    SLIC_ASSERT(shapeID > 0);
#endif
#endif

    return ShapeType(shapeID, shapeData);
  }

  /*!
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

    // Build a ShapeMap from the Conduit shape map.
    axom::Array<IndexType> values, ids;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    buildShapeMap(m_topo, values, ids, allocatorID);
    const ShapeMap shapeMap(values.view(), ids.view());

    const ConnectivityView connectivityView(m_connectivity);
    const ConnectivityView shapes(m_shapes);
    const ConnectivityView sizes(m_sizes);
    const ConnectivityView offsets(m_offsets);
    axom::for_all<ExecSpace>(
      0,
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const ConnectivityView shapeData(
          connectivityView.data() + offsets[zoneIndex],
          sizes[zoneIndex]);
        const auto shapeID = shapeMap[shapes[zoneIndex]];
#if defined(AXOM_DEBUG)
#if defined(AXOM_DEVICE_CODE)
        assert(shapeID > 0);
#else
        SLIC_ASSERT(shapeID > 0);
#endif
#endif
        const ShapeType shape(shapeID, shapeData);
        func(zoneIndex, shape);
      });
  }

  /*!
   * \brief Execute a function for each zone in the mesh.
   *
   * \tparam ExecSpace The execution space for the function body.
   * \tparam ViewType  A view type that contains zone indices.
   * \tparam FuncType  The type for the function/lambda to execute. It will accept a zone index and shape.
   *
   * \param selectedIdsView A view that contains a list of zones to operate on.
   * \param func The function/lambda that will be executed for each zone in the mesh.
   */
  template <typename ExecSpace, typename ViewType, typename FuncType>
  void for_selected_zones(const ViewType &selectedIdsView, FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();

    // Build a ShapeMap from the Conduit shape map.
    axom::Array<IndexType> values, ids;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    buildShapeMap(m_topo, values, ids, allocatorID);
    const ShapeMap shapeMap(values.view(), ids.view());

    // Make views that can be captured.
    const ConnectivityView connectivityView(m_connectivity);
    const ConnectivityView shapes(m_shapes);
    const ConnectivityView sizes(m_sizes);
    const ConnectivityView offsets(m_offsets);
    const ViewType deviceSelectedIdsView(selectedIdsView);
    axom::for_all<ExecSpace>(
      0,
      nSelectedZones,
      AXOM_LAMBDA(axom::IndexType selectIndex) {
        const auto zoneIndex = deviceSelectedIdsView[selectIndex];
        const ConnectivityView shapeData(
          connectivityView.data() + offsets[zoneIndex],
          sizes[zoneIndex]);
        const auto shapeID = shapeMap[shapes[zoneIndex]];
#if defined(AXOM_DEBUG)
#if defined(AXOM_DEVICE_CODE)
        assert(shapeID > 0);
#else
        SLIC_ASSERT(shapeID > 0);
#endif
#endif
        const ShapeType shape(shapeID, shapeData);
        func(selectIndex, zoneIndex, shape);
      });
  }

private:

  const conduit::Node &m_topo;
  ConnectivityView m_connectivity;
  ConnectivityView m_shapes;
  ConnectivityView m_sizes;
  ConnectivityView m_offsets;
  ShapeMap         m_shapeMap;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
