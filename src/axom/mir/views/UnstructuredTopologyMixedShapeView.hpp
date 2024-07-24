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

/**
 * \brief Given a shape value, we can get the Shape::id() that is used internally.
 */
template <typename IndexT>
class ShapeMap
{
public:
  using IndexType = IndexT;

  /**
   * \brief Constructor
   */
  AXOM_HOST_DEVICE ShapeMap() : m_shape_values(), m_shape_ids()
  {
  }

  /**
   * \brief Constructor
   *
   * \param shape_values A view of sorted values used in the Conduit data.
   * \param shape_ids A view of shape ids that correspond to the values from Conduit.
   */
  AXOM_HOST_DEVICE ShapeMap(const axom::ArrayView<IndexType> &shape_values, const axom::ArrayView<IndexType> &shape_ids) : m_shape_values(shape_values), m_shape_ids(shape_ids)
  {
  }

  /**
   * \brief Return the size of the shape map.
   * \return The number of entries in the shape map.
   */
  AXOM_HOST_DEVICE IndexType size() const
  {
    return m_shape_values.size();
  }

  /**
   * \brief Return whether the shape map is empty.
   * \return True if the map is empty; False otherwise.
   */
  AXOM_HOST_DEVICE bool empty() const
  {
    return m_shape_values.empty();
  }

  /**
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


/**
 * \brief This class provides a view for Conduit/Blueprint mixed shape unstructured grids.
 *
 * \tparam IndexT The index type that will be used for connectivity, etc.
 * \tparam ShapeT The shape type.
 */
template <typename ConnT>
class UnstructuredTopologyMixedShapeView
{
public:
  using ConnectivityType = ConnT;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;
  using ShapeType = VariableShape<IndexType>;

  /**
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
                                     const ConnectivityView &offsets) :
    m_topo(topo), m_connectivity(conn), m_shapes(shapes), m_sizes(sizes), m_offsets(offsets)
  {
    SLIC_ASSERT(m_shapes.size() != 0);
    SLIC_ASSERT(m_sizes.size() != 0);
    SLIC_ASSERT(m_offsets.size() != 0);
    SLIC_ASSERT(m_offsets.size() == m_sizes.size() && m_offsets.size() == m_shapes.size());
  }

  /**
   * \brief Return the dimension of the shape.
   *
   * \return -1 for unknown dimension. We'd have to look at the shapes.
   */
  static constexpr int dimension() { return -1; }

  /**
   * \brief Return the number of zones.
   *
   * \return The number of zones.
   */
  IndexType numberOfZones() const
  {
    return m_sizes.size();
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

    // Build a ShapeMap from the Conduit shape map.
    axom::Array<IndexType> values, ids;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    buildShapeMap(values, ids, allocatorID);
    const ShapeMap<IndexType> shapeMap(values.view(), ids.view());

    const ConnectivityView connectivityView(m_connectivity);
    const ConnectivityView shapes(m_shapes);
    const ConnectivityView sizes(m_sizes);
    const ConnectivityView offsets(m_offsets);
    axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(auto zoneIndex)
    {
      const ConnectivityView shapeData(connectivityView.data() + offsets[zoneIndex], sizes[zoneIndex]);
      const auto shapeID = shapeMap[shapes[zoneIndex]];
      // TODO: SLIC_ASSERT(shapeID > 0);
      const ShapeType shape(shapeID, shapeData);
      func(zoneIndex, shape);
    });
  }

  /**
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
  void for_selected_zones(const ViewType &selectedIdsView, const FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();

    // Build a ShapeMap from the Conduit shape map.
    axom::Array<IndexType> values, ids;
    const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    buildShapeMap(values, ids, allocatorID);
    const ShapeMap<IndexType> shapeMap(values.view(), ids.view());

    const ConnectivityView connectivityView(m_connectivity);
    const ConnectivityView shapes(m_shapes);
    const ConnectivityView sizes(m_sizes);
    const ConnectivityView offsets(m_offsets);
    axom::for_all<ExecSpace>(0, nSelectedZones, AXOM_LAMBDA(int selectIndex)
    {
      const auto zoneIndex = selectedIdsView[selectIndex];
      const ConnectivityView shapeData(connectivityView.data() + offsets[zoneIndex], sizes[zoneIndex]);
      const auto shapeID = shapeMap[shapes[zoneIndex]];
      // TODO: SLIC_ASSERT(shapeID > 0);
      const ShapeType shape(shapeID, shapeData);
      func(zoneIndex, shape);
    });
  }

private:
  /**
   * \brief Populate the shape map values/ids arrays using data in the topology's shape_map.
   *
   * \param[out] values The sorted values used for shapes in the topology.
   * \param[out] ids The Shape ids that correspond to the shape values.
   * \param allocatorID The allocator to use when creating the arrays.
   */
  void buildShapeMap(axom::Array<IndexType> &values, axom::Array<IndexType> &ids, int allocatorID) const
  {
    // Make the map from the Conduit shape_map. Use std::map to sort the key values.
    std::map<IndexType, IndexType> sm;
    const conduit::Node &m_shape_map = m_topo.fetch_existing("elements/shape_map");
    for(conduit::index_t i = 0; i < m_shape_map.number_of_children(); i++)
    {
      const auto value = static_cast<IndexType>(m_shape_map[i].to_int());
      sm[value] = axom::mir::views::shapeNameToID(m_shape_map[i].name());
    }

    // Store the map in 2 vectors so data are contiguous.
    const auto n = sm.size();
    std::vector<IndexType> valuesvec, idsvec;
    valuesvec.reserve(n);
    idsvec.reserve(n);
    for(auto it = sm.begin(); it != sm.end(); it++)
    {
      valuesvec.push_back(it->first);
      idsvec.push_back(it->second);
    }

    // Copy the map values to the device memory.
    values = axom::Array<IndexType>(n, n, allocatorID);
    ids = axom::Array<IndexType>(n, n, allocatorID);
    axom::copy(values.data(), valuesvec.data(), n * sizeof(IndexType));
    axom::copy(ids.data(), idsvec.data(), n * sizeof(IndexType));
  }

  const conduit::Node &      m_topo;
  axom::ArrayView<IndexType> m_connectivity;
  axom::ArrayView<IndexType> m_shapes;
  axom::ArrayView<IndexType> m_sizes;
  axom::ArrayView<IndexType> m_offsets;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
