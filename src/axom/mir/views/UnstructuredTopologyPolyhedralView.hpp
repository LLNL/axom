// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_POLYHEDRAL_VIEW_HPP_
#define AXOM_MIR_VIEWS_UNSTRUCTURED_TOPOLOGY_POLYHEDRAL_VIEW_HPP_

#include "axom/mir/views/Shapes.hpp"

namespace axom
{
namespace mir
{
namespace views
{
/**
 * \brief This class implements a view for Blueprint polyhedral topologies.
 */
template <typename ConnType>
class UnstructuredTopologyPolyhedralView
{
public:
  using ConnectivityType = ConnType;
  using ConnectivityView = axom::ArrayView<ConnectivityType>;

  /**
   * \brief This struct contains views that hold polyhedral connectivity.
   */
  struct PolyhedronData
  {
    /// Constructor
    AXOM_HOST_DEVICE
    PolyhedronData(const ConnectivityView &subelement_conn,
                   const ConnectivityView &subelement_sizes,
                   const ConnectivityView &subelement_offsets,
                   const ConnectivityView &element_conn,
                   const ConnectivityView &element_sizes,
                   const ConnectivityView &element_offsets)
      : m_subelement_conn(subelement_conn)
      , m_subelement_sizes(subelement_sizes)
      , m_subelement_offsets(subelement_offsets)
      , m_element_conn(element_conn)
      , m_element_sizes(element_sizes)
      , m_element_offsets(element_offsets)
    { }

    /// Constructor
    AXOM_HOST_DEVICE
    PolyhedronData(const PolyhedronData &obj)
      : m_subelement_conn(obj.m_subelement_conn)
      , m_subelement_sizes(obj.m_subelement_sizes)
      , m_subelement_offsets(obj.m_subelement_offsets)
      , m_element_conn(obj.m_element_conn)
      , m_element_sizes(obj.m_element_sizes)
      , m_element_offsets(obj.m_element_offsets)
    { }

    ConnectivityView m_subelement_conn;
    ConnectivityView m_subelement_sizes;
    ConnectivityView m_subelement_offsets;
    ConnectivityView m_element_conn;
    ConnectivityView m_element_sizes;
    ConnectivityView m_element_offsets;
  };

  /**
   * \brief This struct provides data about Zone i's shape.
   */
  struct PolyhedronShape
  {
    constexpr static IndexType MaximumNumberOfIds = 20 * 3;

    AXOM_HOST_DEVICE constexpr static bool is_polyhedral() { return true; }
    AXOM_HOST_DEVICE constexpr static int id() { return Polyhedron_ShapeID; }

    /// Constructor.
    AXOM_HOST_DEVICE PolyhedronShape(const PolyhedronData &obj, axom::IndexType zi)
      : m_data(obj)
      , m_zoneIndex(zi)
      , m_ids()
    { }

    /// This implementation does not return unique number of nodes.
    AXOM_HOST_DEVICE IndexType numberOfNodes() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        nnodes += getFace(f).size();
      }
      return nnodes;
    }

    AXOM_HOST_DEVICE IndexType numberOfFaces() const
    {
      return m_data.m_element_sizes[m_zoneIndex];
    }

    AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
    {
      return getFace(faceIndex).size();
    }

    AXOM_HOST_DEVICE ConnectivityView getIds() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        const auto faceIds = getFace(f);
        for(axom::IndexType i = 0; i < faceIds.size(); i++)
        {
          if(nnodes + 1 <= MaximumNumberOfIds)
            m_ids[nnodes++] = faceIds[i];
          else
          {
#if !defined(AXOM_DEVICE_CODE)
            SLIC_ERROR("m_ids is not large enough to hold all node ids.");
#endif
            break;
          }
        }
      }
      return ConnectivityView(m_ids.m_data, nnodes);
    }

    AXOM_HOST_DEVICE ConnectivityView getUniqueIds() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        const auto faceIds = getFace(f);
        for(axom::IndexType i = 0; i < faceIds.size(); i++)
        {
          if(!find(m_ids.m_data, nnodes, faceIds[i]))
          {
            if(nnodes + 1 <= MaximumNumberOfIds)
              m_ids[nnodes++] = faceIds[i];
            else
            {
#if !defined(AXOM_DEVICE_CODE)
              SLIC_ERROR("m_ids is not large enough to hold all node ids.");
#endif
              break;
            }
          }
        }
      }
      return ConnectivityView(m_ids.m_data, nnodes);
    }

    AXOM_HOST_DEVICE ConnectivityView getFace(int faceIndex) const
    {
      const ConnectivityView element_face_ids(
        m_data.m_element_conn.data() + m_data.m_element_offsets[m_zoneIndex],
        m_data.m_element_sizes[m_zoneIndex]);
      const auto faceId = element_face_ids[faceIndex];

      return ConnectivityView(
        m_data.m_subelement_conn.data() + m_data.m_subelement_offsets[faceId],
        m_data.m_subelement_sizes[faceId]);
    }

  private:
    AXOM_HOST_DEVICE bool find(const ConnectivityType *arr,
                               axom::IndexType n,
                               ConnectivityType value) const
    {
      bool found = false;
      for(axom::IndexType i = 0; i < n && !found; i++)
      {
        found = arr[i] == value;
      }
      return found;
    }

    PolyhedronData m_data;
    IndexType m_zoneIndex {0};
    mutable axom::StackArray<ConnectivityType, MaximumNumberOfIds> m_ids;
  };
  //----------------------------------------------------------------------------

  using ShapeType = PolyhedronShape;

  /**
   * \brief Constructor.
   */
  UnstructuredTopologyPolyhedralView(const ConnectivityView &subelement_conn,
                                     const ConnectivityView &subelement_sizes,
                                     const ConnectivityView &subelement_offsets,
                                     const ConnectivityView &element_conn,
                                     const ConnectivityView &element_sizes,
                                     const ConnectivityView &element_offsets)
    : m_data(subelement_conn,
             subelement_sizes,
             subelement_offsets,
             element_conn,
             element_sizes,
             element_offsets)
  { }

  /**
   * \brief Return the number of zones in the mesh.
   *
   * \return The number of zones.
   */
  IndexType numberOfZones() const { return m_data.m_element_sizes.size(); }

  /**
   * \brief Return the size of the connectivity.
   *
   * \return The size of the connectivity.
   */
  IndexType connectivitySize() const
  {
    return m_data.element_conn.size();
  }

  /**
   * \brief Return the dimension of the shape.
   *
   * \return The dimension of the shape.
   */
  static constexpr int dimension() { return 3; }

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

    const PolyhedronData sd(m_data);
    axom::for_all<ExecSpace>(
      0,
      nzones,
      AXOM_LAMBDA(int zoneIndex) {
        const PolyhedronShape shape(sd, zoneIndex);
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
  void for_selected_zones(const ViewType &selectedIdsView, FuncType &&func) const
  {
    const auto nSelectedZones = selectedIdsView.size();

    ViewType idsView(selectedIdsView);
    const PolyhedronData sd(m_data);
    axom::for_all<ExecSpace>(
      0,
      nSelectedZones,
      AXOM_LAMBDA(auto selectIndex) {
        const auto zoneIndex = idsView[selectIndex];
        const PolyhedronShape shape(sd, zoneIndex);
        func(selectIndex, zoneIndex, shape);
      });
  }

private:
  PolyhedronData m_data;
};

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
