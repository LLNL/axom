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

template <typename IndexT>
class UnstructuredTopologyPolyhedralView
{
public:
  using IndexType = IndexT;

  struct PolyhedronData
  {
    AXOM_HOST_DEVICE
    PolyhedronData(axom::ArrayView<IndexType> subelement_conn,
              axom::ArrayView<IndexType> subelement_sizes,
              axom::ArrayView<IndexType> subelement_offsets,
              axom::ArrayView<IndexType> element_conn,
              axom::ArrayView<IndexType> element_sizes,
              axom::ArrayView<IndexType> element_offsets) :
      m_subelement_conn(subelement_conn), m_subelement_sizes(subelement_sizes), m_subelement_offsets(subelement_offsets),
      m_element_conn(element_conn), m_element_sizes(element_sizes), m_element_offsets(element_offsets)
    {
    }

    AXOM_HOST_DEVICE
    PolyhedronData(const PolyhedronData &obj) :
      m_subelement_conn(obj.m_subelement_conn), m_subelement_sizes(obj.m_subelement_sizes), m_subelement_offsets(obj.m_subelement_offsets),
      m_element_conn(obj.m_element_conn), m_element_sizes(obj.m_element_sizes), m_element_offsets(obj.m_element_offsets)
    {
    }

    axom::ArrayView<IndexType> m_subelement_conn;
    axom::ArrayView<IndexType> m_subelement_sizes;
    axom::ArrayView<IndexType> m_subelement_offsets;
    axom::ArrayView<IndexType> m_element_conn;
    axom::ArrayView<IndexType> m_element_sizes;
    axom::ArrayView<IndexType> m_element_offsets;
  };

  // Can we provide a way to provide data about Zone i's shape?
  struct PolyhedronShape
  {
    constexpr static IndexType MaximumNumberOfIds = 20 * 3;

    AXOM_HOST_DEVICE constexpr bool is_polyhedral() { return true; }

    AXOM_HOST_DEVICE PolyhedronShape(const PolyhedronData &obj, axom::IndexType zi) : m_data(obj), m_zoneIndex(zi), m_ids()
    {
    }

    /// This implementation does not return unique number of nodes.
    AXOM_HOST_DEVICE IndexType numberOfNodes() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
        nnodes += getFace(f).size();
    }

    AXOM_HOST_DEVICE IndexType numberOfFaces() const
    {
      return m_data.m_element_sizes[m_zoneIndex];
    }

    AXOM_HOST_DEVICE IndexType numberOfNodesInFace(int faceIndex) const
    {
      return getFace(faceIndex).size();
    }

    AXOM_HOST_DEVICE axom::ArrayView<IndexType> getIds() const
    {
      axom::IndexType nnodes = 0;
      const auto nFaces = numberOfFaces();
      for(axom::IndexType f = 0; f < nFaces; f++)
      {
        const auto faceIds = getFace(f);
        for(axom::IndexType i = 0; i < faceIds.size(); i++)
        {
          if(nnodes < MaximumNumberOfIds)
            m_ids[nnodes++] = faceIds[i];
        }
      }
      return axom::ArrayView<IndexType>(m_ids.m_data, nnodes);
    }

    AXOM_HOST_DEVICE axom::ArrayView<IndexType> getFace(int faceIndex) const
    {
      const axom::ArrayView<IndexType> element_face_ids(m_data.m_element_conn.data() + m_data.m_element_offsets[m_zoneIndex], m_data.m_element_sizes[m_zoneIndex]);
      const auto faceId = element_face_ids[faceIndex];

      return axom::ArrayView<IndexType>(m_data.m_subelement_conn.data() + m_data.m_subelement_offsets[faceId], m_data.m_subelement_sizes[faceId]);
    }

  private:
    PolyhedronData m_data;
    IndexType m_zoneIndex {0};
    mutable axom::StackArray<IndexType,MaximumNumberOfIds> m_ids;
  };
  //----------------------------------------------------------------------------

  using ShapeType = PolyhedronShape;

  UnstructuredTopologyPolyhedralView(const axom::ArrayView<IndexType> &subelement_conn,
                                     const axom::ArrayView<IndexType> &subelement_sizes,
                                     const axom::ArrayView<IndexType> &subelement_offsets,
                                     const axom::ArrayView<IndexType> &element_conn,
                                     const axom::ArrayView<IndexType> &element_sizes,
                                     const axom::ArrayView<IndexType> &element_offsets) :
    m_data(subelement_conn, subelement_sizes, subelement_offsets,
           element_conn, element_sizes, element_offsets)
  {
  }

  IndexType numberOfZones() const
  {
    return m_data.m_element_sizes.size();
  }

  template <typename ExecSpace, typename FuncType>
  void for_all_zones(FuncType &&func) const
  {
    const auto nzones = numberOfZones();

    const PolyhedronData sd(m_data);
    axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(int zoneIndex)
    {
      const PolyhedronShape shape(sd, zoneIndex);
      func(zoneIndex, shape);
    });
  }

private:
  PolyhedronData m_data;
};

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
