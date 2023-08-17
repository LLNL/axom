// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_MESH_VIEW_UTIL_H_
#define QUEST_MESH_VIEW_UTIL_H_

#include "axom/config.hpp"
#include "axom/core.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint_mesh.hpp"

#include <memory>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <vector>

namespace axom
{
namespace quest
{

/**
   \brief Utility for accessing views into a blueprint mesh,
   for structured mesh with explicit coordinates.

   Views are single-domain-specific.  They don't apply to multi-domain meshes.
   They are also topology specific, with the topology name given
   to the constructor.

   TODO: Figure out a more appropriate place for this file.
   It's only in axom/quest because the initial need was there.
*/
template<int DIM>
class MeshViewUtil
{
public:

  //@{
  //@name Setting up
  /*!
    @brief Constructor
    @param [in] bpDomain Blueprint single domain.
    @param [in] topologyName Name of topology in the domain.

    If \a topologyName is omitted, use the first topology,

    The topology dimension must match DIM.
  */
  MeshViewUtil(conduit::Node& bpDomain, const std::string& topologyName = "")
    : m_dom(&bpDomain)
    , m_topologyName(topologyName)
    {
      if(m_topologyName.empty())
      {
        m_topologyName = bpDomain.fetch_existing("topologies")[0].name();
      }

      const std::string topoPath = "topologies/" + m_topologyName;
      SLIC_ASSERT(m_dom->has_path(topoPath));
      m_topology = &m_dom->fetch_existing(topoPath);

      SLIC_ASSERT(m_topology->fetch_existing("type").as_string() == "structured");

      const std::string coordsetPath =
        "coordsets/" + m_dom->fetch_existing(topoPath + "/coordset").as_string();
      m_coordset = &m_dom->fetch_existing(coordsetPath);

      SLIC_ASSERT( m_coordset->fetch_existing("type").as_string() == "explicit");

      int dim = conduit::blueprint::mesh::coordset::dims(*m_coordset);
      if(dim != DIM)
      {
        SLIC_ERROR(axom::fmt::format("MeshViewUtil template parameter DIM={} doesn't match mesh topology dimension ({})",
                                     DIM, dim));
      }

      m_cdom = m_dom;
      m_ctopology = m_topology;
      m_ccoordset = m_coordset;
    }

  /*!
    @brief Constructor
    @param [in] bpDomain const Blueprint single domain.
    @param [in] topologyName Name op topology in the domain.
  */
  MeshViewUtil(const conduit::Node& bpDomain, const std::string& topologyName)
    : m_cdom(&bpDomain)
    , m_topologyName(topologyName)
    {
      const std::string topoPath = "topologies/" + m_topologyName;
      SLIC_ASSERT(m_dom->has_path(topoPath));
      m_ctopology = &m_cdom->fetch_existing(topoPath);

      SLIC_ASSERT(m_ctopology->fetch_existing("type").as_string() == "structured");

      const std::string coordsetPath =
        "coordsets/" + m_dom->fetch_existing(topoPath + "/coordset").as_string();
      m_ccoordset = &m_cdom->fetch_existing(coordsetPath);

      SLIC_ASSERT( m_ccoordset->fetch_existing("type").as_string() == "explicit");

      int dim = conduit::blueprint::mesh::coordset::dims(*m_ccoordset);
      if(dim != DIM)
      {
        SLIC_ERROR(axom::fmt::format("MeshViewUtil template parameter DIM={} doesn't match mesh topology dimension ({})",
                                     DIM, dim));
      }
    }
  //@}


  //@{
  //!@name Blueprint data organization
  /*!
    @brief Get the named topology.

    @return the Conduit node at "topologies/<topologyName>".
  */
  const conduit::Node& getTopology()
  {
    return *m_topology;
  }
  const conduit::Node& getTopology() const
  {
    return *m_ctopology;
  }

  /*!
    @brief Get the coordinates conduit::Node for the named topology.
  */
  conduit::Node& getCoordSet()
  {
    return *m_coordset;
  }
  const conduit::Node& getCoordSet() const
  {
    return *m_ccoordset;
  }

  /*!
    @brief Get the spatial dimension of the named topology.
  */
  conduit::index_t getTopologyDim() const
  {
    return DIM;
  }
  //@}


  //@{
  //!@name Sizes and shapes
  /*!
    @brief Get the number of cells in each direction of a blueprint single domain.

    @param domId Index of domain
    @lengths Space for dimension() numbers.
  */
  axom::StackArray<axom::IndexType, DIM> getDomainShape() const
  {
    const conduit::Node& dimsNode = m_ctopology->fetch_existing("elements/dims");
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int i = 0; i < DIM; ++i)
    {
      rval[i] = dimsNode[i].as_int();
    }
    return rval;
  }

  axom::StackArray<axom::IndexType, DIM> getRealExtents(const std::string& association)
  {
    axom::StackArray<axom::IndexType, DIM> rval = getDomainShape();
    if(association == "vertex")
    {
      for(int d=0; d<DIM; ++d)
      {
        ++rval[d];
      }
    }
    else if(association == "element")
    {
      // Nothing to do.
    }
    else
    {
      SLIC_ERROR("MeshVieuUtil only supports element and vertex data association for now.");
    }

    return rval;
  }

  /*!
    @brief Return number of points, excluding ghosts.
  */
  axom::IndexType getCoordsCount() const
  {
    auto domainShape = getDomainShape();
    axom::IndexType rval = 0;
    for(int d=0; d<DIM; ++d)
    {
      rval += 1 + domainShape[d];
    }
    return rval;
  }

  /*!
    @brief Return number of points, including ghosts.
  */
  axom::IndexType getCoordsCountWithGhosts() const
  {
    const conduit::Node& valuesNode = m_ccoordset->fetch_existing("values");
    axom::IndexType rval = valuesNode[0].dtype().number_of_elements();
    return rval;
  }

  /*!
    @brief Return coordinates data allocated shape of coords data,
    (includes ghosts).
  */
  axom::StackArray<axom::IndexType, DIM> getCoordsShapeWithGhosts() const
  {
    const conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    const axom::IndexType coordValuesCount =
      valuesNode[0].dtype().number_of_elements();

    // Shape of allocated memory, including ghosts.
    axom::StackArray<axom::IndexType, DIM> memShape;

    auto stridesPtr = getCoordsStridesPtr();
    if(stridesPtr)
    {
      axom::StackArray<axom::IndexType, DIM> strides;
      for(int d=0; d<DIM; ++d)
      {
        strides[d] = stridesPtr[d];
        memShape[d] = (d < DIM-1 ? stridesPtr[d+1] : coordValuesCount) /
          stridesPtr[d];
      }
    }
    else
    {
      // No strides implies no ghosts, so memory shape is domain shape.
      memShape = getDomainShape();
      for(int d=0; d<DIM; ++d)
      {
        memShape[d] += 1;
      }
    }
    return memShape;
  }

  /*!
    @brief Get the strides of the coordindates data for the
    named topology.  If no strides, return null.
  */
  const conduit::int32* getCoordsStridesPtr() const
  {
    const conduit::Node& topologyDims = m_ctopology->fetch_existing("elements/dims");
    const conduit::int32* rval = nullptr;
    if(topologyDims.has_child("strides"))
    {
      rval = topologyDims.fetch_existing("strides").as_int32_ptr();
    }
    return rval;
  }

  /*!
    @brief Get the strides of the coordinates data array.
    If strides are not specified, assume direction 0 is
    fastest (Conduit's default striding).
  */
  axom::StackArray<axom::IndexType, DIM> getCoordsStrides() const
  {
    axom::StackArray<axom::IndexType, DIM> rval;
    auto stridesPtr = getCoordsStridesPtr();
    if(stridesPtr)
    {
      for(int d=0; d<DIM; ++d)
      {
        rval[d] = stridesPtr[d];
      }
    }
    else
    {
      auto domainShape = getDomainShape();
      axom::IndexType tmpStride = 1;
      for(int d=0; d<DIM; ++d)
      {
        rval[d] = tmpStride;
        tmpStride *= 1 + domainShape[d];
      }
    }
    return rval;
  }

  /*!
    @brief Get the offsets of the coordindates data for the
    named topology.  If no offsets, return null.
  */
  const conduit::int32* getCoordsOffsetsPtr() const
  {
    const conduit::Node& topologyDims = m_ctopology->fetch_existing("elements/dims");
    const conduit::int32* rval = nullptr;
    if(topologyDims.has_child("offsets"))
    {
      rval = topologyDims.fetch_existing("offsets").as_int32_ptr();
    }
    return rval;
  }

  axom::StackArray<axom::IndexType, DIM> getCoordsOffsets() const
  {
    auto offsetsPtr = getCoordsOffsetsPtr();
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d=0; d<DIM; ++d)
    {
      rval[d] = offsetsPtr[d];
    }
    return rval;
  }

  /*!
    @brief Return number of points, including ghosts.
  */
  axom::IndexType getFieldCountWithGhosts(const std::string& fieldName) const
  {
    const conduit::Node& valuesNode = m_cdom->fetch_existing("field/" + fieldName + "values");
    axom::IndexType rval = valuesNode.dtype().number_of_elements();
    return rval;
  }

  /*!
    @brief Get the strides of a named field data.
    If the field has no strides, return null.
  */
  const conduit::int32* getFieldStridesPtr(const std::string& fieldName) const
  {
    const conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    const conduit::int32* rval = fieldNode.has_child("strides") ?
      fieldNode.fetch_existing("strides").as_int32_ptr() : nullptr;
    return rval;
  }

  /*!
    @brief Get the strides of a named field data.
    If strides are not specified, assume direction 0 is
    fastest (Conduit's default striding).
  */
  axom::StackArray<axom::IndexType, DIM> getFieldStrides(const std::string& fieldName) const
  {
    const conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    const bool atVertex = fieldNode.fetch_existing("association").as_string() == "vertex";
    const conduit::int32* stridesPtr = fieldNode.has_child("strides") ?
      fieldNode.fetch_existing("strides").as_int32_ptr() : nullptr;

    axom::StackArray<axom::IndexType, DIM> rval;
    if(stridesPtr)
    {
      for(int d=0; d<DIM; ++d)
      {
        rval[d] = stridesPtr[d];
      }
    }
    else
    {
      auto domainShape = getDomainShape();
      axom::IndexType tmpStride = 1;
      for(int d=0; d<DIM; ++d)
      {
        rval[d] = tmpStride;
        tmpStride *= domainShape[d] + atVertex;
      }
    }
    return rval;
  }
  //@}


  //@{
  //!@name Data views

  /*!
    @brief Return views to all the coordinates.
  */
  axom::Array<axom::ArrayView<double, DIM>> getCoordsViews()
  {
    conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    const axom::IndexType coordValuesCount =
      valuesNode[0].dtype().number_of_elements();

    axom::Array<axom::ArrayView<double, DIM>> rval(DIM);

    // Shape of allocated memory, including ghosts.
    axom::StackArray<axom::IndexType, DIM> memShape;

    auto stridesPtr = getCoordsStridesPtr();
    if(stridesPtr)
    {
      axom::StackArray<axom::IndexType, DIM> strides;
      for(int d=0; d<DIM; ++d)
      {
        strides[d] = stridesPtr[d];
        memShape[d] = (d < DIM-1 ? stridesPtr[d+1] : coordValuesCount) /
          stridesPtr[d];
      }
      for(int d=0; d<DIM; ++d)
      {
        double* dataValues = valuesNode[d].as_double_ptr();
        rval[d] = axom::ArrayView<double, DIM>(dataValues, memShape, strides);
      }
      assert(strides == getCoordsStrides());
    }
    else
    {
      // No strides implies no ghosts, so memory shape is domain shape.
      memShape = getDomainShape();
      axom::StackArray<axom::IndexType, DIM> strides = getCoordsStrides();
      for(int d=0; d<DIM; ++d)
      {
        double* dataValues = valuesNode[d].as_double_ptr();
        rval[d] = axom::ArrayView<double, DIM>(dataValues, memShape, strides);
      }
    }

    return rval;
  }

  /*!
    @brief Return view to a scalar field variable.
  */
  axom::ArrayView<double, DIM> getFieldView(const std::string& fieldName)
  {
    conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    conduit::Node& valuesNode = fieldNode.fetch_existing("values");
    const axom::IndexType valuesCount = valuesNode.dtype().number_of_elements();

    axom::ArrayView<double, DIM> rval;

    // Shape of allocated memory, including ghosts.
    axom::StackArray<axom::IndexType, DIM> memShape;

    auto stridesPtr = fieldNode.has_child("strides") ?
      fieldNode.fetch_existing("strides").as_int32_ptr() : nullptr;
    axom::StackArray<axom::IndexType, DIM> strides = getFieldStrides(fieldName);

    double* dataValues = valuesNode.as_double_ptr();// TODO: Use template parameter T.

    if(stridesPtr)
    {
      for(int d=0; d<DIM; ++d)
      {
        memShape[d] = (d < DIM-1 ? strides[d+1] : valuesCount) /
          strides[d];
      }
    }
    else
    {
      // No strides implies no ghosts, so memory shape is domain shape.
      memShape = getDomainShape();
      for(int d=0; d<DIM; ++d)
      {
        ++memShape[d];
      }
    }
    rval = axom::ArrayView<double, DIM>(dataValues, memShape, strides);

    return rval;
  }
  //@}

  //@{
  //!@name Creating data

  /*!
    @brief Create a new scalar nodal data field.
    @param [in] fieldName
    @param [in] dtype Conduit data type to put in the field.  Must be at least
                big enough for the strides specified.
    @param [in] association "vertex" or "element"
    @param [in] strides Data strides.  Set to zero for no ghosts and default strides.
    @param [in] offsets Data index offsets.
  */
  void createNodalField( const std::string& fieldName, const std::string& association,
                         const conduit::DataType& dtype,
                         const axom::StackArray<axom::IndexType, DIM>& strides,
                         const axom::StackArray<axom::IndexType, DIM>& offsets)
    {
      if(m_dom->has_path("fields/" + fieldName))
      {
        SLIC_ERROR(axom::fmt::format("Cannot create field {}.  It already exists.", fieldName));
      }
      if(association != "vertex" && association != "element")
      {
        SLIC_ERROR(axom::fmt::format("Not yet supporting association '{}'.", association));
      }

      auto realExtents = getRealExtents(association);

      conduit::Node& fieldNode = m_dom->fetch_existing("fields")[fieldName];
      fieldNode["association"] = association;
      fieldNode["topology"] = m_topologyName;

      {
        fieldNode["strides"].set(conduit::DataType::int32(DIM));
        std::int32_t* tmpPtr = fieldNode["strides"].as_int32_ptr();
        for(int d=0; d<DIM; ++d)
        {
          tmpPtr[d] = strides[d];
        }
      }

      {
        fieldNode["offsets"].set(conduit::DataType::int32(DIM));
        std::int32_t* tmpPtr = fieldNode["offsets"].as_int32_ptr();
        for(int d=0; d<DIM; ++d)
        {
          tmpPtr[d] = offsets[d];
        }
      }

      axom::IndexType slowDir = 0;
      for(int d=0; d<DIM; ++d)
      {
        if(strides[slowDir] < strides[d])
        {
          slowDir = d;
        }
      }
      auto extras = dtype.number_of_elements() -
        strides[slowDir]*(offsets[slowDir] + realExtents[slowDir]);
      if(extras < 0)
      {
        SLIC_ERROR(axom::fmt::format("Insufficient space allocated for data, given the offsets and strides."));
      }

      fieldNode["values"].set(dtype);
    }
  //@}

private:
  conduit::Node* m_dom = nullptr;
  conduit::Node* m_topology = nullptr;
  conduit::Node* m_coordset = nullptr;
  const conduit::Node* m_cdom = nullptr;
  const conduit::Node* m_ctopology = nullptr;
  const conduit::Node* m_ccoordset = nullptr;
  std::string m_topologyName;
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_MESH_VIEW_UTIL_H_
