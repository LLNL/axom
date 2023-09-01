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
#include "conduit_blueprint_mcarray.hpp"
#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
  #include "conduit_relay_mpi_io_blueprint.hpp"
#endif

#include <memory>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <vector>

namespace axom
{
namespace quest
{
namespace internal
{
template <typename T>
constexpr T BADINDEX = std::numeric_limits<T>::max();
template <typename T, int DIM>
inline axom::StackArray<T, DIM> makeStackArray(T v = std::numeric_limits<T>::max())
{
  axom::StackArray<T, DIM> rval;
  for(int d = 0; d < DIM; ++d)
  {
    rval[d] = v;
  }
  return rval;
}
template <typename T, int DIM>
inline axom::StackArray<T, DIM> makeStackArray(const T* v)
{
  axom::StackArray<T, DIM> rval;
  for(int d = 0; d < DIM; ++d)
  {
    rval[d] = v[d];
  }
  return rval;
}
}  // namespace internal

/**
   \brief Utility for high-level access into a blueprint mesh,
   for structured mesh with explicit coordinates.

   Blueprint mesh data is sufficient but sparse, leaving users
   to compute a number of intermediate data to get to high-level
   data of interest, such as views into array data.  This class
   encapsulates those common utilities.

   Views are single-domain-specific.  They don't apply to multi-domain
   meshes.  They are also topology specific, with the topology name
   given to the constructor.

   TODO: Figure out if there's a better place for this utility.
   It's only in axom/quest because the initial need was there.
*/
template <int DIM, axom::MemorySpace MemSpace>
class MeshViewUtil
{
public:
  using CoordsViewsType =
    axom::StackArray<axom::ArrayView<double, DIM, MemSpace>, DIM>;
  using ConstCoordsViewsType =
    axom::StackArray<axom::ArrayView<const double, DIM, MemSpace>, DIM>;

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
    m_topology = &m_dom->fetch_existing(topoPath);

    const std::string coordsetPath =
      "coordsets/" + m_dom->fetch_existing(topoPath + "/coordset").as_string();
    m_coordset = &m_dom->fetch_existing(coordsetPath);

    m_cdom = m_dom;
    m_ctopology = m_topology;
    m_ccoordset = m_coordset;

    if(!isValid())
    {
      SLIC_ERROR(
        "Invalid domain in MeshViewUtil.  Domain must be blueprint-valid,"
        " structured, with exlicit coordinates and dimension equal to "
        " the template DIM parameter.");
    }

    computeCoordsDataLayout();
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
    SLIC_ASSERT(m_cdom->has_path(topoPath));
    m_ctopology = &m_cdom->fetch_existing(topoPath);

    const std::string coordsetPath =
      "coordsets/" + m_cdom->fetch_existing(topoPath + "/coordset").as_string();
    m_ccoordset = &m_cdom->fetch_existing(coordsetPath);

    if(!isValid())
    {
      SLIC_ERROR(
        "Invalid domain in MeshViewUtil.  Domain must be blueprint-valid,"
        " structured, with exlicit coordinates. and dimension equal to "
        " the template DIM parameter.");
    }

    computeCoordsDataLayout();
  }

  /*!
    @brief Whether the domain is a valid one for this utility.

    In addition to checking for blueprint validity (unless \a
    checkBlueprint is false), this method checks against its own
    requirements.
  */
  bool isValid(bool checkBlueprint = false) const
  {
    bool rval = true;
    if(checkBlueprint)
    {
      conduit::Node info;
#ifdef AXOM_USE_MPI
      rval = rval &&
        conduit::blueprint::mpi::verify("mesh", *m_cdom, info, MPI_COMM_WORLD);
#else
      rval = rval && conduit::blueprint::verify("mesh", *m_cdom, info);
#endif
    }
    rval =
      rval && (m_ctopology->fetch_existing("type").as_string() == "structured");
    rval =
      rval && (m_ccoordset->fetch_existing("type").as_string() == "explicit");
    rval =
      rval && (conduit::blueprint::mesh::coordset::dims(*m_ccoordset) == DIM);
    return true;
  }
  //@}

  //@{
  //!@name Blueprint data organization
  /*!
    @brief Get the named topology.

    @return the Conduit node at "topologies/<topologyName>".
  */
  const conduit::Node& getTopology() { return *m_topology; }
  const conduit::Node& getTopology() const { return *m_ctopology; }

  /*!
    @brief Get the coordinates conduit::Node for the named topology.
  */
  conduit::Node& getCoordSet() { return *m_coordset; }
  const conduit::Node& getCoordSet() const { return *m_ccoordset; }

  /*!
    @brief Get the spatial dimension of the named topology.
  */
  conduit::index_t getTopologyDim() const { return DIM; }
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
    const conduit::Node& dimsNode =
      m_ctopology->fetch_existing("elements/dims");
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int i = 0; i < DIM; ++i)
    {
      rval[i] = dimsNode[i].as_int();
    }
    return rval;
  }

  //! @brief Return the real (ghost-free) extents of mesh data.
  axom::StackArray<axom::IndexType, DIM> getRealExtents(const std::string& association)
  {
    axom::StackArray<axom::IndexType, DIM> rval = getDomainShape();
    if(association == "vertex")
    {
      for(int d = 0; d < DIM; ++d)
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
      SLIC_ERROR(
        axom::fmt::format("MeshVieuUtil only supports element and vertex data "
                          "association for now, not '{}'.",
                          association));
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
    for(int d = 0; d < DIM; ++d)
    {
      rval += 1 + domainShape[d];
    }
    return rval;
  }

  /*!
    @brief Return the array strides for ghost-free nodal coordinates.
  */
  axom::StackArray<axom::IndexType, DIM> getCoordsStrides() const
  {
    return m_coordsStrides;
  }

  /*!
    @brief Return the array index offsets for nodal coordinates.
  */
  axom::StackArray<axom::IndexType, DIM> getCoordsOffsets() const
  {
    return m_coordsOffsets;
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
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = stridesPtr[d];
        memShape[d] =
          (d < DIM - 1 ? stridesPtr[d + 1] : coordValuesCount) / stridesPtr[d];
      }
    }
    else
    {
      // No strides implies no ghosts, so memory shape is domain shape.
      memShape = getDomainShape();
      for(int d = 0; d < DIM; ++d)
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
    const conduit::Node& topologyDims =
      m_ctopology->fetch_existing("elements/dims");
    const conduit::int32* rval = nullptr;
    if(topologyDims.has_child("strides"))
    {
      rval = topologyDims.fetch_existing("strides").as_int32_ptr();
    }
    return rval;
  }

  /*!
    @brief Return number of points, including ghosts.
  */
  axom::IndexType getFieldCountWithGhosts(const std::string& fieldName) const
  {
    const conduit::Node& valuesNode =
      m_cdom->fetch_existing("field/" + fieldName + "values");
    axom::IndexType rval = valuesNode.dtype().number_of_elements();
    return rval;
  }

  /*!
    @brief Get the strides of a named Blueprint field.
    If strides are not specified, assume direction 0 is
    fastest (Conduit's default striding).
  */
  axom::StackArray<axom::IndexType, DIM> getFieldStrides(
    const std::string& fieldName) const
  {
    const conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    const bool atVertex =
      fieldNode.fetch_existing("association").as_string() == "vertex";
    const conduit::int32* stridesPtr = fieldNode.has_child("strides")
      ? fieldNode.fetch_existing("strides").as_int32_ptr()
      : nullptr;

    axom::StackArray<axom::IndexType, DIM> rval;
    if(stridesPtr)
    {
      for(int d = 0; d < DIM; ++d)
      {
        rval[d] = stridesPtr[d];
      }
    }
    else
    {
      auto domainShape = getDomainShape();
      axom::IndexType tmpStride = 1;
      for(int d = 0; d < DIM; ++d)
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

  //!@brief Return the views of the DIM coordinates component data.
  CoordsViewsType getCoordsViews(bool withGhosts = false)
  {
    conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    axom::StackArray<axom::ArrayView<double, DIM, MemSpace>, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      auto* dataPtr = valuesNode[d].as_double_ptr();
      rval[d] = axom::ArrayView<double, DIM, MemSpace>(dataPtr,
                                                       m_coordsMemShape,
                                                       m_coordsStrides);
    }

    if(withGhosts == false)
    {
      axom::StackArray<axom::IndexType, DIM> offsets =
        conduitIndicesToMultidimIndices(
          m_ctopology->fetch_existing("elements/dims"),
          "offsets",
          0);
      axom::StackArray<axom::IndexType, DIM> counts = m_domainShape;
      for(int d = 0; d < DIM; ++d)
      {
        ++counts[d];
      }
      for(int d = 0; d < DIM; ++d)
      {
        auto rval1 = rval[d];
        rval[d] = rval1.subspan(offsets, counts);
      }
    }

    return rval;
  }

  //!@brief Return the views of the DIM coordinates component data.
  ConstCoordsViewsType getConstCoordsViews(bool withGhosts = false) const
  {
    const conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    axom::StackArray<axom::ArrayView<const double, DIM, MemSpace>, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      auto* dataPtr = valuesNode[d].as_double_ptr();
      rval[d] = axom::ArrayView<const double, DIM, MemSpace>(dataPtr,
                                                             m_coordsMemShape,
                                                             m_coordsStrides);
    }

    if(withGhosts == false)
    {
      axom::StackArray<axom::IndexType, DIM> offsets =
        conduitIndicesToMultidimIndices(
          m_ctopology->fetch_existing("elements/dims"),
          "offsets",
          0);
      axom::StackArray<axom::IndexType, DIM> counts = m_domainShape;
      for(int d = 0; d < DIM; ++d)
      {
        ++counts[d];
      }
      for(int d = 0; d < DIM; ++d)
      {
        auto rval1 = rval[d];
        rval[d] = rval1.subspan(offsets, counts);
      }
    }

    return rval;
  }

  /*!
    @brief Return view to a scalar field variable.

    WARNING: The view returned has an allocator id determined by
    \a MemSpace, regardless of the memory type.

    WARNING: Assuming, without checking, that the field contains
    data of type \a T.  User is responsible for using the correct
    type.
  */
  template <typename T>
  axom::ArrayView<T, DIM, MemSpace> getFieldView(const std::string& fieldName,
                                                 bool withGhosts = false)
  {
    if(!m_dom)
    {
      SLIC_ERROR(
        "Cannot use getFieldView from MeshViewUtil initialized with a const "
        "conduit::Node.  Use getConstFieldView.");
    }
    conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    const std::string association =
      fieldNode.fetch_existing("association").as_string();
    conduit::Node& valuesNode = fieldNode.fetch_existing("values");
    const axom::IndexType valuesCount = valuesNode.dtype().number_of_elements();

    SLIC_ASSERT_MSG(
      association == "vertex" || association == "element",
      "MeshViewUtil only supports vertex and element-based fields right now.");
    const bool onVertex = association == "vertex";

    axom::StackArray<axom::IndexType, DIM> strides =
      conduitIndicesToMultidimIndices(fieldNode, "strides");
    if(fieldNode.has_child("strides"))
    {
      const conduit::int32* stridesPtr =
        fieldNode.fetch_existing("strides").as_int32_ptr();
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = stridesPtr[d];
      }
    }
    else
    {
      axom::IndexType tmpStride = 1;
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = tmpStride;
        tmpStride *= m_domainShape[d] + onVertex;
      }
    }

    axom::StackArray<axom::IndexType, DIM> shape;
    for(int d = 0; d < DIM - 1; ++d)
    {
      shape[d] = strides[d + 1] / strides[d];
    }
    shape[DIM - 1] = valuesCount / strides[DIM - 1];

    T* dataPtr = static_cast<T*>(valuesNode.data_ptr());
    axom::ArrayView<double, DIM, MemSpace> rval(dataPtr, shape, strides);

    if(withGhosts == false)
    {
      axom::StackArray<axom::IndexType, DIM> offsets =
        conduitIndicesToMultidimIndices(fieldNode, "offsets", 0);
      axom::StackArray<axom::IndexType, DIM> counts = m_domainShape;
      for(int d = 0; d < DIM; ++d)
      {
        counts[d] += onVertex;
      }
      auto rval1 = rval;
      rval = rval1.subspan(offsets, counts);
    }

    return rval;
  }

  /*!
    @brief Return view to a scalar field variable.

    WARNING: The view returned has an allocator id determined by
    \a MemSpace, regardless of the memory type.

    WARNING: Assuming, without checking, that the field contains
    data of type \a T.  User is responsible for using the correct
    type.
  */
  template <typename T>
  axom::ArrayView<const T, DIM, MemSpace> getConstFieldView(
    const std::string& fieldName,
    bool withGhosts = false)
  {
    const conduit::Node& fieldNode =
      m_cdom->fetch_existing("fields/" + fieldName);
    const std::string association =
      fieldNode.fetch_existing("association").as_string();
    const conduit::Node& valuesNode = fieldNode.fetch_existing("values");
    const axom::IndexType valuesCount = valuesNode.dtype().number_of_elements();

    SLIC_ASSERT_MSG(
      association == "vertex" || association == "element",
      "MeshViewUtil only supports vertex and element-based fields right now.");
    const bool onVertex = association == "vertex";

    axom::StackArray<axom::IndexType, DIM> strides =
      conduitIndicesToMultidimIndices(fieldNode, "strides");
    if(fieldNode.has_child("strides"))
    {
      const conduit::int32* stridesPtr =
        fieldNode.fetch_existing("strides").as_int32_ptr();
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = stridesPtr[d];
      }
    }
    else
    {
      axom::IndexType tmpStride = 1;
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = tmpStride;
        tmpStride *= m_domainShape[d] + onVertex;
      }
    }

    axom::StackArray<axom::IndexType, DIM> shape;
    for(int d = 0; d < DIM - 1; ++d)
    {
      shape[d] = strides[d + 1] / strides[d];
    }
    shape[DIM - 1] = valuesCount / strides[DIM - 1];

    const T* dataPtr = static_cast<const T*>(valuesNode.data_ptr());
    axom::ArrayView<const T, DIM, MemSpace> rval(dataPtr, shape, strides);

    if(withGhosts == false)
    {
      axom::StackArray<axom::IndexType, DIM> offsets =
        conduitIndicesToMultidimIndices(fieldNode, "offsets", 0);
      axom::StackArray<axom::IndexType, DIM> counts = m_domainShape;
      for(int d = 0; d < DIM; ++d)
      {
        counts[d] += onVertex;
      }
      auto rval1 = rval;
      rval = rval1.subspan(offsets, counts);
    }

    return rval;
  }
  //@}

  //@{
  //!@name Creating data

  /*!
    @brief Create a new scalar nodal data field.
    @param [in] fieldName
    @param [in] association "vertex" or "element"
    @param [in] dtype Conduit data type to put in the field.  Must be at least
                big enough for the strides specified.
    @param [in] strides Data strides.  Set to zero for no ghosts and default strides.
    @param [in] offsets Data index offsets.

    Field data allocation is done by Conduit, so the data lives in
    host memory.  Conduit currently doesn't provide a means to allocate
    the array in device memory.
  */
  void createField(const std::string& fieldName,
                   const std::string& association,
                   const conduit::DataType& dtype,
                   const axom::StackArray<axom::IndexType, DIM>& strides,
                   const axom::StackArray<axom::IndexType, DIM>& offsets)
  {
    if(m_dom->has_path("fields/" + fieldName))
    {
      SLIC_ERROR(
        axom::fmt::format("Cannot create field {}.  It already exists.",
                          fieldName));
    }
    if(association != "vertex" && association != "element")
    {
      SLIC_ERROR(
        axom::fmt::format("Not yet supporting association '{}'.", association));
    }

    auto realExtents = getRealExtents(association);

    conduit::Node& fieldNode = m_dom->fetch("fields/" + fieldName);
    fieldNode["association"] = association;
    fieldNode["topology"] = m_topologyName;

    {
      fieldNode["strides"].set(conduit::DataType::int32(DIM));
      std::int32_t* tmpPtr = fieldNode["strides"].as_int32_ptr();
      for(int d = 0; d < DIM; ++d)
      {
        tmpPtr[d] = strides[d];
      }
    }

    {
      fieldNode["offsets"].set(conduit::DataType::int32(DIM));
      std::int32_t* tmpPtr = fieldNode["offsets"].as_int32_ptr();
      for(int d = 0; d < DIM; ++d)
      {
        tmpPtr[d] = offsets[d];
      }
    }

    axom::IndexType slowDir = 0;
    for(int d = 0; d < DIM; ++d)
    {
      if(strides[slowDir] < strides[d])
      {
        slowDir = d;
      }
    }
    auto extras = dtype.number_of_elements() -
      strides[slowDir] * (offsets[slowDir] + realExtents[slowDir]);
    if(extras < 0)
    {
      SLIC_ERROR(
        axom::fmt::format("Insufficient space allocated for data, given the "
                          "offsets and strides."));
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

  axom::StackArray<axom::IndexType, DIM> m_domainShape;
  axom::StackArray<axom::IndexType, DIM> m_coordsStrides;
  axom::StackArray<axom::IndexType, DIM> m_coordsOffsets;
  axom::StackArray<axom::IndexType, DIM> m_coordsMemShape;

  axom::StackArray<axom::IndexType, DIM> conduitIndicesToMultidimIndices(
    const conduit::Node& node,
    const std::string& path,
    axom::IndexType defaultVal = std::numeric_limits<axom::IndexType>::max()) const
  {
    if(node.has_path(path))
    {
      const conduit::int32* ptr = node.fetch_existing(path).as_int32_ptr();
      return internal::makeStackArray<axom::IndexType, DIM>(ptr);
    }
    return internal::makeStackArray<axom::IndexType, DIM>(defaultVal);
  }

  void computeCoordsDataLayout()
  {
    const conduit::Node& topologyDims =
      m_ctopology->fetch_existing("elements/dims");

    const std::string coordsetName =
      m_ctopology->fetch_existing("coordset").as_string();
    const conduit::Node& coordsValues = m_cdom->fetch_existing(
      axom::fmt::format("coordsets/{}/values", coordsetName));

    const axom::IndexType coordValuesCount =
      coordsValues[0].dtype().number_of_elements();

    const bool coordsAreInterleaved =
      conduit::blueprint::mcarray::is_interleaved(coordsValues);

    for(int i = 0; i < DIM; ++i)
    {
      m_domainShape[i] = topologyDims[i].as_int();
    }

    m_coordsOffsets = conduitIndicesToMultidimIndices(topologyDims, "offsets", 0);

    if(topologyDims.has_child("strides"))
    {
      const conduit::int32* stridesPtr =
        topologyDims.fetch_existing("strides").as_int32_ptr();
      for(int d = 0; d < DIM; ++d)
      {
        m_coordsStrides[d] = stridesPtr[d];
      }
    }
    else
    {
      axom::IndexType tmpStride = coordsAreInterleaved ? DIM : 1;
      for(int d = 0; d < DIM; ++d)
      {
        m_coordsStrides[d] = tmpStride;
        tmpStride *= 1 + m_domainShape[d];
      }
    }

    for(int d = 0; d < DIM - 1; ++d)
    {
      m_coordsMemShape[d] = m_coordsStrides[d + 1] / m_coordsStrides[d];
    }
    m_coordsMemShape[DIM - 1] = coordValuesCount / m_coordsStrides[DIM - 1];
  }
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_MESH_VIEW_UTIL_H_
