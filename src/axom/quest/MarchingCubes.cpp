// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubes.cpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
const axom::StackArray<axom::IndexType, 2> twoZeros{0, 0};

MarchingCubes::MarchingCubes(RuntimePolicy runtimePolicy,
                             int allocatorID,
                             MarchingCubesDataParallelism dataParallelism)
  : m_runtimePolicy(runtimePolicy)
  , m_allocatorID(allocatorID)
  , m_dataParallelism(dataParallelism)
  , m_singles()
  , m_topologyName()
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName()
  , m_maskPath()
  , m_facetIndexOffsets(0, 0)
  , m_facetCount(0)
  , m_facetNodeIds(twoZeros, m_allocatorID)
  , m_facetNodeCoords(twoZeros, m_allocatorID)
  , m_facetParentIds(0, 0, m_allocatorID)
{
}

// Set the object up for a blueprint mesh state.
void MarchingCubes::initialize(
  const conduit::Node& bpMesh,
  const std::string& topologyName,
  const std::string& maskField)
{
  SLIC_ASSERT_MSG(
    conduit::blueprint::mesh::is_multi_domain(bpMesh),
    "MarchingCubes class input mesh must be in multidomain format.");

  clear();

  m_topologyName = topologyName;
  m_maskFieldName = maskField;

  /*
    To avoid slow memory allocations (especially on GPUs) keep the
    single-domain objects around and just re-initialize them.  Arrays
    will be cleared, but not deallocated.  To really deallocate
    memory, deallocate the MarchingCubes object.  The actual number
    of domains is m_domainCount, not m_singles.size().
  */
  auto newDomainCount = conduit::blueprint::mesh::number_of_domains(bpMesh);

  while( m_singles.size() < newDomainCount )
  {
    m_singles.emplace_back(new MarchingCubesSingleDomain(m_runtimePolicy,
                                                         m_allocatorID,
                                                         m_dataParallelism));
  }

  for (int d = 0; d < newDomainCount; ++d)
  {
    const auto& dom = bpMesh.child(d);
    m_singles[d]->initialize(dom, m_topologyName, maskField);
  }

  m_domainCount = newDomainCount;
}

void MarchingCubes::setFunctionField(const std::string& fcnField)
{
  m_fcnFieldName = fcnField;
  m_fcnPath = "fields/" + fcnField;
  for(auto& s : m_singles)
  {
    s->setFunctionField(fcnField);
  }
}

void MarchingCubes::computeIsocontour(double contourVal)
{
  AXOM_PERF_MARK_FUNCTION("MarchingCubes::computeIsoContour");
  // Mark and scan domains while adding up their
  // facet counts to get the total facet counts.
  m_facetIndexOffsets.resize(m_singles.size());
  m_facetCount = 0;
  for(axom::IndexType d = 0; d < m_singles.size(); ++d)
  {
    auto& single = *m_singles[d];
    single.setContourValue(contourVal);
    single.markCrossings();
    single.scanCrossings();
    m_facetIndexOffsets[d] = m_facetCount;
    m_facetCount += single.getContourCellCount();
  }

  allocateOutputBuffers();

  // Tell singles where to put contour data.
  auto facetNodeIdsView = m_facetNodeIds.view();
  auto facetNodeCoordsView = m_facetNodeCoords.view();
  auto facetParentIdsView = m_facetParentIds.view();
  for(axom::IndexType d = 0; d < m_singles.size(); ++d)
  {
    m_singles[d]->getImpl().setOutputBuffers(facetNodeIdsView,
                                             facetNodeCoordsView,
                                             facetParentIdsView,
                                             m_facetIndexOffsets[d]);
  }

  for(axom::IndexType d = 0; d < m_singles.size(); ++d)
  {
    m_singles[d]->computeContour();
  }
}

axom::IndexType MarchingCubes::getContourNodeCount() const
{
  axom::IndexType contourNodeCount = 0;
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    contourNodeCount += m_singles[dId]->getContourNodeCount();
  }
  return contourNodeCount;
}

/*
  Domain ids are provided as a new Array instead of ArrayView because
  we don't store it internally.
*/
axom::Array<MarchingCubes::DomainIdType> MarchingCubes::getContourFacetDomainIds(
  int allocatorID) const
{
  // Put parent domain ids into a new Array.
  const axom::IndexType len = getContourCellCount();
  axom::Array<MarchingCubes::DomainIdType> rval(
    len,
    len,
    allocatorID != axom::INVALID_ALLOCATOR_ID ? allocatorID : m_allocatorID);
  for(int d = 0; d < m_singles.size(); ++d)
  {
    MarchingCubes::DomainIdType domainId = m_singles[d]->getDomainId(d);
    axom::IndexType contourCellCount = m_singles[d]->getContourCellCount();
    axom::IndexType offset = m_facetIndexOffsets[d];
    axom::detail::ArrayOps<MarchingCubes::DomainIdType, MemorySpace::Dynamic>::fill(
      rval.data(),
      offset,
      contourCellCount,
      allocatorID,
      domainId);
  }
  return rval;
}

void MarchingCubes::clear()
{
  for(int d = 0; d < m_singles.size(); ++d)
  {
    m_singles[d]->getImpl().clear();
  }
  m_domainCount = 0;
  m_facetNodeIds.clear();
  m_facetNodeCoords.clear();
  m_facetParentIds.clear();
}

void MarchingCubes::populateContourMesh(
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
  const std::string& cellIdField,
  const std::string& domainIdField) const
{
  AXOM_PERF_MARK_FUNCTION("MarchingCubes::populateContourMesh");
  if(!cellIdField.empty() &&
     !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
  {
    // Create cellId field, currently the multidimensional index of the parent cell.
    mesh.createField<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
  }

  if(!domainIdField.empty() &&
     !mesh.hasField(domainIdField, axom::mint::CELL_CENTERED))
  {
    mesh.createField<DomainIdType>(domainIdField, axom::mint::CELL_CENTERED);
  }

  // Reserve space once for all local domains.
  const axom::IndexType contourCellCount = getContourCellCount();
  const axom::IndexType contourNodeCount = getContourNodeCount();
  // Temporarily disable reservation due to unknown bug.
  // See https://github.com/LLNL/axom/pull/1271
  mesh.reserveCells(contourCellCount);
  mesh.reserveNodes(contourNodeCount);

  if(m_facetCount)
  {
    // Put nodes and cells into the mesh.
    // If data is not in host memory, copy to temporary host memory first.
    axom::MemorySpace internalMemorySpace =
      axom::detail::getAllocatorSpace(m_allocatorID);
    const bool hostAndInternalMemoriesAreSeparate =
      internalMemorySpace != axom::MemorySpace::Dynamic
#ifdef AXOM_USE_UMPIRE
      && internalMemorySpace != axom::MemorySpace::Host
#endif
      ;
    const int hostAllocatorId = hostAndInternalMemoriesAreSeparate
      ?
#ifdef AXOM_USE_UMPIRE
      axom::detail::getAllocatorID<axom::MemorySpace::Host>()
#else
      axom::detail::getAllocatorID<axom::MemorySpace::Dynamic>()
#endif
      : m_allocatorID;

    if(hostAndInternalMemoriesAreSeparate)
    {
      axom::Array<double, 2> tmpfacetNodeCoords(m_facetNodeCoords,
                                                hostAllocatorId);
      axom::Array<axom::IndexType, 2> tmpfacetNodeIds(m_facetNodeIds,
                                                      hostAllocatorId);
      mesh.appendNodes(tmpfacetNodeCoords.data(), contourNodeCount);
      mesh.appendCells(tmpfacetNodeIds.data(), contourCellCount);
    }
    else
    {
      mesh.appendNodes(m_facetNodeCoords.data(), contourNodeCount);
      mesh.appendCells(m_facetNodeIds.data(), contourCellCount);
    }

    if(!cellIdField.empty())
    {
      // Put parent cell ids into the mesh.
      axom::IndexType* cellIdPtr =
        mesh.getFieldPtr<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
      axom::copy(cellIdPtr,
                 m_facetParentIds.data(),
                 m_facetCount * sizeof(axom::IndexType));
    }

    if(!domainIdField.empty())
    {
      // Put parent domain ids into the mesh.
      auto* domainIdPtr =
        mesh.getFieldPtr<DomainIdType>(domainIdField, axom::mint::CELL_CENTERED);
      auto tmpContourFacetDomainIds = getContourFacetDomainIds(hostAllocatorId);
      axom::copy(domainIdPtr,
                 tmpContourFacetDomainIds.data(),
                 m_facetCount * sizeof(DomainIdType));
    }
  }
}

void MarchingCubes::allocateOutputBuffers()
{
  AXOM_PERF_MARK_FUNCTION("MarchingCubes::allocateOutputBuffers");
  if(!m_singles.empty())
  {
    int ndim = m_singles[0]->spatialDimension();
    const auto nodeCount = m_facetCount * ndim;
    m_facetNodeIds.resize(axom::StackArray<axom::IndexType, 2>{m_facetCount, ndim}, 0);
    m_facetNodeCoords.resize(axom::StackArray<axom::IndexType, 2>{nodeCount, ndim}, 0.0);
    m_facetParentIds.resize(axom::StackArray<axom::IndexType, 1>{m_facetCount}, 0);
  }
}

MarchingCubesSingleDomain::MarchingCubesSingleDomain(
  RuntimePolicy runtimePolicy,
  int allocatorID,
  MarchingCubesDataParallelism dataPar)
  : m_runtimePolicy(runtimePolicy)
  , m_allocatorID(allocatorID)
  , m_dataParallelism(dataPar)
  , m_dom(nullptr)
  , m_ndim(0)
  , m_topologyName()
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName()
  , m_maskPath()
{
  return;
}

void MarchingCubesSingleDomain::initialize(
  const conduit::Node& dom,
  const std::string& topologyName,
  const std::string& maskField)
{
  m_topologyName = topologyName;

  SLIC_ASSERT_MSG(
    !conduit::blueprint::mesh::is_multi_domain(dom),
    "MarchingCubesSingleDomain is single-domain only.  Try MarchingCubes.");
  SLIC_ASSERT(
    dom.fetch_existing("topologies/" + m_topologyName + "/type").as_string() ==
    "structured");

  const std::string coordsetPath = "coordsets/" +
    dom.fetch_existing("topologies/" + m_topologyName + "/coordset").as_string();
  SLIC_ASSERT(dom.has_path(coordsetPath));

  if(!m_maskPath.empty())
  {
    m_maskPath = maskField.empty() ? std::string() : "fields/" + maskField;
    SLIC_ASSERT(dom.has_path(m_maskPath + "/values"));
  }
  else
  {
    m_maskPath.clear();
  }

  m_dom = &dom;

  m_ndim = conduit::blueprint::mesh::topology::dims(
    dom.fetch_existing(axom::fmt::format("topologies/{}", m_topologyName)));
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);

  SLIC_ASSERT_MSG(
    !conduit::blueprint::mcarray::is_interleaved(
      dom.fetch_existing(coordsetPath + "/values")),
    "MarchingCubes currently requires contiguous coordinates layout.");

  m_impl =
    axom::quest::detail::marching_cubes::newMarchingCubesImpl(m_runtimePolicy,
                                                              m_allocatorID,
                                                              m_ndim);

  m_impl->initialize(dom, topologyName, maskField);
  m_impl->setDataParallelism(m_dataParallelism);
}

MarchingCubes::DomainIdType MarchingCubesSingleDomain::getDomainId(
  MarchingCubes::DomainIdType defaultId) const
{
  MarchingCubes::DomainIdType rval = defaultId;
  if(m_dom->has_path("state/domain_id"))
  {
    rval = m_dom->fetch_existing("state/domain_id").as_int();
  }
  return rval;
}

}  // end namespace quest
}  // end namespace axom
