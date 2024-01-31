// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
// #include "axom/quest/detail/MarchingCubesHybridParallel.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
MarchingCubes::MarchingCubes(RuntimePolicy runtimePolicy,
                             int allocatorID,
                             MarchingCubesDataParallelism dataParallelism,
                             const conduit::Node& bpMesh,
                             const std::string& topologyName,
                             const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
  , m_allocatorID(allocatorID)
  , m_dataParallelism(dataParallelism)
  , m_singles()
  , m_topologyName(topologyName)
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName(maskField)
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  SLIC_ASSERT_MSG(
    conduit::blueprint::mesh::is_multi_domain(bpMesh),
    "MarchingCubes class input mesh must be in multidomain format.");

  m_singles.reserve(conduit::blueprint::mesh::number_of_domains(bpMesh));
  for(auto& dom : bpMesh.children())
  {
    m_singles.emplace_back(new MarchingCubesSingleDomain(m_runtimePolicy,
                                                         m_allocatorID,
                                                         m_dataParallelism,
                                                         dom,
                                                         m_topologyName,
                                                         maskField));
  }
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
    m_singles[d]->setOutputBuffers(facetNodeIdsView,
                                   facetNodeCoordsView,
                                   facetParentIdsView,
                                   m_facetIndexOffsets[d]);
  }

  for(axom::IndexType d = 0; d < m_singles.size(); ++d)
  {
    m_singles[d]->computeContour();
  }
}

axom::IndexType MarchingCubes::getContourCellCount() const
{
  axom::IndexType contourCellCount = 0;
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    contourCellCount += m_singles[dId]->getContourCellCount();
  }
  return contourCellCount;
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

void MarchingCubes::populateContourMesh(
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
  const std::string& cellIdField,
  const std::string& domainIdField) const
{
  if(!cellIdField.empty() &&
     !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
  {
    // Create cellId field, currently the multidimensional index of the parent cell.
    mesh.createField<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
  }

  if(!domainIdField.empty() &&
     !mesh.hasField(domainIdField, axom::mint::CELL_CENTERED))
  {
    mesh.createField<axom::IndexType>(domainIdField, axom::mint::CELL_CENTERED);
  }

  // Reserve space once for all local domains.
  const axom::IndexType contourCellCount = getContourCellCount();
  const axom::IndexType contourNodeCount = getContourNodeCount();
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
      ? axom::detail::getAllocatorID<axom::MemorySpace::Dynamic>()
      : m_allocatorID;

    if(hostAndInternalMemoriesAreSeparate)
    {
      axom::Array<double, 2> tmpfacetNodeCoords(m_facetNodeCoords,
                                                hostAllocatorId);
      axom::Array<axom::IndexType, 2> tmpfacetNodeIds(m_facetNodeIds,
                                                      hostAllocatorId);
      mesh.appendNodes(tmpfacetNodeCoords.data(),
                       mesh.getDimension() * m_facetCount);
      mesh.appendCells(tmpfacetNodeIds.data(), m_facetCount);
    }
    else
    {
      mesh.appendNodes(m_facetNodeCoords.data(),
                       mesh.getDimension() * m_facetCount);
      mesh.appendCells(m_facetNodeIds.data(), m_facetCount);
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
        mesh.getFieldPtr<axom::IndexType>(domainIdField,
                                          axom::mint::CELL_CENTERED);
      auto tmpContourFacetDomainIds = getContourFacetDomainIds(hostAllocatorId);
      axom::copy(domainIdPtr,
                 tmpContourFacetDomainIds.data(),
                 m_facetCount * sizeof(axom::IndexType));
    }
  }
}

void MarchingCubes::allocateOutputBuffers()
{
  if(!m_singles.empty())
  {
    int ndim = m_singles[0]->spatialDimension();
    const auto nodeCount = m_facetCount * ndim;
    m_facetNodeIds =
      axom::Array<axom::IndexType, 2>({m_facetCount, ndim}, m_allocatorID);
    m_facetNodeCoords = axom::Array<double, 2>({nodeCount, ndim}, m_allocatorID);
    axom::StackArray<axom::IndexType, 1> t1 {m_facetCount};
    m_facetParentIds = axom::Array<axom::IndexType, 1>(t1, m_allocatorID);
  }
}

MarchingCubesSingleDomain::MarchingCubesSingleDomain(
  RuntimePolicy runtimePolicy,
  int allocatorID,
  MarchingCubesDataParallelism dataPar,
  const conduit::Node& dom,
  const std::string& topologyName,
  const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
  , m_allocatorID(allocatorID)
  , m_dataParallelism(dataPar)
  , m_dom(nullptr)
  , m_ndim(0)
  , m_topologyName(topologyName)
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName(maskField)
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  // Set domain first, to get m_ndim, which is required to allocate m_impl.
  setDomain(dom);

  m_impl =
    axom::quest::detail::marching_cubes::newMarchingCubesImpl(m_runtimePolicy,
                                                              m_allocatorID,
                                                              m_ndim);

  m_impl->initialize(*m_dom, m_topologyName, m_maskFieldName);
  m_impl->setDataParallelism(m_dataParallelism);
  return;
}

void MarchingCubesSingleDomain::setDomain(const conduit::Node& dom)
{
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
    SLIC_ASSERT(dom.has_path(m_maskPath + "/values"));
  }

  m_dom = &dom;

  m_ndim = conduit::blueprint::mesh::topology::dims(
    dom.fetch_existing(axom::fmt::format("topologies/{}", m_topologyName)));
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);

  SLIC_ASSERT_MSG(
    !conduit::blueprint::mcarray::is_interleaved(
      dom.fetch_existing(coordsetPath + "/values")),
    "MarchingCubes currently requires contiguous coordinates layout.");
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
