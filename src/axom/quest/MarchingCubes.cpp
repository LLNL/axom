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
#include "axom/quest/detail/MarchingCubesFullParallel.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
MarchingCubes::MarchingCubes(RuntimePolicy runtimePolicy,
                             const conduit::Node& bpMesh,
                             const std::string& topologyName,
                             const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
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
#if 1
  // Mark and scan domains while adding up their
  // facet counts to get the total facet counts.
  m_facetIndexOffsets.resize(m_singles.size());
  m_facetCount = 0;
  for(axom::IndexType d=0; d<m_singles.size(); ++d)
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
  for(axom::IndexType d=0; d<m_singles.size(); ++d)
  {
    m_singles[d]->setOutputBuffers(facetNodeIdsView,
                                   facetNodeCoordsView,
                                   facetParentIdsView,
                                   m_facetIndexOffsets[d]);
  }

  for(axom::IndexType d=0; d<m_singles.size(); ++d)
  {
    m_singles[d]->computeContour();
  }
#else
  SLIC_ASSERT_MSG(!m_fcnFieldName.empty(),
                  "You must call setFunctionField before computeIsocontour.");

  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];
    single->setDataParallelism(m_dataParallelism);
    single->computeIsocontour(contourVal);
  }
#endif
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

void MarchingCubes::populateContourMesh(
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
  const std::string& cellIdField,
  const std::string& domainIdField)
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

#if 1
  if (m_facetCount) {
    mesh.appendCells(m_facetNodeIds.data(), m_facetCount);
    mesh.appendNodes(m_facetNodeCoords.data(), mesh.getDimension()*m_facetCount);

    axom::IndexType* cellIdPtr =
      mesh.getFieldPtr<axom::IndexType>(cellIdField,
                                        axom::mint::CELL_CENTERED);
    axom::copy(cellIdPtr,
               m_facetParentIds.data(),
               m_facetCount*sizeof(axom::IndexType));

    // TODO: Move domain id stuff into a separate function.
    auto* domainIdPtr =
      mesh.getFieldPtr<axom::IndexType>(domainIdField,
                                        axom::mint::CELL_CENTERED);
    for(int d = 0; d < m_singles.size(); ++d)
    {
      axom::detail::ArrayOps<axom::IndexType, MemorySpace::Dynamic>::fill(
        domainIdPtr,
        m_facetIndexOffsets[d],
        m_singles[d]->getContourCellCount(),
        execution_space<axom::SEQ_EXEC>::allocatorID(),
        m_singles[d]->getDomainId(d));
    }
  }
#else
  // Populate mesh from single domains and add domain id if requested.
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];

    auto nPrev = mesh.getNumberOfCells();
    single->populateContourMesh(mesh, cellIdField);
    auto nNew = mesh.getNumberOfCells();

    if(nNew > nPrev && !domainIdField.empty())
    {
      auto* domainIdPtr =
        mesh.getFieldPtr<axom::IndexType>(domainIdField,
                                          axom::mint::CELL_CENTERED);

      int userDomainId = single->getDomainId(dId);

      axom::detail::ArrayOps<axom::IndexType, MemorySpace::Dynamic>::fill(
        domainIdPtr,
        nPrev,
        nNew - nPrev,
        execution_space<axom::SEQ_EXEC>::allocatorID(),
        userDomainId);
    }
  }
#endif
  SLIC_ASSERT(mesh.getNumberOfNodes() == contourNodeCount);
  SLIC_ASSERT(mesh.getNumberOfCells() == contourCellCount);
}

void MarchingCubes::allocateOutputBuffers()
{
  int allocatorId = -1;
  if(m_runtimePolicy == MarchingCubes::RuntimePolicy::seq)
  {
    allocatorId = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::omp)
  {
    allocatorId = axom::execution_space<axom::OMP_EXEC>::allocatorID();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::cuda)
  {
    allocatorId = axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(m_runtimePolicy == MarchingCubes::RuntimePolicy::hip)
  {
    allocatorId = axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
  }
#endif
  else
  {
    SLIC_ERROR(axom::fmt::format(
      "MarchingCubes doesn't recognize runtime policy {}",
      m_runtimePolicy));
  }

  if (!m_singles.empty())
  {
    int ndim = m_singles[0]->spatialDimension();
    const auto nodeCount = m_facetCount * ndim;
    m_facetNodeIds =
      axom::Array<axom::IndexType, 2>({m_facetCount, ndim}, allocatorId);
    m_facetNodeCoords =
      axom::Array<double, 2>({nodeCount, ndim},  allocatorId);
    axom::StackArray<axom::IndexType, 1> t1{m_facetCount};
    m_facetParentIds =
      axom::Array<axom::IndexType, 1>(t1, allocatorId);
  }
}

MarchingCubesSingleDomain::MarchingCubesSingleDomain(RuntimePolicy runtimePolicy,
                                                     MarchingCubesDataParallelism dataPar,
                                                     const conduit::Node& dom,
                                                     const std::string& topologyName,
                                                     const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
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

  /*
    We have 2 implementations.  MarchingCubesHybridParallel is faster on the host
    and MarchingCubesFullParallel is faster on GPUs.  Both work in all cases.
    We can choose based on runtime policy or by user choice
  */
  if(m_dataParallelism ==
     axom::quest::MarchingCubesDataParallelism::hybridParallel ||
     (m_dataParallelism == axom::quest::MarchingCubesDataParallelism::byPolicy &&
      m_runtimePolicy == RuntimePolicy::seq))
  {
    SLIC_WARNING("Not really using hybrid while developing.  Using full parallel.");
    m_impl = axom::quest::detail::marching_cubes::newMarchingCubesFullParallel(
      m_runtimePolicy,
      m_ndim);
  }
  else
  {
    m_impl = axom::quest::detail::marching_cubes::newMarchingCubesFullParallel(
      m_runtimePolicy,
      m_ndim);
  }

  m_impl->initialize(*m_dom, m_topologyName, m_maskFieldName);
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

int MarchingCubesSingleDomain::getDomainId(int defaultId) const
{
  int rval = defaultId;
  if(m_dom->has_path("state/domain_id"))
  {
    rval = m_dom->fetch_existing("state/domain_id").as_int();
  }
  return rval;
}

#if 0
void MarchingCubesSingleDomain::computeIsocontour(double contourVal)
{
  SLIC_ASSERT_MSG(!m_fcnFieldName.empty(),
                  "You must call setFunctionField before computeIsocontour.");

  // We have 2 implementations.  MarchingCubesHybridParallel is faster on the host
  // and MarchingCubesFullParallel is faster on GPUs.  Both work in all cases.
  // We can choose based on runtime policy or by user choice
  if(m_dataParallelism ==
       axom::quest::MarchingCubesDataParallelism::hybridParallel ||
     (m_dataParallelism == axom::quest::MarchingCubesDataParallelism::byPolicy &&
      m_runtimePolicy == RuntimePolicy::seq))
  {
    m_impl = axom::quest::detail::marching_cubes::newMarchingCubesHybridParallel(
      m_runtimePolicy,
      m_ndim);
  }
  else
  {
    m_impl = axom::quest::detail::marching_cubes::newMarchingCubesFullParallel(
      m_runtimePolicy,
      m_ndim);
  }
  m_impl->initialize(*m_dom, m_topologyName, m_fcnFieldName, m_maskFieldName);
  m_impl->setContourValue(contourVal);
  m_impl->computeContourMesh();
}
#endif

}  // end namespace quest
}  // end namespace axom
