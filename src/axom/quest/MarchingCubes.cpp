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
#include "axom/quest/detail/MarchingCubesHybridParallel.hpp"
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
  AXOM_PERF_MARK_FUNCTION("MarchingCubes::computeIsoContour");
  SLIC_ASSERT_MSG(!m_fcnFieldName.empty(),
                  "You must call setFunctionField before computeIsocontour.");

  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];
    single->setDataParallelism(m_dataParallelism);
    single->computeIsocontour(contourVal);
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

void MarchingCubes::populateContourMesh(
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
  const std::string& cellIdField,
  const std::string& domainIdField)
{
  AXOM_PERF_MARK_FUNCTION("MarchingCubes::populateContourMesh");
  if(!cellIdField.empty() &&
     !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
  {
    mesh.createField<axom::IndexType>(cellIdField,
                                      axom::mint::CELL_CENTERED,
                                      mesh.getDimension());
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
  SLIC_ASSERT(mesh.getNumberOfNodes() == contourNodeCount);
  SLIC_ASSERT(mesh.getNumberOfCells() == contourCellCount);
}

MarchingCubesSingleDomain::MarchingCubesSingleDomain(RuntimePolicy runtimePolicy,
                                                     const conduit::Node& dom,
                                                     const std::string& topologyName,
                                                     const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
  , m_dom(nullptr)
  , m_ndim(0)
  , m_topologyName(topologyName)
  , m_fcnFieldName()
  , m_fcnPath()
  , m_maskFieldName(maskField)
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  setDomain(dom);
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

void MarchingCubesSingleDomain::setFunctionField(const std::string& fcnField)
{
  m_fcnFieldName = fcnField;
  m_fcnPath = "fields/" + fcnField;
  SLIC_ASSERT(m_dom->has_path(m_fcnPath));
  SLIC_ASSERT(m_dom->fetch_existing(m_fcnPath + "/association").as_string() ==
              "vertex");
  SLIC_ASSERT(m_dom->has_path(m_fcnPath + "/values"));
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

}  // end namespace quest
}  // end namespace axom
