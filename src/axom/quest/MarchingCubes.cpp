// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "conduit_blueprint.hpp"
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
  , m_fcnPath()
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
  m_fcnPath = "fields/" + fcnField;
  for(auto& s : m_singles)
  {
    s->setFunctionField(fcnField);
  }
}

void MarchingCubes::computeIsocontour(double contourVal)
{
  SLIC_ASSERT_MSG(!m_fcnPath.empty(),
                  "You must call setFunctionField before computeIsocontour.");

  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];
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
      // TODO: Verify that UnstructuredMesh only supports host memory.
      axom::detail::ArrayOps<axom::IndexType, MemorySpace::Dynamic>::fill(
        domainIdPtr,
        nPrev,
        nNew - nPrev,
        execution_space<axom::SEQ_EXEC>::allocatorID(),
        dId);
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
  , m_fcnPath()
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  SLIC_ASSERT_MSG(
    isValidRuntimePolicy(runtimePolicy),
    fmt::format("Policy '{}' is not a valid runtime policy", runtimePolicy));

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

  m_ndim = conduit::blueprint::mesh::coordset::dims(
    dom.fetch_existing("coordsets/coords"));
  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);

  SLIC_ASSERT_MSG(
    !conduit::blueprint::mcarray::is_interleaved(
      dom.fetch_existing(coordsetPath + "/values")),
    "MarchingCubes currently requires contiguous coordinates layout.");
}

void MarchingCubesSingleDomain::setFunctionField(const std::string& fcnField)
{
  m_fcnPath = "fields/" + fcnField;
  SLIC_ASSERT(m_dom->has_path(m_fcnPath));
  SLIC_ASSERT(m_dom->fetch_existing(m_fcnPath + "/association").as_string() ==
              "vertex");
  SLIC_ASSERT(m_dom->has_path(m_fcnPath + "/values"));
}

void MarchingCubesSingleDomain::computeIsocontour(double contourVal)
{
  SLIC_ASSERT_MSG(!m_fcnPath.empty(),
                  "You must call setFunctionField before computeIsocontour.");

  allocateImpl();
  const std::string coordsetPath = "coordsets/" +
    m_dom->fetch_existing("topologies/" + m_topologyName + "/coordset").as_string();
  m_impl->initialize(*m_dom, coordsetPath, m_fcnPath, m_maskPath);
  m_impl->setContourValue(contourVal);
  m_impl->markCrossings();
  m_impl->scanCrossings();
  m_impl->computeContour();
}

void MarchingCubesSingleDomain::allocateImpl()
{
  using namespace detail::marching_cubes;
  if(m_runtimePolicy == RuntimePolicy::seq)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::SEQ_EXEC, axom::SEQ_EXEC>)
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::SEQ_EXEC, axom::SEQ_EXEC>);
  }
#ifdef _AXOM_MC_USE_OPENMP
  else if(m_runtimePolicy == RuntimePolicy::omp)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::OMP_EXEC, axom::SEQ_EXEC>)
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::OMP_EXEC, axom::SEQ_EXEC>);
  }
#endif
#ifdef _AXOM_MC_USE_CUDA
  else if(m_runtimePolicy == RuntimePolicy::cuda)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>)
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>);
  }
#endif
#ifdef _AXOM_MC_USE_HIP
  else if(m_runtimePolicy == RuntimePolicy::hip)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>)
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>);
  }
#endif
  else
  {
    SLIC_ERROR(axom::fmt::format(
      "MarchingCubesSingleDomain has no implementation for runtime policy {}",
      m_runtimePolicy));
  }
}

}  // end namespace quest
}  // end namespace axom
