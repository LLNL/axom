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
                             const std::string& coordsetName,
                             const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
  , m_singles()
  , m_coordsetPath("coordsets/" + coordsetName)
  , m_fcnPath()
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  m_singles.reserve(conduit::blueprint::mesh::number_of_domains(bpMesh));
  for(auto& dom : bpMesh.children())
  {
    m_singles.emplace_back(
      new MarchingCubesSingleDomain(m_runtimePolicy, dom, coordsetName, maskField));
  }
}

void MarchingCubes::set_function_field(const std::string& fcnField)
{
  m_fcnPath = "fields/" + fcnField;
  for(auto& s : m_singles)
  {
    s->set_function_field(fcnField);
  }
}

void MarchingCubes::compute_isocontour(double contourVal)
{
  SLIC_ASSERT_MSG(
    !m_fcnPath.empty(),
    "You must call set_function_field before compute_isocontour.");

  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];
    single->compute_isocontour(contourVal);
  }
}

axom::IndexType MarchingCubes::get_contour_cell_count() const
{
  axom::IndexType contourCellCount = 0;
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    contourCellCount += m_singles[dId]->get_contour_cell_count();
  }
  return contourCellCount;
}

axom::IndexType MarchingCubes::get_contour_node_count() const
{
  axom::IndexType contourNodeCount = 0;
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    contourNodeCount += m_singles[dId]->get_contour_node_count();
  }
  return contourNodeCount;
}

void MarchingCubes::populate_contour_mesh(
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
  const axom::IndexType contourCellCount = get_contour_cell_count();
  const axom::IndexType contourNodeCount = get_contour_node_count();
  mesh.reserveCells(contourCellCount);
  mesh.reserveNodes(contourNodeCount);

  // Populate mesh from single domains and add domain id if requested.
  for(int dId = 0; dId < m_singles.size(); ++dId)
  {
    std::unique_ptr<MarchingCubesSingleDomain>& single = m_singles[dId];

    auto nPrev = mesh.getNumberOfCells();
    single->populate_contour_mesh(mesh, cellIdField);
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
                                                     const std::string& coordsetName,
                                                     const std::string& maskField)
  : m_runtimePolicy(runtimePolicy)
  , m_dom(nullptr)
  , m_ndim(0)
  , m_coordsetPath("coordsets/" + coordsetName)
  , m_fcnPath()
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
{
  SLIC_ASSERT_MSG(
    isValidRuntimePolicy(runtimePolicy),
    fmt::format("Policy '{}' is not a valid runtime policy", runtimePolicy));

  set_domain(dom);
  return;
}

void MarchingCubesSingleDomain::set_domain(const conduit::Node& dom)
{
  SLIC_ASSERT_MSG(
    !conduit::blueprint::mesh::is_multi_domain(dom),
    "MarchingCubesSingleDomain is single-domain only.  Try MarchingCubes.");

  SLIC_ASSERT(dom.has_path(m_coordsetPath));
  SLIC_ASSERT(dom["topologies/mesh/type"].as_string() == "structured");

  if(!m_maskPath.empty())
  {
    SLIC_ASSERT(dom.has_path(m_maskPath + "/values"));
  }

  m_dom = &dom;

  const conduit::Node& dimsNode =
    m_dom->fetch_existing("topologies/mesh/elements/dims");

  m_ndim = dimsNode.number_of_children();

  SLIC_ASSERT(m_ndim >= 2 && m_ndim <= 3);

  const conduit::Node& coordsValues = dom[m_coordsetPath + "/values"];
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);
  SLIC_ASSERT_MSG(
    !isInterleaved,
    "MarchingCubes currently requires contiguous coordinates layout.");
}

void MarchingCubesSingleDomain::set_function_field(const std::string& fcnField)
{
  m_fcnPath = "fields/" + fcnField;
  SLIC_ASSERT(m_dom->has_path(m_fcnPath));
  SLIC_ASSERT(m_dom->fetch_existing(m_fcnPath + "/association").as_string() ==
              "vertex");
  SLIC_ASSERT(m_dom->has_path(m_fcnPath + "/values"));
}

void MarchingCubesSingleDomain::compute_isocontour(double contourVal)
{
  SLIC_ASSERT_MSG(
    !m_fcnPath.empty(),
    "You must call set_function_field before compute_isocontour.");

  allocate_impl();
  m_impl->initialize(*m_dom, m_coordsetPath, m_fcnPath, m_maskPath);
  m_impl->set_contour_value(contourVal);
  m_impl->mark_crossings();
  m_impl->scan_crossings();
  m_impl->compute_contour();
}

void MarchingCubesSingleDomain::allocate_impl()
{
  using namespace detail::marching_cubes;
  if(m_runtimePolicy == RuntimePolicy::seq)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(new MarchingCubesImpl<2, axom::SEQ_EXEC>)
      : std::unique_ptr<ImplBase>(new MarchingCubesImpl<3, axom::SEQ_EXEC>);
  }
#ifdef _AXOM_MC_USE_OPENMP
  else if(m_runtimePolicy == RuntimePolicy::omp)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(new MarchingCubesImpl<2, axom::OMP_EXEC>)
      : std::unique_ptr<ImplBase>(new MarchingCubesImpl<3, axom::OMP_EXEC>);
  }
#endif
#ifdef _AXOM_MC_USE_CUDA
  else if(m_runtimePolicy == RuntimePolicy::cuda)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(new MarchingCubesImpl<2, axom::CUDA_EXEC<256>>)
      : std::unique_ptr<ImplBase>(new MarchingCubesImpl<3, axom::CUDA_EXEC<256>>);
  }
#endif
#ifdef _AXOM_MC_USE_HIP
  else if(m_runtimePolicy == RuntimePolicy::hip)
  {
    m_impl = m_ndim == 2
      ? std::unique_ptr<ImplBase>(new MarchingCubesImpl<2, axom::HIP_EXEC<256>>)
      : std::unique_ptr<ImplBase>(new MarchingCubesImpl<3, axom::HIP_EXEC<256>>);
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
