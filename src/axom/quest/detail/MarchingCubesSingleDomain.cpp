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
#include "axom/quest/detail/MarchingCubesSingleDomain.hpp"
#include "axom/quest/detail/MarchingCubesImpl.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

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

int32_t MarchingCubesSingleDomain::getDomainId(
  int32_t defaultId) const
{
  int rval = defaultId;
  if(m_dom->has_path("state/domain_id"))
  {
    rval = m_dom->fetch_existing("state/domain_id").to_int32();
  }
  return rval;
}

}  // end namespace quest
}  // end namespace axom
