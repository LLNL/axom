// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/quest/DistributedClosestPoint.hpp"
#include "axom/quest/detail/DistributedClosestPointImpl.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"

#include <limits>
#include <cstdlib>

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

namespace axom
{
namespace quest
{
DistributedClosestPoint::DistributedClosestPoint()
  : m_mpiComm(MPI_COMM_WORLD)
  , m_mpiCommIsPrivate(false)
  , m_allocatorID(axom::INVALID_ALLOCATOR_ID)
  , m_sqDistanceThreshold(std::numeric_limits<double>::max())
{
  setDefaultAllocatorID();
  setMpiCommunicator(MPI_COMM_WORLD);
}

DistributedClosestPoint::~DistributedClosestPoint()
{
  if(m_mpiCommIsPrivate)
  {
    int mpiIsFinalized = 0;
    MPI_Finalized(&mpiIsFinalized);
    if(!mpiIsFinalized)
    {
      MPI_Comm_free(&m_mpiComm);
    }
  }
}

void DistributedClosestPoint::setDefaultAllocatorID()
{
  int defaultAllocatorID = axom::INVALID_ALLOCATOR_ID;
  switch(m_runtimePolicy)
  {
  case RuntimePolicy::seq:
    defaultAllocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    break;

#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  case RuntimePolicy::omp:
    defaultAllocatorID = axom::execution_space<axom::OMP_EXEC>::allocatorID();
    break;
#endif

#ifdef __CUDACC__
  #ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  case RuntimePolicy::cuda:

    // Use unified memory
    defaultAllocatorID =
      axom::getUmpireResourceAllocatorID(umpire::resource::Unified);
    break;
  #endif
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  case RuntimePolicy::hip:

    // Use unified memory
    defaultAllocatorID =
      axom::getUmpireResourceAllocatorID(umpire::resource::Unified);

    break;
#endif
  }
  if(defaultAllocatorID == axom::INVALID_ALLOCATOR_ID)
  {
    SLIC_ERROR(
      axom::fmt::format("There is no default allocator for runtime policy {}",
                        m_runtimePolicy));
  }
  setAllocatorID(defaultAllocatorID);
}

void DistributedClosestPoint::setAllocatorID(int allocatorID)
{
  SLIC_ASSERT_MSG(allocatorID != axom::INVALID_ALLOCATOR_ID,
                  "Invalid allocator id.");
  m_allocatorID = allocatorID;

  if(m_dcp_2 != nullptr)
  {
    m_dcp_2->setAllocatorID(m_allocatorID);
  }
  if(m_dcp_3 != nullptr)
  {
    m_dcp_3->setAllocatorID(m_allocatorID);
  }
}

void DistributedClosestPoint::setMpiCommunicator(MPI_Comm mpiComm, bool duplicate)
{
  if(m_mpiCommIsPrivate)
  {
    MPI_Comm_free(&m_mpiComm);
  }

  if(duplicate)
  {
    MPI_Comm_dup(mpiComm, &m_mpiComm);
  }
  else
  {
    m_mpiComm = mpiComm;
  }
  m_mpiCommIsPrivate = duplicate;
}

void DistributedClosestPoint::setDimension(int dim)
{
  SLIC_ERROR_IF(dim < 2 || dim > 3,
                "DistributedClosestPoint query only supports 2D or 3D queries");
  m_dimension = dim;
}

void DistributedClosestPoint::setDistanceThreshold(double threshold)
{
  SLIC_ERROR_IF(threshold < 0.0, "Distance threshold must be non-negative.");
  m_sqDistanceThreshold = threshold * threshold;
}

void DistributedClosestPoint::setOutput(const std::string& field, bool on)
{
  // clang-format off
  bool* f
    = field == "cp_rank" ? &m_outputRank
    : field == "cp_index" ? &m_outputIndex
    : field == "cp_distance" ? &m_outputDistance
    : field == "cp_coords" ? &m_outputCoords
    : field == "cp_domain_index" ? &m_outputDomainIndex
    : nullptr;
  // clang-format on
  SLIC_ERROR_IF(f == nullptr,
                axom::fmt::format(
                  "Invalid field '{}' should be one of these: "
                  "cp_rank, cp_index, cp_distance, cp_coords, cp_domain_index",
                  field));
  *f = on;
}

void DistributedClosestPoint::setObjectMesh(const conduit::Node& meshNode,
                                            const std::string& topologyName)
{
  SLIC_ASSERT(this->isValidBlueprint(meshNode));

  const bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(meshNode);

  // If meshNode isn't multidomain, create a temporary multidomain representation.
  std::shared_ptr<conduit::Node> tmpNode;
  if(!isMultidomain)
  {
    tmpNode = std::make_shared<conduit::Node>();
    conduit::blueprint::mesh::to_multi_domain(meshNode, *tmpNode);
  }
  const conduit::Node& mdMeshNode(isMultidomain ? meshNode : *tmpNode);
  verifyTopologyName(mdMeshNode, topologyName);

  auto domainCount = conduit::blueprint::mesh::number_of_domains(mdMeshNode);

  // Extract the dimension from the coordinate values group
  // use allreduce since some ranks might be empty
  {
    int localDim = -1;
    if(domainCount > 0)
    {
      const conduit::Node& domain0 = mdMeshNode.child(0);
      const conduit::Node& topology =
        domain0.fetch_existing("topologies/" + topologyName);
      const std::string coordsetName =
        topology.fetch_existing("coordset").as_string();
      const conduit::Node& coordset =
        domain0.fetch_existing("coordsets/" + coordsetName);
      const conduit::Node& coordsetValues = coordset.fetch_existing("values");
      localDim = internal::extractDimension(coordsetValues);
    }
    int dim = -1;
    MPI_Allreduce(&localDim, &dim, 1, MPI_INT, MPI_MAX, m_mpiComm);
    setDimension(dim);
  }

  allocateQueryInstance();

  switch(m_dimension)
  {
  case 2:
    m_dcp_2->importObjectPoints(mdMeshNode, topologyName);
    break;
  case 3:
    m_dcp_3->importObjectPoints(mdMeshNode, topologyName);
    break;
  }

  return;
}

bool DistributedClosestPoint::generateBVHTree()
{
  SLIC_ASSERT_MSG(m_objectMeshCreated,
                  "Must call 'setObjectMesh' before calling generateBVHTree");

  bool success = false;

  // dispatch to implementation class over dimension
  switch(m_dimension)
  {
  case 2:
    success = m_dcp_2->generateBVHTree();
    break;
  case 3:
    success = m_dcp_3->generateBVHTree();
    break;
  }

  return success;
}

void DistributedClosestPoint::computeClosestPoints(conduit::Node& query_node,
                                                   const std::string& topology)
{
  SLIC_ASSERT_MSG(m_objectMeshCreated,
                  "Must call 'setObjectMesh' before calling generateBVHTree");

  SLIC_ASSERT(this->isValidBlueprint(query_node));

  // dispatch to implementation class over dimension
  switch(m_dimension)
  {
  case 2:
    m_dcp_2->setSquaredDistanceThreshold(m_sqDistanceThreshold);
    m_dcp_2->setMpiCommunicator(m_mpiComm);
    m_dcp_2->setOutputSwitches(m_outputRank,
                               m_outputIndex,
                               m_outputDistance,
                               m_outputCoords,
                               m_outputDomainIndex);
    m_dcp_2->computeClosestPoints(query_node, topology);
    break;
  case 3:
    m_dcp_3->setSquaredDistanceThreshold(m_sqDistanceThreshold);
    m_dcp_3->setMpiCommunicator(m_mpiComm);
    m_dcp_3->setOutputSwitches(m_outputRank,
                               m_outputIndex,
                               m_outputDistance,
                               m_outputCoords,
                               m_outputDomainIndex);
    m_dcp_3->computeClosestPoints(query_node, topology);
    break;
  }
}

void DistributedClosestPoint::allocateQueryInstance()
{
  SLIC_ASSERT_MSG(m_objectMeshCreated == false, "Object mesh already created");

  switch(m_dimension)
  {
  case 2:
    m_dcp_2 =
      std::make_unique<internal::DistributedClosestPointImpl<2>>(m_runtimePolicy,
                                                                 m_allocatorID,
                                                                 m_isVerbose);
    m_objectMeshCreated = true;
    break;
  case 3:
    m_dcp_3 =
      std::make_unique<internal::DistributedClosestPointImpl<3>>(m_runtimePolicy,
                                                                 m_allocatorID,
                                                                 m_isVerbose);
    m_objectMeshCreated = true;
    break;
  }

  SLIC_ASSERT_MSG(
    m_objectMeshCreated,
    "Called allocateQueryInstance, but did not create an instance");
}

bool DistributedClosestPoint::isValidBlueprint(const conduit::Node& mesh_node) const
{
  bool success = true;
  conduit::Node info;
  if(!conduit::blueprint::mpi::verify("mesh", mesh_node, info, m_mpiComm))
  {
    SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
    success = false;
  }

  return success;
}

void DistributedClosestPoint::verifyTopologyName(const conduit::Node& meshNode,
                                                 const std::string& topologyName)
{
  std::string coordsetPath;
  const std::string topologyPath =
    axom::fmt::format("topologies/{}", topologyName);
  for(axom::IndexType d = 0; d < meshNode.number_of_children(); ++d)
  {
    const auto& domain = meshNode.child(d);
    if(!domain.has_path(topologyPath))
    {
      auto errMsg = fmt::format("No such topology '{}' found.", topologyName);
      if(domain.has_path("coordsets/" + topologyName))
      {
        errMsg += fmt::format(
          "  You may have mistakenly specified a coordset name."
          "  The interface has changed to use topology name"
          " instead of coordset.");
      }
      SLIC_ERROR(errMsg);
    }
  }
}

}  // end namespace quest
}  // end namespace axom
