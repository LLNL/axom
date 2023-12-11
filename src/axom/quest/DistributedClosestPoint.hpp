// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/quest/detail/DistributedClosestPointImpl.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"

#include <memory>
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
/**
 * \brief Encapsulated the Distributed closest point query for a collection of query points
 * over an "object mesh"
 *
 * The object mesh and the query mesh are provided as conduit nodes
 * using the mesh blueprint schema.  Each of these are distributed
 * over the same MPI rank space.  Ranks are allowed to have any number
 * of domains, including zero.  This class orchestrates passing the
 * query points to all ranks whose object meshes might contain a
 * closest point.
 *
 * \note The class currently supports object meshes that are comprised
 * of a collection of points.  In the future, we'd like to consider
 * more general object meshes, e.g. triangle meshes.
 *
 * To use this class, first set some parameters, such as the runtime
 * execution policy, then pass in the object mesh and build a spatial
 * index over this mesh.  Finally, compute the closest points in the
 * object mesh to each point in a query mesh using the \a
 * computeClosestPoints() function.
 *
 * \note To prevent mixing unrelated MPI communications, you can set a
 * custom MPI Communicator using setMpiCommunicator().
 */
class DistributedClosestPoint
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

public:
  DistributedClosestPoint()
    : m_mpiComm(MPI_COMM_WORLD)
    , m_mpiCommIsPrivate(false)
  {
    setDefaultAllocatorID();
    setMpiCommunicator(MPI_COMM_WORLD);
  }

  ~DistributedClosestPoint()
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

  /*!
    @brief Set runtime execution policy for local queries

    See axom::runtime_policy.
  */
  void setRuntimePolicy(RuntimePolicy policy) { m_runtimePolicy = policy; }

  /*!  @brief Sets the allocator ID to the default associated with the
    execution policy
  */
  void setDefaultAllocatorID()
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

#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
    case RuntimePolicy::cuda:
      defaultAllocatorID =
        axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
      break;
#endif

#ifdef AXOM_RUNTIME_POLICY_USE_HIP
    case RuntimePolicy::hip:
      defaultAllocatorID =
        axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
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

  /*!  @brief Sets the allocator ID to the default associated with the
    execution policy
  */
  void setAllocatorID(int allocatorID)
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

  /**
   * \brief Set the MPI communicator.
   *
   * By default, the communicator is MPI_COMM_WORLD.
   *
   * \param mpiComm The MPI communicator to use.
   * \param duplicate Whether to duplicate mpiComm for exclusive use
   */
  void setMpiCommunicator(MPI_Comm mpiComm, bool duplicate = false)
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

  /**
   * \brief Sets the dimension for the query
   *
   * \note Users do not need to call this function explicitly. The dimension
   * is set by the \a setObjectMesh function
   */
  void setDimension(int dim)
  {
    SLIC_ERROR_IF(
      dim < 2 || dim > 3,
      "DistributedClosestPoint query only supports 2D or 3D queries");
    m_dimension = dim;
  }

  /**
   * \brief Sets the threshold for the query
   *
   * \param [in] threshold Ignore distances greater than this value.
   */
  void setDistanceThreshold(double threshold)
  {
    SLIC_ERROR_IF(threshold < 0.0, "Distance threshold must be non-negative.");
    m_sqDistanceThreshold = threshold * threshold;
  }

  /*!
    @brief Set what fields to output.

    @param [i] field Must be one of these:
      "cp_rank", "cp_index", "cp_distance", "cp_coords", "cp_domain_index".
    @param [i] on Whether to enable outputing @c field.

    By default, all are on.
  */
  void setOutput(const std::string& field, bool on)
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
    SLIC_ERROR_IF(
      f == nullptr,
      axom::fmt::format(
        "Invalid field '{}' should be one of these: "
        "cp_rank, cp_index, cp_distance, cp_coords, cp_domain_index",
        field));
    *f = on;
  }

  /// Sets the logging verbosity of the query. By default the query is not verbose
  void setVerbosity(bool isVerbose) { m_isVerbose = isVerbose; }

  /**
   * \brief Sets the object mesh for the query
   *
   * \param [in] meshNode Conduit node for the object mesh
   * \param [in] topologyName The name of the topology for the object mesh
   *
   * \pre \a meshNode must follow the mesh blueprint convention.
   * \pre Dimension of the mesh must be 2D or 3D
   */
  void setObjectMesh(const conduit::Node& meshNode,
                     const std::string& topologyName)
  {
    SLIC_ASSERT(this->isValidBlueprint(meshNode));

    const bool isMultidomain =
      conduit::blueprint::mesh::is_multi_domain(meshNode);

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

  /**
   * \brief Generates a BVH tree over the object mesh using the runtime execution policy
   *
   * \pre Users must set the object mesh before generating the BVH tree
   * \sa setObjectMesh()
   */
  bool generateBVHTree()
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

  /**
   * \brief Computes the closest point on the object mesh for each point
   * on the provided query mesh
   *
   * \param [in] query_node The root node of a mesh blueprint for the query points,
   * which can be empty if there are no query points for the calling rank
   * \param [in] topology The name of the topology within query_node
   *
   * @c queryMesh should have data on the host, regardless of the runtime
   * policy setting.  Data will be copied to device as needed.
   *
   * On completion, the query mesh contains the following fields:
   *   - cp_rank: will hold the rank of the object point containing the closest point
   *   - cp_domain_index: will hold the index of the object domain containing
   *     the closest points.
   *   - cp_index: Will hold the index of the closest object points.
   *     For multiple object mesh domains on a rank, cp_index is relative to
   *     each domain.
   *   - cp_coords: Will hold the coordinates of the closest points
   *     interleaved in a 1D array.
   *   - cp_distance: Will hold the distances to the closest points.
   * See setOutput() to toggle these outputs.
   *
   * \note The current implementation assumes that the mesh coordinates
   * are interleaved or contiguous.  The output cp_coords will be contiguous.
   */
  void computeClosestPoints(conduit::Node& query_node, const std::string& topology)
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
      m_dcp_2->setOutputSwitches(m_outputRank,
                                 m_outputIndex,
                                 m_outputDistance,
                                 m_outputCoords,
                                 m_outputDomainIndex);
      m_dcp_3->computeClosestPoints(query_node, topology);
      break;
    }
  }

private:
  void allocateQueryInstance()
  {
    SLIC_ASSERT_MSG(m_objectMeshCreated == false, "Object mesh already created");

    switch(m_dimension)
    {
    case 2:
      m_dcp_2 = std::make_unique<internal::DistributedClosestPointImpl<2>>(
        m_runtimePolicy,
        m_allocatorID,
        m_isVerbose);
      m_objectMeshCreated = true;
      break;
    case 3:
      m_dcp_3 = std::make_unique<internal::DistributedClosestPointImpl<3>>(
        m_runtimePolicy,
        m_allocatorID,
        m_isVerbose);
      m_objectMeshCreated = true;
      break;
    }

    SLIC_ASSERT_MSG(
      m_objectMeshCreated,
      "Called allocateQueryInstance, but did not create an instance");
  }

  /// Check validity of blueprint group
  bool isValidBlueprint(const conduit::Node& mesh_node) const
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

  void verifyTopologyName(const conduit::Node& meshNode,
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

private:
  RuntimePolicy m_runtimePolicy {RuntimePolicy::seq};
  MPI_Comm m_mpiComm;
  bool m_mpiCommIsPrivate;
  int m_allocatorID {axom::INVALID_ALLOCATOR_ID};
  int m_dimension {-1};
  bool m_isVerbose {false};
  double m_sqDistanceThreshold {std::numeric_limits<double>::max()};

  bool m_objectMeshCreated {false};

  bool m_outputRank = true;
  bool m_outputIndex = true;
  bool m_outputDistance = true;
  bool m_outputCoords = true;
  bool m_outputDomainIndex = true;

  // One instance per dimension
  std::unique_ptr<internal::DistributedClosestPointImpl<2>> m_dcp_2;
  std::unique_ptr<internal::DistributedClosestPointImpl<3>> m_dcp_3;
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
