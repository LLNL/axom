// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISTRIBUTED_CLOSEST_POINT_H_
#define QUEST_DISTRIBUTED_CLOSEST_POINT_H_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mcarray.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"
#include "conduit_relay_io.hpp"

#include <memory>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <vector>

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

// Add some helper preprocessor defines for using OPENMP, CUDA, and HIP policies
// within the distributed closest point query.
// These are only used when building with RAJA and Umpire
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #ifdef AXOM_USE_OPENMP
    #define _AXOM_DCP_USE_OPENMP
  #endif
  #ifdef AXOM_USE_CUDA
    #define _AXOM_DCP_USE_CUDA
  #endif
  #ifdef AXOM_USE_HIP
    #define _AXOM_DCP_USE_HIP
  #endif
#endif

namespace axom
{
namespace quest
{
/// Enum for runtime execution policy
enum class DistributedClosestPointRuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2,
  hip = 3
};

namespace internal
{
// Utility function to dump a conduit node on each rank, e.g. for debugging
inline void dump_node(const conduit::Node& n,
                      const std::string&& fname,
                      const std::string& protocol = "json")
{
  conduit::relay::io::save(n, fname, protocol);
};

/**
 * \brief Utility function to get a typed pointer to the beginning of an array
 * stored by a conduit::Node
 */
template <typename T>
T* getPointer(conduit::Node& node)
{
  T* ptr = node.value();
  return ptr;
}

/**
 * \brief Utility function to create an axom::ArrayView over the array
 * of native types stored by a conduit::Node
 */
template <typename T>
axom::ArrayView<T> ArrayView_from_Node(conduit::Node& node, int sz)
{
  T* ptr = node.value();
  return axom::ArrayView<T>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,2>
 *
 * \warning Assumes the underlying data is an MCArray with stride 2 access
 */
template <>
inline axom::ArrayView<primal::Point<double, 2>> ArrayView_from_Node(
  conduit::Node& node,
  int sz)
{
  using PointType = primal::Point<double, 2>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,3>
 *
 * \warning Assumes the underlying data is an MCArray with stride 3 access
 */
template <>
inline axom::ArrayView<primal::Point<double, 3>> ArrayView_from_Node(
  conduit::Node& node,
  int sz)
{
  using PointType = primal::Point<double, 3>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Put BoundingBox into a Conduit Node.
 */
template <int NDIMS>
void put_bounding_box_to_conduit_node(const primal::BoundingBox<double, NDIMS>& bb,
                                      conduit::Node& node)
{
  node["dim"].set(bb.dimension());
  if(bb.isValid())
  {
    node["lo"].set(bb.getMin().data(), bb.dimension());
    node["hi"].set(bb.getMax().data(), bb.dimension());
  }
}

/**
 * \brief Get BoundingBox from a Conduit Node.
 */
template <int NDIMS>
void get_bounding_box_from_conduit_node(primal::BoundingBox<double, NDIMS>& bb,
                                        const conduit::Node& node)
{
  using PointType = primal::Point<double, NDIMS>;

  SLIC_ASSERT(NDIMS == node.fetch_existing("dim").as_int());

  bb.clear();

  if(node.has_child("lo"))
  {
    bb.addPoint(PointType(node.fetch_existing("lo").as_double_ptr(), NDIMS));
    bb.addPoint(PointType(node.fetch_existing("hi").as_double_ptr(), NDIMS));
  }
}

/// Helper function to extract the dimension from the coordinate values group
/// of a mesh blueprint coordset
inline int extractDimension(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node.has_child("z") ? 3 : (values_node.has_child("y") ? 2 : 1);
}

/// Helper function to extract the number of points from the coordinate values group
/// of a mesh blueprint coordset
inline int extractSize(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node["x"].dtype().number_of_elements();
}

namespace relay
{
namespace mpi
{
/**
 * \brief Sends a conduit node along with its schema using MPI_Isend
 *
 * \param [in] node node to send
 * \param [in] dest ID of MPI rank to send to
 * \param [in] tag tag for MPI message
 * \param [in] comm MPI communicator to use
 * \param [in] request object holding state for the sent data
 * \note Adapted from conduit's relay::mpi's \a send_using_schema and \a isend to use
 * non-blocking \a MPI_Isend instead of blocking \a MPI_Send
 */
inline int isend_using_schema(conduit::Node& node,
                              int dest,
                              int tag,
                              MPI_Comm comm,
                              conduit::relay::mpi::Request* request)
{
  conduit::Schema s_data_compact;

  // schema will only be valid if compact and contig
  if(node.is_compact() && node.is_contiguous())
  {
    s_data_compact = node.schema();
  }
  else
  {
    node.schema().compact_to(s_data_compact);
  }
  const std::string snd_schema_json = s_data_compact.to_json();

  conduit::Schema s_msg;
  s_msg["schema_len"].set(conduit::DataType::int64());
  s_msg["schema"].set(conduit::DataType::char8_str(snd_schema_json.size() + 1));
  s_msg["data"].set(s_data_compact);

  // create a compact schema to use
  conduit::Schema s_msg_compact;
  s_msg.compact_to(s_msg_compact);
  request->m_buffer.reset();
  request->m_buffer.set_schema(s_msg_compact);

  // set up the message's node using this schema
  request->m_buffer["schema_len"].set((int64)snd_schema_json.length());
  request->m_buffer["schema"].set(snd_schema_json);
  request->m_buffer["data"].update(node);

  // for wait_all,  this must always be NULL except for
  // the irecv cases where copy out is necessary
  // isend case must always be NULL
  request->m_rcv_ptr = nullptr;

  auto msg_data_size = request->m_buffer.total_bytes_compact();
  int mpi_error = MPI_Isend(const_cast<void*>(request->m_buffer.data_ptr()),
                            static_cast<int>(msg_data_size),
                            MPI_BYTE,
                            dest,
                            tag,
                            comm,
                            &(request->m_request));

  // Error checking -- Note: expansion of CONDUIT_CHECK_MPI_ERROR
  if(static_cast<int>(mpi_error) != MPI_SUCCESS)
  {
    char check_mpi_err_str_buff[MPI_MAX_ERROR_STRING];
    int check_mpi_err_str_len = 0;
    MPI_Error_string(mpi_error, check_mpi_err_str_buff, &check_mpi_err_str_len);

    SLIC_ERROR(
      fmt::format("MPI call failed: error code = {} error message = {}",
                  mpi_error,
                  check_mpi_err_str_buff));
  }

  return mpi_error;
}

/// A modified version of the conduit method.
// This version works correctly when src is MPI_ANY_SOURCE
// and tag is MPI_ANY_TAG.  When conduit supports this,
// this version can be removed.
inline int recv_using_schema(conduit::Node& node, int src, int tag, MPI_Comm comm)
{
  MPI_Status status;

  int mpi_error = MPI_Probe(src, tag, comm, &status);

  // CONDUIT_CHECK_MPI_ERROR(mpi_error);
  // Expand the conduit macro:
  if(static_cast<int>(mpi_error) != MPI_SUCCESS)
  {
    char check_mpi_err_str_buff[MPI_MAX_ERROR_STRING];
    int check_mpi_err_str_len = 0;
    MPI_Error_string(mpi_error, check_mpi_err_str_buff, &check_mpi_err_str_len);

    SLIC_ERROR(
      fmt::format("MPI call failed: error code = {} error message = {}",
                  mpi_error,
                  check_mpi_err_str_buff));
  }

  int buffer_size = 0;
  MPI_Get_count(&status, MPI_BYTE, &buffer_size);

  conduit::Node n_buffer(conduit::DataType::uint8(buffer_size));

  mpi_error = MPI_Recv(n_buffer.data_ptr(),
                       buffer_size,
                       MPI_BYTE,
                       status.MPI_SOURCE,
                       status.MPI_TAG,
                       comm,
                       &status);

  uint8* n_buff_ptr = (uint8*)n_buffer.data_ptr();

  conduit::Node n_msg;
  // length of the schema is sent as a 64-bit signed int
  // NOTE: we aren't using this value  ...
  n_msg["schema_len"].set_external((int64*)n_buff_ptr);
  n_buff_ptr += 8;
  // wrap the schema string
  n_msg["schema"].set_external_char8_str((char*)(n_buff_ptr));
  // create the schema
  conduit::Schema rcv_schema;
  conduit::Generator gen(n_msg["schema"].as_char8_str());
  gen.walk(rcv_schema);

  // advance by the schema length
  n_buff_ptr += n_msg["schema"].total_bytes_compact();

  // apply the schema to the data
  n_msg["data"].set_external(rcv_schema, n_buff_ptr);

  // copy out to our result node
  node.update(n_msg["data"]);

  return mpi_error;
}

}  // namespace mpi
}  // namespace relay

/**
 * \brief Implements the DistributedClosestPoint query for a specified dimension
 * using a provided execution policy (e.g. sequential, openmp, cuda, hip)
 *
 * \tparam NDIMS The dimension of the object mesh and query points
 */
template <int NDIMS>
class DistributedClosestPointImpl
{
public:
  static constexpr int DIM = NDIMS;
  using RuntimePolicy = DistributedClosestPointRuntimePolicy;
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using PointArray = axom::Array<PointType>;
  using BoxArray = axom::Array<BoxType>;

  using SeqBVHTree = spin::BVH<DIM, axom::SEQ_EXEC>;
#ifdef _AXOM_DCP_USE_OPENMP
  using OmpBVHTree = spin::BVH<DIM, axom::OMP_EXEC>;
#endif
#ifdef _AXOM_DCP_USE_CUDA
  using CudaBVHTree = spin::BVH<DIM, axom::CUDA_EXEC<256>>;
#endif
#ifdef _AXOM_DCP_USE_HIP
  using HipBVHTree = spin::BVH<DIM, axom::HIP_EXEC<256>>;
#endif

private:
  struct MinCandidate
  {
    /// Squared distance to query point
    double minSqDist {numerics::floating_point_limits<double>::max()};
    /// Index within mesh of closest element
    int minElem {-1};
    /// MPI rank of closest element
    int minRank {-1};
  };

public:
  DistributedClosestPointImpl(RuntimePolicy runtimePolicy,
                              int allocatorID,
                              bool isVerbose)
    : m_runtimePolicy(runtimePolicy)
    , m_isVerbose(isVerbose)
    , m_sqDistanceThreshold(std::numeric_limits<double>::max())
    , m_allocatorID(allocatorID)
    , m_mpiComm(MPI_COMM_NULL)
    , m_rank(-1)
    , m_nranks(-1)
    , m_objectPts(0, 0, allocatorID)
  {
    SLIC_ASSERT(allocatorID != axom::INVALID_ALLOCATOR_ID);

    setMpiCommunicator(MPI_COMM_WORLD);
  }

  /**
   * \brief Set the MPI communicator.
   */
  void setMpiCommunicator(MPI_Comm mpiComm)
  {
    m_mpiComm = mpiComm;
    MPI_Comm_rank(m_mpiComm, &m_rank);
    MPI_Comm_size(m_mpiComm, &m_nranks);
  }

  /**
   * \brief Sets the threshold for the query
   *
   * \param [in] threshold Ignore distances greater than this value.
   */
  void setSquaredDistanceThreshold(double sqThreshold)
  {
    SLIC_ERROR_IF(sqThreshold < 0.0,
                  "Squared distance-threshold must be non-negative.");
    m_sqDistanceThreshold = sqThreshold;
  }

  /*!  @brief Sets the allocator ID to the default associated with the
    execution policy
  */
  void setAllocatorID(int allocatorID)
  {
    SLIC_ASSERT(allocatorID != axom::INVALID_ALLOCATOR_ID);
    // TODO: If appropriate, how to check for compatibility with runtime policy?
    m_allocatorID = allocatorID;
  }

public:
  /**
   * Import object mesh points from the object blueprint mesh into internal memory.
   *
   * \param [in] mdMeshNode The blueprint mesh containing the object points.
   * \param [in] valuesPath The path to the mesh points.
   * \note This function currently supports mesh blueprints with the "point" topology
   *
   * TODO: See if some of the copies in this method can be optimized out.
   */
  void importObjectPoints(const conduit::Node& mdMeshNode,
                          const std::string& valuesPath)
  {
    SLIC_ASSERT(sizeof(double) * DIM == sizeof(PointType));

    // Count points in the mesh.
    int ptCount = 0;
    for(const conduit::Node& domain : mdMeshNode.children())
    {
      auto& values = domain.fetch_existing(valuesPath);
      const int N = internal::extractSize(values);
      ptCount += N;
    }

    // Copy points to internal memory
    PointArray pts(ptCount, ptCount);
    std::size_t copiedCount = 0;
    conduit::Node tmpValues;
    for(const conduit::Node& domain : mdMeshNode.children())
    {
      auto& values = domain.fetch_existing(valuesPath);

      bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(values);
      if(!isInterleaved)
      {
        conduit::blueprint::mcarray::to_interleaved(values, tmpValues);
      }
      const conduit::Node& copySrc = isInterleaved ? values : tmpValues;

      const int N = internal::extractSize(copySrc);
      const std::size_t nBytes = sizeof(double) * DIM * N;

      axom::copy(pts.data() + copiedCount,
                 copySrc.fetch_existing("x").data_ptr(),
                 nBytes);
      tmpValues.reset();
      copiedCount += N;
    }
    m_objectPts = PointArray(pts, m_allocatorID);  // copy point array to ExecSpace
  }

  /// Predicate to check if the BVH tree has been initialized
  bool isBVHTreeInitialized() const
  {
    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      return m_bvh_seq.get() != nullptr;

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      return m_bvh_omp.get() != nullptr;
#else
      break;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      return m_bvh_cuda.get() != nullptr;
#else
      break;
#endif
    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      return m_bvh_hip.get() != nullptr;
#else
      break;
#endif
    }

    return false;
  }

  /// Generates the BVH tree for the classes execution space
  bool generateBVHTree()
  {
    // Delegates to generateBVHTreeImpl<> which uses
    // the execution space templated bvh tree

    SLIC_ASSERT_MSG(!isBVHTreeInitialized(), "BVH tree already initialized");

    // In case user changed the allocator after setObjectMesh,
    // move the object point data to avoid repetitive page faults.
    if(m_objectPts.getAllocatorID() != m_allocatorID)
    {
      PointArray tmpPoints(m_objectPts, m_allocatorID);
      m_objectPts.swap(tmpPoints);
    }

    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      m_bvh_seq = std::make_unique<SeqBVHTree>();
      return generateBVHTreeImpl<SeqBVHTree>(m_bvh_seq.get());

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      m_bvh_omp = std::make_unique<OmpBVHTree>();
      return generateBVHTreeImpl<OmpBVHTree>(m_bvh_omp.get());
#else
      break;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      m_bvh_cuda = std::make_unique<CudaBVHTree>();
      return generateBVHTreeImpl<CudaBVHTree>(m_bvh_cuda.get());
#else
      break;
#endif

    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      m_bvh_hip = std::make_unique<HipBVHTree>();
      return generateBVHTreeImpl<HipBVHTree>(m_bvh_hip.get());
#else
      break;
#endif
    }

    // Fail safe -- we should never reach this line!
    SLIC_ERROR("Failed to initialize the BVH tree");

    return false;
  }

  /// Get local copy of all ranks BVH root bounding boxes.
  void gatherBVHRoots()
  {
    SLIC_ASSERT_MSG(
      isBVHTreeInitialized(),
      "BVH tree must be initialized before calling 'gatherBVHRoots");

    BoxType local_bb;
    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      local_bb = m_bvh_seq->getBounds();
      break;

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      local_bb = m_bvh_omp->getBounds();
      break;
#else
      break;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      local_bb = m_bvh_cuda->getBounds();
      break;
#else
      break;
#endif

    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      local_bb = m_bvh_hip->getBounds();
      break;
#else
      break;
#endif
    }

    gatherBoundingBoxes(local_bb, m_objectPartitionBbs);
  }

  /// Allgather one bounding box from each rank.
  void gatherBoundingBoxes(const BoxType& aabb, BoxArray& all_aabbs) const
  {
    axom::Array<double> sendbuf(2 * DIM);
    aabb.getMin().to_array(&sendbuf[0]);
    aabb.getMax().to_array(&sendbuf[DIM]);
    axom::Array<double> recvbuf(m_nranks * sendbuf.size());
    // Note: Using axom::Array<double,2> may reduce clutter a tad.
    int errf = MPI_Allgather(sendbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             recvbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             m_mpiComm);
    SLIC_ASSERT(errf == MPI_SUCCESS);

    all_aabbs.clear();
    all_aabbs.reserve(m_nranks);
    for(int i = 0; i < m_nranks; ++i)
    {
      PointType lower(&recvbuf[i * 2 * DIM]);
      PointType upper(&recvbuf[i * 2 * DIM + DIM]);
      all_aabbs.emplace_back(BoxType(lower, upper, false));
    }
  }

  /// Compute bounding box for local part of a mesh.
  BoxType computeMeshBoundingBox(conduit::Node& xferNode) const
  {
    BoxType rval;

    conduit::Node& xferDoms = xferNode.fetch_existing("xferDoms");
    for(conduit::Node& xferDom : xferDoms.children())
    {
      const int qPtCount = xferDom.fetch_existing("qPtCount").value();

      /// Extract fields from the input node as ArrayViews
      auto queryPts =
        ArrayView_from_Node<PointType>(xferDom.fetch_existing("coords"),
                                       qPtCount);
      for(const auto& p : queryPts)
      {
        rval.addPoint(p);
      }
    }

    return rval;
  }

  /*!
   * Copy parts of query mesh partition to a conduit::Node for
   * computation and communication.
   * queryNode must be a blueprint multidomain mesh.
   */
  void node_copy_query_to_xfer(conduit::Node& queryNode,
                               conduit::Node& xferNode,
                               const std::string& coordset) const
  {
    xferNode["homeRank"] = m_rank;
    xferNode["is_first"] = 1;

    conduit::Node& xferDoms = xferNode["xferDoms"];
    for(auto& queryDom : queryNode.children())
    {
      const std::string& domName = queryDom.name();
      conduit::Node& xferDom = xferDoms[domName];
      conduit::Node& fields = queryDom.fetch_existing("fields");
      conduit::Node& queryCoords =
        queryDom.fetch_existing(fmt::format("coordsets/{}", coordset));
      conduit::Node& queryCoordsValues = queryCoords.fetch_existing("values");

      // clang-format off
      copy_components_to_interleaved(queryCoordsValues, xferDom["coords"]);

      xferDom["cp_index"].set_external(fields.fetch_existing("cp_index/values"));
      xferDom["cp_rank"].set_external(fields.fetch_existing("cp_rank/values"));
      copy_components_to_interleaved(fields.fetch_existing("cp_coords/values"),
                                     xferDom["cp_coords"]);

      if(fields.has_path("cp_distance"))
      {
        xferDom["debug/cp_distance"].set_external(fields.fetch_existing("cp_distance/values"));
      }
      // clang-format on

      const int dim = internal::extractDimension(queryCoordsValues);
      const int qPtCount = internal::extractSize(queryCoordsValues);
      xferDom["qPtCount"] = qPtCount;
      xferDom["dim"] = dim;
    }
  }

  /// Copy xferNode back to query mesh partition.
  void node_copy_xfer_to_query(const conduit::Node& xferNode,
                               conduit::Node& queryNode) const
  {
    const conduit::Node& xferDoms = xferNode.fetch_existing("xferDoms");
    conduit::index_t childCount = queryNode.number_of_children();
    SLIC_ASSERT(xferDoms.number_of_children() == childCount);
    for(conduit::index_t ci = 0; ci < childCount; ++ci)
    {
      const conduit::Node& xferDom = xferDoms.child(ci);
      conduit::Node& queryDom = queryNode.child(ci);
      conduit::Node& fields = queryDom.fetch_existing("fields");

      {
        auto& src = xferDom.fetch_existing("cp_rank");
        auto& dst = fields.fetch_existing("cp_rank/values");
        if(dst.data_ptr() != src.data_ptr())
        {
          dst.update_compatible(src);
        }
      }

      {
        auto& src = xferDom.fetch_existing("cp_index");
        auto& dst = fields.fetch_existing("cp_index/values");
        if(dst.data_ptr() != src.data_ptr())
        {
          dst.update_compatible(src);
        }
      }

      copy_interleaved_to_components(xferDom.fetch_existing("cp_coords"),
                                     fields.fetch_existing("cp_coords/values"));

      if(xferDom.has_path("debug/cp_distance"))
      {
        auto& src = xferDom.fetch_existing("debug/cp_distance");
        auto& dst = fields.fetch_existing("cp_distance/values");
        if(dst.data_ptr() != src.data_ptr())
        {
          dst.update_compatible(src);
        }
      }
    }
  }

  /*
    Special copy from coordinates (in a format that's not
    necessarily interleaved) to a 1D array of interleaved values).
    If coordinates are already interleaved, copy pointer.
  */
  void copy_components_to_interleaved(conduit::Node& components,
                                      conduit::Node& interleaved) const
  {
    const int dim = internal::extractDimension(components);
    const int qPtCount = internal::extractSize(components);
    bool interleavedSrc = conduit::blueprint::mcarray::is_interleaved(components);
    if(interleavedSrc)
    {
      interleaved.set_external(internal::getPointer<double>(components.child(0)),
                               dim * qPtCount);
    }
    else
    {
      // Copy from component-wise src to 1D-interleaved dst.
      interleaved.reset();
      interleaved.set_dtype(conduit::DataType::float64(dim * qPtCount));
      for(int d = 0; d < dim; ++d)
      {
        auto src = components.child(d).as_float64_array();
        double* dst = interleaved.as_float64_ptr() + d;
        for(int i = 0; i < qPtCount; ++i)
        {
          dst[i * dim] = src[i];
        }
      }
    }
  }

  /*
    Special copy from 1D interleaved coordinate values back to
    component-wise storage.
    This is a nop if they point to the same data.
  */
  void copy_interleaved_to_components(const conduit::Node& interleaved,
                                      conduit::Node& components) const
  {
    if(interleaved.data_ptr() != components.child(0).data_ptr())
    {
      const int dim = internal::extractDimension(components);
      const int qPtCount = internal::extractSize(components);
      // Copy from 1D-interleaved src to component-wise dst.
      for(int d = 0; d < dim; ++d)
      {
        const double* src = interleaved.as_float64_ptr() + d;
        auto dst = components.child(d).as_float64_array();
        for(int i = 0; i < qPtCount; ++i)
        {
          dst[i] = src[i * dim];
        }
      }
    }
  }

  /**
   * \brief Computes the closest point within the objects for each query point
   * in the provided particle mesh, provided in the mesh blueprint rooted at \a queryMesh
   *
   * \param queryMesh The root node of a mesh blueprint for the query points
   * Can be empty if there are no query points for the calling rank
   * \param coordset The coordinate set for the query points
   *
   * When the query mesh contains query points, it uses the \a coordset coordinate set
   * of the provided blueprint mesh and  contains the following fields:
   *   - cp_rank: Will hold the rank of the object point containing the closest point
   *   - cp_index: Will hold the index of the object point containing the closest point.
   *     For multiple object mesh domains on a rank, cp_index increases consecutively
   *     from one domain to the next.
   *   - cp_coords: Will hold the coordinates of the closest points
   *     interleaved in a 1D array.
   *
   * \note The current implementation assumes that the coordinates and
   * cp_coords are interleaved or contiguous.
   *
   * \note We're temporarily also using a cp_distance field while debugging this class.
   * The code will use this field if it is present in \a queryMesh.
   *
   * We use non-blocking sends for performance and deadlock avoidance.
   * The worst case could incur nranks^2 sends.  To avoid excessive
   * buffer usage, we occasionally check the sends for completion,
   * using check_send_requests().
   */
  void computeClosestPoints(conduit::Node& queryMesh_,
                            const std::string& coordset) const
  {
    SLIC_ASSERT_MSG(
      isBVHTreeInitialized(),
      "BVH tree must be initialized before calling 'computeClosestPoints");

    // If query mesh isn't multidomain, create a temporary multidomain representation.
    const bool qmIsMultidomain =
      conduit::blueprint::mesh::is_multi_domain(queryMesh_);
    std::shared_ptr<conduit::Node> tmpNode;
    if(!qmIsMultidomain)
    {
      tmpNode = std::make_shared<conduit::Node>();
      conduit::blueprint::mesh::to_multi_domain(queryMesh_, *tmpNode);
    }
    conduit::Node& queryMesh = qmIsMultidomain ? queryMesh_ : *tmpNode;

    std::map<int, std::shared_ptr<conduit::Node>> xferNodes;

    // create conduit Node containing data that has to xfer between ranks.
    // The node will be mostly empty if there are no domains on this rank
    {
      xferNodes[m_rank] = std::make_shared<conduit::Node>();
      conduit::Node& xferNode = *xferNodes[m_rank];
      node_copy_query_to_xfer(queryMesh, xferNode, coordset);
      xferNode["homeRank"] = m_rank;
    }

    BoxType myQueryBb = computeMeshBoundingBox(*xferNodes[m_rank]);
    put_bounding_box_to_conduit_node(myQueryBb, xferNodes[m_rank]->fetch("aabb"));
    BoxArray allQueryBbs;
    gatherBoundingBoxes(myQueryBb, allQueryBbs);

    {
      conduit::Node& xferNode = *xferNodes[m_rank];
      computeLocalClosestPointsByPolicy(xferNode);
    }

    const auto& myObjectBb = m_objectPartitionBbs[m_rank];
    int remainingRecvs = 0;
    for(int r = 0; r < m_nranks; ++r)
    {
      if(r != m_rank)
      {
        const auto& otherQueryBb = allQueryBbs[r];
        double sqDistance =
          axom::primal::squared_distance(otherQueryBb, myObjectBb);
        if(sqDistance <= m_sqDistanceThreshold)
        {
          ++remainingRecvs;
        }
      }
    }

    // arbitrary tags for send/recv xferNode.
    const int tag = 987342;

    std::list<conduit::relay::mpi::Request> isendRequests;

    {
      /*
        Send local query mesh to next rank with close-enough object
        partition, if any.  Increase remainingRecvs, because this data
        will come back.
      */
      int firstRecipForMyQuery = next_recipient(*xferNodes[m_rank]);
      if(m_nranks == 1)
      {
        SLIC_ASSERT(firstRecipForMyQuery == -1);
      }

      if(firstRecipForMyQuery == -1)
      {
        // No need to send anywhere.  Put computed data back into queryMesh.
        node_copy_xfer_to_query(*xferNodes[m_rank], queryMesh);
        xferNodes.erase(m_rank);
      }
      else
      {
        isendRequests.emplace_back(conduit::relay::mpi::Request());
        auto& req = isendRequests.back();
        relay::mpi::isend_using_schema(*xferNodes[m_rank],
                                       firstRecipForMyQuery,
                                       tag,
                                       m_mpiComm,
                                       &req);
        ++remainingRecvs;
      }
    }

    while(remainingRecvs > 0)
    {
      SLIC_INFO_IF(
        m_isVerbose,
        fmt::format("=======  {} receives remaining =======", remainingRecvs));

      // Receive the next xferNode
      std::shared_ptr<conduit::Node> recvXferNodePtr =
        std::make_shared<conduit::Node>();
      relay::mpi::recv_using_schema(*recvXferNodePtr,
                                    MPI_ANY_SOURCE,
                                    tag,
                                    m_mpiComm);

      const int homeRank = recvXferNodePtr->fetch_existing("homeRank").as_int();
      --remainingRecvs;
      xferNodes[homeRank] = recvXferNodePtr;
      conduit::Node& xferNode = *xferNodes[homeRank];

      if(homeRank == m_rank)
      {
        node_copy_xfer_to_query(xferNode, queryMesh);
      }
      else
      {
        computeLocalClosestPointsByPolicy(xferNode);

        isendRequests.emplace_back(conduit::relay::mpi::Request());
        auto& isendRequest = isendRequests.back();
        int nextRecipient = next_recipient(xferNode);
        SLIC_ASSERT(nextRecipient != -1);
        relay::mpi::isend_using_schema(xferNode,
                                       nextRecipient,
                                       tag,
                                       m_mpiComm,
                                       &isendRequest);

        // Check non-blocking sends to free memory.
        check_send_requests(isendRequests, false);
      }

    }  // remainingRecvs loop

    // Complete remaining non-blocking sends.
    while(!isendRequests.empty())
    {
      check_send_requests(isendRequests, true);
    }

    MPI_Barrier(m_mpiComm);
    slic::flushStreams();
  }

private:
  /// Distance search using local object partition and xferNode.
  void computeLocalClosestPointsByPolicy(conduit::Node& xferNode) const
  {
    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      computeLocalClosestPoints<SeqBVHTree>(m_bvh_seq.get(), xferNode);
      break;

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      computeLocalClosestPoints<OmpBVHTree>(m_bvh_omp.get(), xferNode);
#endif
      break;

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      computeLocalClosestPoints<CudaBVHTree>(m_bvh_cuda.get(), xferNode);
#endif
      break;

    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      computeLocalClosestPoints<HipBVHTree>(m_bvh_hip.get(), xferNode);
#endif
      break;
    }
  }

  /**
    Determine the next rank (in ring order) with an object partition
    close to the query points in xferNode.  The intent is to send
    xferNode there next.
  */
  int next_recipient(const conduit::Node& xferNode) const
  {
    int homeRank = xferNode.fetch_existing("homeRank").value();
    BoxType bb;
    get_bounding_box_from_conduit_node(bb, xferNode.fetch_existing("aabb"));
    for(int i = 1; i < m_nranks; ++i)
    {
      int maybeNextRecip = (m_rank + i) % m_nranks;
      if(maybeNextRecip == homeRank)
      {
        return maybeNextRecip;
      }
      double sqDistance =
        primal::squared_distance(bb, m_objectPartitionBbs[maybeNextRecip]);
      if(sqDistance <= m_sqDistanceThreshold)
      {
        return maybeNextRecip;
      }
    }
    return -1;
  }

  /// Wait for some non-blocking sends (if any) to finish.
  void check_send_requests(std::list<conduit::relay::mpi::Request>& isendRequests,
                           bool atLeastOne) const
  {
    std::vector<MPI_Request> reqs;
    for(auto& isr : isendRequests) reqs.push_back(isr.m_request);

    int inCount = static_cast<int>(reqs.size());
    int outCount = 0;
    std::vector<int> indices(reqs.size(), -1);
    if(atLeastOne)
      MPI_Waitsome(inCount,
                   reqs.data(),
                   &outCount,
                   indices.data(),
                   MPI_STATUSES_IGNORE);
    else
      MPI_Testsome(inCount,
                   reqs.data(),
                   &outCount,
                   indices.data(),
                   MPI_STATUSES_IGNORE);
    indices.resize(outCount);

    auto reqIter = isendRequests.begin();
    int prevIdx = 0;
    for(const int idx : indices)
    {
      for(; prevIdx < idx; ++prevIdx) ++reqIter;
      reqIter = isendRequests.erase(reqIter);
      ++prevIdx;
    }
  }

  /**
   * \brief Extracts a field \a fieldName from the mesh blueprint
   *
   * \tparam T The type for the underlying array
   * \param mesh_node The conduit node at the root of the mesh blueprint
   * \param field_name The name of the field
   * \param field_template Template string for the path to the field
   * \param num_objectPts The size of the field
   * \return An arrayview over the field data
   */
  template <typename T>
  axom::ArrayView<T> extractField(conduit::Node& mesh_node,
                                  std::string&& field_name,
                                  std::string&& path_template,
                                  int num_objectPts) const
  {
    const std::string path = axom::fmt::format(path_template, field_name);
    SLIC_ASSERT_MSG(
      mesh_node.has_path(path),
      fmt::format(
        "Input to `computeClosestPoint()` must have a field named `{}`",
        field_name));

    return internal::ArrayView_from_Node<T>(mesh_node[path], num_objectPts);
  }

  // Note: following should be private, but nvcc complains about lambdas in private scope
public:
  /// Templated implementation of generateBVHTree function
  template <typename BVHTreeType>
  bool generateBVHTreeImpl(BVHTreeType* bvh)
  {
    using ExecSpace = typename BVHTreeType::ExecSpaceType;

    SLIC_ASSERT(bvh != nullptr);

    const int npts = m_objectPts.size();
    axom::Array<BoxType> boxesArray(npts, npts, m_allocatorID);
    auto boxesView = boxesArray.view();
    auto pointsView = m_objectPts.view();

    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { boxesView[i] = BoxType {pointsView[i]}; });

    // Build bounding volume hierarchy
    bvh->setAllocatorID(m_allocatorID);
    int result = bvh->initialize(boxesView, npts);

    gatherBVHRoots();

    return (result == spin::BVH_BUILD_OK);
  }

  /**
   * This method assumes xferNode is a blueprint single-domain mesh.
   */
  template <typename BVHTreeType>
  void computeLocalClosestPoints(const BVHTreeType* bvh,
                                 conduit::Node& xferNode) const
  {
    using ExecSpace = typename BVHTreeType::ExecSpaceType;
    using axom::primal::squared_distance;
    using int32 = axom::int32;

    // Note: There is some additional computation the first time this function
    // is called for a query node, even if the local object mesh is empty
    const bool hasObjectPoints = m_objectPts.size() > 0;
    const bool is_first = xferNode.has_path("is_first");
    if(!hasObjectPoints && !is_first)
    {
      return;
    }

    conduit::Node& xferDoms = xferNode["xferDoms"];
    for(conduit::Node& xferDom : xferDoms.children())
    {
      // --- Set up arrays and views in the execution space
      // Arrays are initialized in that execution space the first time they are processed
      // and are copied in during subsequent processing

      // Check dimension and extract the number of points
      SLIC_ASSERT(xferDom.fetch_existing("dim").as_int() == NDIMS);
      const int qPtCount = xferDom.fetch_existing("qPtCount").value();

      /// Extract fields from the input node as ArrayViews
      auto queryPts =
        ArrayView_from_Node<PointType>(xferDom.fetch_existing("coords"),
                                       qPtCount);
      auto cpIndexes =
        ArrayView_from_Node<axom::IndexType>(xferDom.fetch_existing("cp_index"),
                                             qPtCount);
      auto cpRanks =
        ArrayView_from_Node<axom::IndexType>(xferDom.fetch_existing("cp_rank"),
                                             qPtCount);
      auto cpCoords =
        ArrayView_from_Node<PointType>(xferDom.fetch_existing("cp_coords"),
                                       qPtCount);

      /// Create ArrayViews in ExecSpace that are compatible with fields
      // This deep-copies host memory in xferDom to device memory.
      // TODO: Avoid copying arrays (here and at the end) if both are on the host
      auto cp_idx = is_first
        ? axom::Array<axom::IndexType>(qPtCount, qPtCount, m_allocatorID)
        : axom::Array<axom::IndexType>(cpIndexes, m_allocatorID);
      auto cp_rank = is_first
        ? axom::Array<axom::IndexType>(qPtCount, qPtCount, m_allocatorID)
        : axom::Array<axom::IndexType>(cpRanks, m_allocatorID);

      /// PROBLEM: The striding does not appear to be retained by conduit relay
      ///          We might need to transform it? or to use a single array w/ pointers into it?
      auto cp_pos = is_first
        ? axom::Array<PointType>(qPtCount, qPtCount, m_allocatorID)
        : axom::Array<PointType>(cpCoords, m_allocatorID);

      // DEBUG
      const bool has_cp_distance = xferDom.has_path("debug/cp_distance");
      auto minDist = has_cp_distance
        ? ArrayView_from_Node<double>(
            xferDom.fetch_existing("debug/cp_distance"),
            qPtCount)
        : ArrayView<double>();

      auto cp_dist = has_cp_distance
        ? (is_first ? axom::Array<double>(qPtCount, qPtCount, m_allocatorID)
                    : axom::Array<double>(minDist, m_allocatorID))
        : axom::Array<double>(0, 0, m_allocatorID);
      auto query_min_dist = cp_dist.view();
      // BTNG Q: Why do we need query_min_dist?  Why not use cp_dist?
      // END DEBUG

      if(is_first)
      {
        cp_rank.fill(-1);
        cp_idx.fill(-1);
        const PointType nowhere(std::numeric_limits<double>::signaling_NaN());
        cp_pos.fill(nowhere);
        if(has_cp_distance)
        {
          cp_dist.fill(std::numeric_limits<double>::signaling_NaN());
        }
      }
      auto query_inds = cp_idx.view();
      auto query_ranks = cp_rank.view();
      auto query_pos = cp_pos.view();

      /// Create an ArrayView in ExecSpace that is compatible with queryPts
      PointArray execPoints(queryPts, m_allocatorID);
      auto query_pts = execPoints.view();

      if(hasObjectPoints)
      {
        // Get a device-useable iterator
        auto it = bvh->getTraverser();
        const int rank = m_rank;

        double* sqDistThresh = axom::allocate<double>(
          1,
          axom::execution_space<ExecSpace>::allocatorID());
        *sqDistThresh = m_sqDistanceThreshold;

        auto pointsView = m_objectPts.view();

        AXOM_PERF_MARK_SECTION(
          "ComputeClosestPoints",
          axom::for_all<ExecSpace>(
            qPtCount,
            AXOM_LAMBDA(int32 idx) mutable {
              PointType qpt = query_pts[idx];

              MinCandidate curr_min {};
              // Preset cur_min to the closest point found so far.
              if(query_ranks[idx] >= 0)
              {
                curr_min.minSqDist = squared_distance(qpt, query_pos[idx]);
                curr_min.minElem = query_inds[idx];
                curr_min.minRank = query_ranks[idx];
              }

              auto checkMinDist = [&](int32 current_node,
                                      const int32* leaf_nodes) {
                const int candidate_idx = leaf_nodes[current_node];
                const PointType candidate_pt = pointsView[candidate_idx];
                const double sq_dist = squared_distance(qpt, candidate_pt);

                if(sq_dist < curr_min.minSqDist)
                {
                  curr_min.minSqDist = sq_dist;
                  curr_min.minElem = candidate_idx;
                  curr_min.minRank = rank;
                }
              };

              auto traversePredicate = [&](const PointType& p,
                                           const BoxType& bb) -> bool {
                auto sqDist = squared_distance(p, bb);
                return sqDist <= curr_min.minSqDist && sqDist <= sqDistThresh[0];
              };

              // Traverse the tree, searching for the point with minimum distance.
              it.traverse_tree(qpt, checkMinDist, traversePredicate);

              // If modified, update the fields that changed
              if(curr_min.minRank == rank)
              {
                query_inds[idx] = curr_min.minElem;
                query_ranks[idx] = curr_min.minRank;
                query_pos[idx] = pointsView[curr_min.minElem];

                //DEBUG
                if(has_cp_distance)
                {
                  query_min_dist[idx] = sqrt(curr_min.minSqDist);
                }
              }
            }););

        axom::deallocate(sqDistThresh);
      }

      axom::copy(cpIndexes.data(),
                 query_inds.data(),
                 cpIndexes.size() * sizeof(axom::IndexType));
      axom::copy(cpRanks.data(),
                 query_ranks.data(),
                 cpRanks.size() * sizeof(axom::IndexType));
      if(1)
        axom::copy(cpCoords.data(),
                   query_pos.data(),
                   cpCoords.size() * sizeof(PointType));

      // DEBUG
      if(has_cp_distance)
      {
        axom::copy(minDist.data(),
                   query_min_dist.data(),
                   minDist.size() * sizeof(double));
      }
    }

    // Data has now been initialized
    if(is_first)
    {
      xferNode.remove_child("is_first");
    }
  }

private:
  RuntimePolicy m_runtimePolicy;
  bool m_isVerbose {false};
  double m_sqDistanceThreshold;
  int m_allocatorID;
  MPI_Comm m_mpiComm;
  bool m_mpiCommIsPrivate;
  int m_rank;
  int m_nranks;

  /*!
    @brief Object point array.

    Points from all object mesh domains are flattened here.
    TODO: Will need a way to map the m_objectPts index back
    to the domain that has the object point.
  */
  PointArray m_objectPts;

  /*!  @brief Object partition bounding boxes, one per rank.
    All are in physical space, not index space.
  */
  BoxArray m_objectPartitionBbs;

  std::unique_ptr<SeqBVHTree> m_bvh_seq;

#ifdef _AXOM_DCP_USE_OPENMP
  std::unique_ptr<OmpBVHTree> m_bvh_omp;
#endif

#ifdef _AXOM_DCP_USE_CUDA
  std::unique_ptr<CudaBVHTree> m_bvh_cuda;
#endif

#ifdef _AXOM_DCP_USE_HIP
  std::unique_ptr<HipBVHTree> m_bvh_hip;
#endif
};  // DistributedClosestPointImpl

}  // namespace internal

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
 * \note The class currently supports object meshes that are comprised of a collection of points.
 * In the future, we'd like to consider more general object meshes, e.g. triangle meshes.
 *
 * To use this class, first set some parameters, such as the runtime execution policy,
 * then pass in the object mesh and build a spatial index over this mesh.
 * Finally, compute the closest points in the object mesh to each point in a query mesh
 * using the \a computeClosestPoint() function.
 *
 * \note The implementation currently assumes that the coordinates for the positions and vector field
 * data are interleaved (i.e. xyzxyzxyz....). We will relax this assumption in the future to support both
 * interleaved and strided data.
 *
 * \note To prevent mixing unrelated MPI communications, you can set a
 * custom MPI Communicator using setMpiCommunicator().
 */
class DistributedClosestPoint
{
public:
  using RuntimePolicy = DistributedClosestPointRuntimePolicy;

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

  /// Set the runtime execution policy for the query
  void setRuntimePolicy(RuntimePolicy policy)
  {
    SLIC_ASSERT_MSG(
      isValidRuntimePolicy(policy),
      fmt::format("Policy '{}' is not a valid runtime policy", policy));
    m_runtimePolicy = policy;
  }

  /// Predicate to determine if a given \a RuntimePolicy is valid for this configuration
  bool isValidRuntimePolicy(RuntimePolicy policy) const
  {
    switch(policy)
    {
    case RuntimePolicy::seq:
      return true;

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      return true;
#else
      return false;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      return true;
#else
      return false;
#endif
    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      return true;
#else
      return false;
#endif
    }

    return false;
  }

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
    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      defaultAllocatorID = axom::execution_space<axom::OMP_EXEC>::allocatorID();
#endif
      break;

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      defaultAllocatorID =
        axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
#endif
      break;

    case RuntimePolicy::hip:
#ifdef _AXOM_DCP_USE_HIP
      defaultAllocatorID =
        axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
#endif
      break;
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

  /// Sets the logging verbosity of the query. By default the query is not verbose
  void setVerbosity(bool isVerbose) { m_isVerbose = isVerbose; }

  /**
   * \brief Sets the object mesh for the query
   *
   * \param [in] meshNode Conduit node for the object mesh
   * \param [in] coordset The name of the coordset for the object mesh's coordinates
   *
   * \pre \a meshNode must follow the mesh blueprint convention.
   * \pre Dimension of the mesh must be 2D or 3D
   */
  void setObjectMesh(const conduit::Node& meshNode, const std::string& coordset)
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

    auto domainCount = conduit::blueprint::mesh::number_of_domains(mdMeshNode);
    const std::string valuesPath = fmt::format("coordsets/{}/values", coordset);

    // Extract the dimension from the coordinate values group
    // use allreduce since some ranks might be empty
    {
      int localDim = -1;
      if(domainCount > 0)
      {
        const conduit::Node& domain0(mdMeshNode[0]);
        SLIC_ASSERT(domain0.has_path(valuesPath));
        auto& values = domain0.fetch_existing(valuesPath);
        localDim = internal::extractDimension(values);
      }
      int dim = -1;
      MPI_Allreduce(&localDim, &dim, 1, MPI_INT, MPI_MAX, m_mpiComm);
      setDimension(dim);
    }

    allocateQueryInstance();

    switch(m_dimension)
    {
    case 2:
      m_dcp_2->importObjectPoints(mdMeshNode, valuesPath);
      break;
    case 3:
      m_dcp_3->importObjectPoints(mdMeshNode, valuesPath);
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
   * \param [in] query_node conduit node containing the query points
   * \param [in] coordset The name of the coordinates within query_node
   *
   * \pre query_node must follow the conduit mesh blueprint convention
   */
  void computeClosestPoints(conduit::Node& query_node, const std::string& coordset)
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
      m_dcp_2->computeClosestPoints(query_node, coordset);
      break;
    case 3:
      m_dcp_3->setSquaredDistanceThreshold(m_sqDistanceThreshold);
      m_dcp_3->setMpiCommunicator(m_mpiComm);
      m_dcp_3->computeClosestPoints(query_node, coordset);
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

private:
  RuntimePolicy m_runtimePolicy {RuntimePolicy::seq};
  MPI_Comm m_mpiComm;
  bool m_mpiCommIsPrivate;
  int m_allocatorID {axom::INVALID_ALLOCATOR_ID};
  int m_dimension {-1};
  bool m_isVerbose {false};
  double m_sqDistanceThreshold {std::numeric_limits<double>::max()};

  bool m_objectMeshCreated {false};

  // One instance per dimension
  std::unique_ptr<internal::DistributedClosestPointImpl<2>> m_dcp_2;
  std::unique_ptr<internal::DistributedClosestPointImpl<3>> m_dcp_3;
};

}  // end namespace quest
}  // end namespace axom

// Cleanup local #defines
#undef _AXOM_DCP_USE_OPENMP
#undef _AXOM_DCP_USE_CUDA
#undef _AXOM_DCP_USE_HIP

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
