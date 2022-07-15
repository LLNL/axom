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
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_mpi.hpp"
#include "conduit_relay_io.hpp"

#include <memory>
#include <limits>
#include <cstdlib>
#include <cmath>

#ifndef AXOM_USE_MPI
  #error This file requires Axom to be configured with MPI
#endif
#include "mpi.h"

// Add some helper preprocessor defines for using OPENMP and CUDA policies
// within the distributed closest point query.
// These are only used when building with RAJA and Umpire
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #ifdef AXOM_USE_OPENMP
    #define _AXOM_DCP_USE_OPENMP
  #endif
  #ifdef AXOM_USE_CUDA
    #define _AXOM_DCP_USE_CUDA
  #endif
#endif

#include <axom/core/utilities/WhereMacro.hpp>
namespace axom
{
namespace quest
{
/// Enum for runtime execution policy
enum class DistributedClosestPointRuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2
};

namespace internal
{
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
axom::ArrayView<primal::Point<double, 2>> ArrayView_from_Node(conduit::Node& node,
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
axom::ArrayView<primal::Point<double, 3>> ArrayView_from_Node(conduit::Node& node,
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
void put_to_conduit_node(const primal::BoundingBox<double, NDIMS>& bb,
                         conduit::Node& node)
{
  node["dim"].set(bb.dimension());
  node["lo"].set(bb.getMin().data(), bb.dimension());
  node["hi"].set(bb.getMax().data(), bb.dimension());
}

/**
 * \brief Get BoundingBox from a Conduit Node.
 */
template <int NDIMS>
void get_from_conduit_node(primal::BoundingBox<double, NDIMS>& bb,
                           const conduit::Node& node)
{
  int dim = node.fetch_existing("dim").as_int();
  SLIC_ASSERT(dim == NDIMS);
  const double* lo = node.fetch_existing("lo").as_double_ptr();
  const double* hi = node.fetch_existing("hi").as_double_ptr();
  bb =
    primal::BoundingBox<double, NDIMS>(primal::Point<double, NDIMS>(lo, NDIMS),
                                       primal::Point<double, NDIMS>(hi, NDIMS));
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
/// Struct to hold state related to the isend_using_schema function
struct ISendRequest
{
  MPI_Request mpi_request;
  conduit::Schema node_schema;
  conduit::Schema msg_schema;
  conduit::Node msg_node;
  ISendRequest()
    : mpi_request(MPI_REQUEST_NULL)
    {}
};

/**
 * \brief Sends a conduit node along with its schema using MPI_Isend
 *
 * \param [in] node node to send
 * \param [in] dest ID of MPI rank to send to
 * \param [in] tag tag for MPI message
 * \param [in] comm MPI communicator to use
 * \param [in] request An instance of ISendRequest that holds state for the sent data
 * \note Adapted from conduit's relay::mpi's \a send_using_schema and \a isend to use
 * non-blocking \a MPI_Isend instead of blocking \a MPI_Send
 */
inline int isend_using_schema(conduit::Node& node,
                              int dest,
                              int tag,
                              MPI_Comm comm,
                              ISendRequest* request)
{
  // schema will only be valid if compact and contig
  if(node.is_compact() && node.is_contiguous())
  {
    request->node_schema = node.schema();
  }
  else
  {
    node.schema().compact_to(request->node_schema);
  }
  const std::string snd_schema_json = request->node_schema.to_json();

  // create a compact schema to use
  conduit::Schema s_msg;
  s_msg["schema_len"].set(conduit::DataType::int64());
  s_msg["schema"].set(conduit::DataType::char8_str(snd_schema_json.size() + 1));
  s_msg["data"].set(request->node_schema);
  s_msg.compact_to(request->msg_schema);
  request->msg_node.reset();
  request->msg_node.set_schema(request->msg_schema);

  // set up the message's node using this schema
  request->msg_node["schema_len"].set((int64)snd_schema_json.length());
  request->msg_node["schema"].set(snd_schema_json);
  request->msg_node["data"].update(node);

  auto msg_data_size = request->msg_node.total_bytes_compact();
  int mpi_error = MPI_Isend(const_cast<void*>(request->msg_node.data_ptr()),
                            static_cast<int>(msg_data_size),
                            MPI_BYTE,
                            dest,
                            tag,
                            comm,
                            &(request->mpi_request));

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

/**
 * \brief Orchestrates sending and receiving conduit nodes among mpi ranks
 *
 * \param [in] send_node rooted node to send
 * \param [in] recv_node rooted node to receive
 * \param [in] send_rank MPI rank to send \a send_node to
 * \param [in] recv_rank MPI rank to receive \a recv_node from
 * \param [in] tag tag for MPI messages
 * \param [in] comm MPI communicator to use
 *
 * Sends/receives a pair of conduit nodes (along with their schemas) from the current rank.
 * Uses non-blocking send (isend) and blocking receives and ensures that the request
 * objects associated with the send are finalized
 */
inline void send_and_recv_node(conduit::Node& send_node,
                               conduit::Node& recv_node,
                               int send_rank,
                               int recv_rank,
                               int tag,
                               MPI_Comm comm)
{
  ISendRequest req;

  // non-blocking send
  isend_using_schema(send_node, send_rank, tag, comm, &req);

  // blocking receive
  conduit::relay::mpi::recv_using_schema(recv_node, recv_rank, tag, comm);

  // sender blocks until receiver is done
  MPI_Wait(&(req.mpi_request), MPI_STATUS_IGNORE);
}

}  // namespace mpi
}  // namespace relay

/**
 * \brief Implements the DistributedClosestPoint query for a specified dimension
 * using a provided execution policy (e.g. sequential, openmp, cuda)
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
  DistributedClosestPointImpl(RuntimePolicy runtimePolicy, bool isVerbose)
    : m_runtimePolicy(runtimePolicy)
    , m_isVerbose(isVerbose)
    , m_sqDistanceThreshold(std::numeric_limits<double>::max())
  {
    setAllocatorID();

    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
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

public:
  /**
   * Utility function to set the array of points
   *
   * \param [in] coords The root group of a mesh blueprint's coordinate values
   * \note This function currently supports mesh blueprints with the "point" topology
   */
  void importObjectPoints(const conduit::Node& coords, int nPts)
  {
    // Extract pointers to the coordinate data
    // Note: The following assumes the coordinates are contiguous with stride NDIMS
    // TODO: Generalize to support other strides

    // TODO: Add error checking for children 'x', 'y' and 'z' and striding

    // Copy the data into the point array of primal points
    PointArray pts(nPts, nPts);
    const std::size_t nbytes = sizeof(double) * DIM * nPts;
    axom::copy(pts.data(), coords["x"].data_ptr(), nbytes);

    m_points = PointArray(pts, m_allocatorID);  // copy point array to ExecSpace
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
    }

    return false;
  }

  /// Generates the BVH tree for the classes execution space
  bool generateBVHTree()
  {
    // Delegates to generateBVHTreeImpl<> which uses
    // the execution space templated bvh tree

    SLIC_ASSERT_MSG(!isBVHTreeInitialized(), "BVH tree already initialized");

    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      m_bvh_seq = std::unique_ptr<SeqBVHTree>(new SeqBVHTree);
      return generateBVHTreeImpl<SeqBVHTree>(m_bvh_seq.get());

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      m_bvh_omp = std::unique_ptr<OmpBVHTree>(new OmpBVHTree);
      return generateBVHTreeImpl<OmpBVHTree>(m_bvh_omp.get());
#else
      break;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      m_bvh_cuda = std::unique_ptr<CudaBVHTree>(new CudaBVHTree);
      return generateBVHTreeImpl<CudaBVHTree>(m_bvh_cuda.get());
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
#else
      break;
#endif

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      local_bb = m_bvh_cuda->getBounds();
#else
      break;
#endif
    }

std::cout << fmt::format("at {}, {} has objectBB {}", __WHERE, m_rank, local_bb) << std::endl;
    gatherBoundingBoxes(local_bb, m_objectPartitionBbs);
  }

  /// Allgather one bounding box from each rank.
  void gatherBoundingBoxes(const BoxType& aabb, BoxArray& all_aabbs) const
  {
    // Using MPI calls.  Should we change to Conduit relay MPI?
    Array<double> sendbuf(2 * DIM);
    aabb.getMin().to_array(&sendbuf[0]);
    aabb.getMax().to_array(&sendbuf[DIM]);
    Array<double> recvbuf(m_nranks * sendbuf.size());
    int errf = MPI_Allgather(sendbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             recvbuf.data(),
                             2 * DIM,
                             mpi_traits<double>::type,
                             MPI_COMM_WORLD);
    SLIC_ASSERT(errf == MPI_SUCCESS);

    all_aabbs.clear();
    all_aabbs.reserve(m_nranks);
    for(int i = 0; i < m_nranks; ++i)
    {
      PointType lower(&recvbuf[i * 2 * DIM]);
      PointType upper(&recvbuf[i * 2 * DIM + DIM]);
      all_aabbs.emplace_back(BoxType(lower, upper));
    }
  }

  /// Compute bounding box for local part of a mesh.
  BoxType computeMeshBoundingBox(conduit::Node& mesh,
                                 const std::string& coordset) const
  {
    auto& coords = mesh[fmt::format("coordsets/{}/values", coordset)];
    SLIC_ASSERT( internal::extractDimension(coords) == NDIMS );
    const int npts = internal::extractSize(coords);
    BoxType rval;
    ArrayView<PointType> queryPts = ArrayView_from_Node<PointType>(mesh.fetch_existing("coordsets/coords/values/x"), npts);
    for(const auto &p : queryPts)
    {
      rval.addPoint(p);
    }
    return rval;
  }

  /**
   * \brief Computes the closest point within the objects for each query point
   * in the provided particle mesh, provided in the mesh blueprint rooted at \a query_mesh
   *
   * \param query_mesh The root node of a mesh blueprint for the query points
   * \param coordset The coordinate set for the query points
   *
   * Uses the \a coordset coordinate set of the provided blueprint mesh
   *
   * The particle mesh must contain the following fields:
   *   - cp_rank: Will hold the rank of the object point containing the closest point
   *   - cp_index: Will hold the index of the object point containing the closest point
   *   - closest_point: Will hold the position of the closest point
   *
   * \note The current implementation assumes that the coordinates and closest_points and contiguous
   * with stride NDIMS. We intend to loosen this restriction in the future
   *
   * \note We're temporarily also using a min_distance field while debugging this class.
   * The code will use this field if it is present in \a query_mesh.
   */
  void computeClosestPoints(conduit::Node& mesh_node,
                            const std::string& coordset) const
  {
    return new_computeClosestPoints(mesh_node, coordset);

    SLIC_ASSERT_MSG(
      isBVHTreeInitialized(),
      "BVH tree must be initialized before calling 'computeClosestPoints");

    // Utility function to dump a conduit node on each rank, e.g. for debugging
    auto dumpNode = [=](const conduit::Node& n,
                        const std::string&& fname,
                        const std::string& protocol = "json") {
      conduit::relay::io::save(n, fname, protocol);
    };

    // create conduit node containing data that has to xfer between ranks
    conduit::Node xfer_node;
    {
      // clang-format off
      auto& coords = mesh_node[fmt::format("coordsets/{}/values", coordset)];
      const int dim = internal::extractDimension(coords);
      const int npts = internal::extractSize(coords);
      SLIC_ASSERT( npts == mesh_node["fields/cp_rank/values"].dtype().number_of_elements() );

      xfer_node["npts"] = npts;
      xfer_node["dim"] = dim;
      xfer_node["src_rank"] = m_rank;
      xfer_node["coords"].set_external(internal::getPointer<double>(coords["x"]), dim * npts);
      xfer_node["cp_index"].set_external(internal::getPointer<axom::IndexType>(mesh_node["fields/cp_index/values"]), npts);
      xfer_node["cp_rank"].set_external(internal::getPointer<axom::IndexType>(mesh_node["fields/cp_rank/values"]), npts);
      xfer_node["closest_point"].set_external(internal::getPointer<double>(mesh_node["fields/closest_point/values/x"]), dim * npts);

      if(mesh_node.has_path("fields/min_distance"))
      {
        xfer_node["debug/min_distance"].set_external(internal::getPointer<double>(mesh_node["fields/min_distance/values"]), npts);
      }
      // clang-format on
    }

    if(m_isVerbose)
    {
      dumpNode(xfer_node, fmt::format("round_{}_r{}_begin.json", 0, m_rank));
    }

    // Find initial values on this rank
    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      computeLocalClosestPoints<SeqBVHTree>(m_bvh_seq.get(), xfer_node, true);
      break;

    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      computeLocalClosestPoints<OmpBVHTree>(m_bvh_omp.get(), xfer_node, true);
#endif
      break;

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      computeLocalClosestPoints<CudaBVHTree>(m_bvh_cuda.get(), xfer_node, true);
#endif
      break;
    }

    if(m_isVerbose)
    {
      dumpNode(xfer_node, fmt::format("round_{}_r{}_end.json", 0, m_rank));
    }

    if(m_nranks > 1)
    {
      // arbitrary tags for sending data to other ranks and getting it back
      const int tag_before = 1234;
      const int tag_after = 4321;

      // NOTE: uses a naive algorithm to computed distributed closest points
      // Every rank sends it data to every other rank.
      // TODO: Devise a more efficient algorithm to only send data to ranks with closer points
      for(int i = 1; i < m_nranks; ++i)
      {
        SLIC_INFO_IF(m_isVerbose,
                     fmt::format("=======  Round {}/{} =======", i, m_nranks));

        const int dst_rank = (m_rank + i) % m_nranks;
        const int rec_rank = (m_rank - i + m_nranks) % m_nranks;

        SLIC_INFO_IF(
          m_isVerbose,
          fmt::format("Rank {} -- sending to dst {}", m_rank, dst_rank));

        if(m_isVerbose)
        {
          dumpNode(xfer_node, fmt::format("round_{}_r{}_begin.json", i, m_rank));
        }

        // send and receive the query point data
        conduit::Node rec_node;
        internal::relay::mpi::send_and_recv_node(xfer_node,
                                                 rec_node,
                                                 dst_rank,
                                                 rec_rank,
                                                 tag_before,
                                                 MPI_COMM_WORLD);

        const int src_rank = rec_node["src_rank"].value();
        if(m_isVerbose)
        {
          dumpNode(
            rec_node,
            fmt::format("round_{}_r{}_comm_from_{}_A.json", i, m_rank, src_rank));
        }

        // compute the local data
        switch(m_runtimePolicy)
        {
        case RuntimePolicy::seq:
          computeLocalClosestPoints<SeqBVHTree>(m_bvh_seq.get(), rec_node, false);
          break;

        case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
          computeLocalClosestPoints<OmpBVHTree>(m_bvh_omp.get(), rec_node, false);
#endif
          break;

        case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
          computeLocalClosestPoints<CudaBVHTree>(m_bvh_cuda.get(), rec_node, false);
#endif
          break;
        }

        if(m_isVerbose)
        {
          dumpNode(
            rec_node,
            fmt::format("round_{}_r{}_comm_from_{}_B.json", i, m_rank, src_rank));
        }

        // update results
        conduit::Node proc_node;
        internal::relay::mpi::send_and_recv_node(rec_node,
                                                 proc_node,
                                                 src_rank,
                                                 dst_rank,
                                                 tag_after,
                                                 MPI_COMM_WORLD);

        if(m_isVerbose)
        {
          dumpNode(
            proc_node,
            fmt::format("round_{}_r{}_comm_from_{}_C.json", i, m_rank, dst_rank));
        }

        // copy data to mesh_node from proc_node
        const int npts = proc_node["npts"].value();
        SLIC_ASSERT_MSG(npts == xfer_node["npts"].as_int(),
                        fmt::format("{} vs {}",
                                    proc_node["npts"].as_int(),
                                    xfer_node["npts"].as_int()));
        // BTNG: These copies don't destroy existing data in xfer_node because they came from xfer_node, got updated remotely and sent back.
        axom::copy(xfer_node["cp_rank"].data_ptr(),
                   proc_node["cp_rank"].data_ptr(),
                   npts * sizeof(axom::IndexType));
        axom::copy(xfer_node["cp_index"].data_ptr(),
                   proc_node["cp_index"].data_ptr(),
                   npts * sizeof(axom::IndexType));
        axom::copy(xfer_node["closest_point"].data_ptr(),
                   proc_node["closest_point"].data_ptr(),
                   npts * sizeof(PointType));

        if(m_isVerbose)
        {
          dumpNode(mesh_node,
                   axom::fmt::format("round_{}_r{}_end.json", i, m_rank));

          SLIC_ASSERT_MSG(
            conduit::blueprint::mcarray::is_interleaved(
              mesh_node["fields/closest_point/values"]),
            fmt::format("After copy on iteration {}, 'closest_point' field of "
                        "'mesh_node' is not interleaved",
                        i));
        }

        MPI_Barrier(MPI_COMM_WORLD);
        slic::flushStreams();
      }
    }
  }
  void new_computeClosestPoints(conduit::Node& query_mesh,
                                const std::string& coordset) const
  {
    SLIC_ASSERT_MSG(
      isBVHTreeInitialized(),
      "BVH tree must be initialized before calling 'computeClosestPoints");

    // Utility function to dump a conduit node on each rank, e.g. for debugging
    auto dumpNode = [=](const conduit::Node& n,
                        const std::string&& fname,
                        const std::string& protocol = "json") {
      conduit::relay::io::save(n, fname, protocol);
    };

    BoxType queryPartitionBb = computeMeshBoundingBox(query_mesh, coordset);
std::cout << fmt::format("at {}, {} has objectBB {}, queryBB{}", __WHERE, m_rank, m_objectPartitionBbs[m_rank], queryPartitionBb) << std::endl;

    const int qmNpts = internal::extractSize(query_mesh[fmt::format("coordsets/{}/values", coordset)]);

    // create conduit node containing data that has to xfer between ranks
    conduit::Node xfer_node;
    {
      // clang-format off
      auto& coords = query_mesh.fetch_existing(fmt::format("coordsets/{}/values", coordset));
      const int dim = internal::extractDimension(coords);
      const int npts = internal::extractSize(coords);

      xfer_node["npts"] = npts;
      xfer_node["dim"] = dim;
      xfer_node["src_rank"] = m_rank;
      xfer_node["is_first"] = true;
      xfer_node["coords"].set_external(internal::getPointer<double>(coords["x"]), dim * npts);
      xfer_node["cp_index"].set_external(internal::getPointer<axom::IndexType>(query_mesh.fetch_existing("fields/cp_index/values")), npts);
      xfer_node["cp_rank"].set_external(internal::getPointer<axom::IndexType>(query_mesh.fetch_existing("fields/cp_rank/values")), npts);
      xfer_node["closest_point"].set_external(internal::getPointer<double>(query_mesh.fetch_existing("fields/closest_point/values/x")), dim * npts);
      put_to_conduit_node(queryPartitionBb, xfer_node["aabb"]);
      {
        auto cpIndexes =
          ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_index"), npts);
        auto cpRanks =
          ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_rank"), npts);
        for ( axom::IndexType ii=0; ii<cpRanks.size(); ++ii ) {
          cpRanks[ii] = -1;
          cpIndexes[ii] = -1;
        }
      }
checkXferNode(xfer_node);

      if(query_mesh.has_path("fields/min_distance"))
      {
        xfer_node["debug/min_distance"].set_external(internal::getPointer<double>(query_mesh["fields/min_distance/values"]), npts);
      }
      // clang-format on
    }

    if(m_isVerbose)
    {
      dumpNode(xfer_node, fmt::format("round_{}_r{}_begin.json", 0, m_rank));
    }

    // arbitrary tags for send/recv xfer_node.
    const int tag = 1234;

    std::vector<relay::mpi::ISendRequest> isendRequests(m_nranks);

    int toRecvCount = m_nranks - 1;
    // bool is_first = true;

    for(int round = 0; round < m_nranks || xfer_node.number_of_children() != 0; ++round)
    {
      // 1. If xfer_node isn't a skip message (it's real data, local or remote),
      //    forward it along the ring.
      // 2. Receive new xfer_node forwarded along the ring.
      // 3. If xfer_node is a skip message, go to next round.
      SLIC_INFO_IF(
        m_isVerbose,
        fmt::format("=======  Starting round {}/{} =======", round, m_nranks));

      // Empty xfer_node means the previous round received our own data back
      // or a skip message, both of which terminates progression on xfer_node.
      if(xfer_node.number_of_children() != 0 && m_nranks > 1)
      {
if (m_rank==1 && xfer_node.fetch_existing("src_rank").as_int() == 2) {
std::cout << fmt::format("at {}, {} has {}'s data in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), round) << std::endl;
}
        get_from_conduit_node(queryPartitionBb, xfer_node.fetch_existing("aabb"));

        // Send the current xfer_node so we can receive the next one.
        int next_dst = (m_rank + 1) % m_nranks;
        get_from_conduit_node(queryPartitionBb, xfer_node["aabb"]);
        while (xfer_node.number_of_children() != 0 && next_dst != m_rank)
        {
          if(m_isVerbose)
          {
            dumpNode(xfer_node,
                     fmt::format("round_{}_r{}_begin.json", round, m_rank));
          }

          double sqDist =
            squared_distance(m_objectPartitionBbs[next_dst], queryPartitionBb);
std::cout << fmt::format("at {}, {} has sqDist {:.4f} between {}'s and object partition {} in round {}", __WHERE, m_rank, sqDist, xfer_node.fetch_existing("src_rank").as_int(), next_dst, round) << std::endl;

          if(next_dst != m_rank)
          {
            SLIC_INFO_IF(
              m_isVerbose,
              fmt::format("Rank {} -- sending to dst {}", m_rank, next_dst));

            // TODO: We can also track the farthest distance partition is from
            // found points.  If this is closer than sqDist, there's no point
            // sending to next_dst.  This could save communication even without
            // a user-specified threshold.
            if(sqDist <= m_sqDistanceThreshold ||
               next_dst == xfer_node.fetch_existing("src_rank").as_int())
            {
              if(xfer_node.fetch_existing("src_rank").as_int() == m_rank) ++toRecvCount;
              xfer_node["sender_rank"] = m_rank;
std::cout << fmt::format("at {}, {} sending {}'s data to {} in round {}", __WHERE, xfer_node.fetch_existing("sender_rank").as_int(), xfer_node.fetch_existing("src_rank").as_int(), next_dst, round) << std::endl;
checkXferNode(xfer_node);
              isend_using_schema(xfer_node,
                                 next_dst,
                                 tag,
                                 MPI_COMM_WORLD,
                                 &isendRequests[next_dst]);
              assert(xfer_node.number_of_children() > 0);
              xfer_node.reset();
              assert(xfer_node.number_of_children() == 0);
              break;
            }
            else
            {
              conduit::Node skip;
              // I think we can use an empty Node, but be explicit for now.
              skip["skip"] = true;
              skip["src_rank"] = xfer_node.fetch_existing("src_rank");
              skip["sender_rank"] = m_rank;
std::cout << fmt::format("at {}, {} sending {}'s skip to {} in round {}", __WHERE, skip["sender_rank"].as_int(), skip["src_rank"].as_int(), next_dst, round) << std::endl;
              isend_using_schema(skip,
                                 next_dst,
                                 tag,
                                 MPI_COMM_WORLD,
                                 &isendRequests[next_dst]);
            }
          }

          next_dst = (next_dst + 1) % m_nranks;
        } // next_dst loop
      } // block to send xfer_node
      else
      {
std::cout << fmt::format("at {}, {} has nothing in round {}", __WHERE, m_rank, round) << std::endl;
      }

      if(xfer_node.number_of_children() == 0 && toRecvCount > 0) // if(m_nranks > 1)
      {
        // TODO: What if the local query mesh is still in xfer_node because it skipped the entire ring?
        conduit::relay::mpi::recv_using_schema(xfer_node,
                                               MPI_ANY_SOURCE,
                                               tag,
                                               MPI_COMM_WORLD);
        --toRecvCount;
std::cout << fmt::format("at {}, {} received {}'s data from {} in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), xfer_node.fetch_existing("sender_rank").as_int(), round) << std::endl;
        if(xfer_node.has_path("skip"))
        {
std::cout << fmt::format("at {}, {} skipping {}'s data from {} in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), xfer_node.fetch_existing("sender_rank").as_int(), round) << std::endl;
          xfer_node.reset();
        }
else {
  checkXferNode(xfer_node);
}
      }

      if(xfer_node.number_of_children() != 0)
      {
if (xfer_node.has_path("sender_rank"))
  std::cout << fmt::format("at {}, {} processing {}'s data from {} in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), xfer_node.fetch_existing("sender_rank").as_int(), round) << std::endl;
else
  std::cout << fmt::format("at {}, {} processing {}'s data in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), round) << std::endl;

        // Distance search using local object partition and the xfer_node.
        switch(m_runtimePolicy)
        {
        case RuntimePolicy::seq:
          computeLocalClosestPoints<SeqBVHTree>(m_bvh_seq.get(), xfer_node, xfer_node.has_path("is_first"));
          break;

        case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
          computeLocalClosestPoints<OmpBVHTree>(m_bvh_omp.get(), xfer_node, xfer_node.has_path("is_first"));
#endif
          break;

        case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
          computeLocalClosestPoints<CudaBVHTree>(m_bvh_cuda.get(), xfer_node, xfer_node.has_path("is_first"));
#endif
          break;
        }
if(xfer_node.has_path("is_first")) xfer_node.remove("is_first"); // is_first = false;

        if(xfer_node.fetch_existing("src_rank").as_int() == m_rank)
        {
if (xfer_node.has_path("sender_rank"))
  std::cout << fmt::format("at {}, {} saving {}'s data from {} in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), xfer_node.fetch_existing("sender_rank").as_int(), round) << std::endl;
else
  std::cout << fmt::format("at {}, {} saving {}'s data in round {}", __WHERE, m_rank, xfer_node.fetch_existing("src_rank").as_int(), round) << std::endl;
          // This is the query partition we started with.
          // Copy it back to query_mesh.
          const int npts = xfer_node["npts"].value();
          assert( npts == qmNpts );
          assert( npts == internal::extractSize(query_mesh[fmt::format("coordsets/{}/values", coordset)]) );
          auto& qmcpr = query_mesh.fetch_existing("fields/cp_rank/values");
          auto& qmcpi = query_mesh.fetch_existing("fields/cp_index/values");
          auto& qmcpcp = query_mesh.fetch_existing("fields/closest_point/values/x");
          // query_mesh doesn't have npts.
          if(xfer_node.fetch_existing("cp_rank").data_ptr() != qmcpr.data_ptr())
          {
            axom::copy(qmcpr.data_ptr(),
                       xfer_node.fetch_existing("cp_rank").data_ptr(),
                       npts * sizeof(axom::IndexType));
          }
          if(xfer_node.fetch_existing("cp_index").data_ptr() != qmcpi.data_ptr())
          {
            axom::copy(qmcpi.data_ptr(),
                       xfer_node.fetch_existing("cp_index").data_ptr(),
                       npts * sizeof(axom::IndexType));
          }
          if(xfer_node.fetch_existing("closest_point").data_ptr() != qmcpcp.data_ptr())
          {
            axom::copy(qmcpcp.data_ptr(),
                       xfer_node.fetch_existing("closest_point").data_ptr(),
                       npts * sizeof(PointType));
          }

          if(m_isVerbose)
          {
            dumpNode(query_mesh,
                     axom::fmt::format("round_{}_r{}_end.json", round, m_rank));

            SLIC_ASSERT_MSG(
              conduit::blueprint::mcarray::is_interleaved(
                query_mesh["fields/closest_point/values"]),
              fmt::format(
                "After copy on iteration {}, 'closest_point' field of "
                "'query_mesh' is not interleaved",
                round));
          }
          // Take xfer_node out of circulation,
          // because it has completed its trip through the ring.
          xfer_node.reset();
        }
      } // Locally process xfer_node

      SLIC_INFO_IF(
        m_isVerbose,
        fmt::format("=======  End of round {}/{} =======", round, m_nranks));
    }  // round loop

    // Complete non-blocking sends.
    std::vector<MPI_Request> allReqs;
    allReqs.reserve(isendRequests.size());
    for(const auto& isr : isendRequests)
    {
      allReqs.push_back(isr.mpi_request);
    }
    std::vector<MPI_Status> allStats(isendRequests.size());
    MPI_Waitall(int(allReqs.size()), allReqs.data(), allStats.data());

    MPI_Barrier(MPI_COMM_WORLD);
    slic::flushStreams();
  }

private:
  /// Sets the allocator ID to the default associated with the execution policy
  void setAllocatorID()
  {
    // This function uses the default allocator ID for the execution space
    // TODO: Add overload to allow the user to set an allocator ID

    switch(m_runtimePolicy)
    {
    case RuntimePolicy::seq:
      m_allocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
      break;
    case RuntimePolicy::omp:
#ifdef _AXOM_DCP_USE_OPENMP
      m_allocatorID = axom::execution_space<axom::OMP_EXEC>::allocatorID();
#endif
      break;

    case RuntimePolicy::cuda:
#ifdef _AXOM_DCP_USE_CUDA
      m_allocatorID = axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
#endif
      break;
    }
  }

  /**
   * \brief Extracts a field \a fieldName from the mesh blueprint
   *
   * \tparam T The type for the underlying array
   * \param mesh_node The conduit node at the root of the mesh blueprint
   * \param field_name The name of the field
   * \param field_template Template string for the path to the field
   * \param num_points The size of the field
   * \return An arrayview over the field data
   */
  template <typename T>
  axom::ArrayView<T> extractField(conduit::Node& mesh_node,
                                  std::string&& field_name,
                                  std::string&& path_template,
                                  int num_points) const
  {
    const std::string path = axom::fmt::format(path_template, field_name);
    SLIC_ASSERT_MSG(
      mesh_node.has_path(path),
      fmt::format(
        "Input to `computeClosestPoint()` must have a field named `{}`",
        field_name));

    return internal::ArrayView_from_Node<T>(mesh_node[path], num_points);
  }

  // Note: following should be private, but nvcc complains about lambdas in private scope
public:
  /// Templated implementation of generateBVHTree function
  template <typename BVHTreeType>
  bool generateBVHTreeImpl(BVHTreeType* bvh)
  {
    using ExecSpace = typename BVHTreeType::ExecSpaceType;

    SLIC_ASSERT(bvh != nullptr);

    const int npts = m_points.size();
    axom::Array<BoxType> boxesArray(npts, npts, m_allocatorID);
    auto boxesView = boxesArray.view();

    /// GOT TO HERE -- fix for templated ExecSpace!
    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { boxesView[i] = BoxType {m_points[i]}; });

    // Build bounding volume hierarchy
    bvh->setAllocatorID(m_allocatorID);
    int result = bvh->initialize(boxesView, npts);

    gatherBVHRoots();

    return (result == spin::BVH_BUILD_OK);
  }

  void checkXferNode(conduit::Node &xfer_node) const
  {
    const int npts = xfer_node.fetch_existing("npts").value();
    auto cpIndexes =
      ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_index"), npts);
    auto cpRanks =
      ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_rank"), npts);
    assert(cpIndexes.size() == cpRanks.size());
    for ( axom::IndexType ii=0; ii<cpRanks.size(); ++ii ) {
      assert((cpRanks[ii] == -1) == (cpIndexes[ii] == -1)); // DBG
    }

    auto cp_idx = axom::Array<axom::IndexType>(cpIndexes, m_allocatorID);
    auto cp_ranks = axom::Array<axom::IndexType>(cpRanks, m_allocatorID);

    auto query_inds = cp_idx.view();
    auto query_ranks = cp_ranks.view();
    for ( axom::IndexType ii=0; ii<query_inds.size(); ++ii ) {
      assert((query_ranks[ii] == -1) == (query_inds[ii] == -1)); // DBG
    }
  }

  /**
   */
  template <typename BVHTreeType>
  void computeLocalClosestPoints(const BVHTreeType* bvh,
                                 conduit::Node& xfer_node,
                                 bool is_first) const
  {
    using ExecSpace = typename BVHTreeType::ExecSpaceType;
    using axom::primal::squared_distance;
    using int32 = axom::int32;

    // Extract the dimension and number of points from the coordinate values group
    const int dim = xfer_node.fetch_existing("dim").value();
    const int npts = xfer_node.fetch_existing("npts").value();
    SLIC_ASSERT(dim == NDIMS);

    /// Extract fields from the input node as ArrayViews
    auto queryPts = ArrayView_from_Node<PointType>(xfer_node.fetch_existing("coords"), npts);
    auto cpIndexes =
      ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_index"), npts);
    auto cpRanks =
      ArrayView_from_Node<axom::IndexType>(xfer_node.fetch_existing("cp_rank"), npts);
    auto closestPts =
      ArrayView_from_Node<PointType>(xfer_node.fetch_existing("closest_point"), npts);

    /// Create ArrayViews in ExecSpace that are compatible with fields
    // This deep-copies host memory in xfer_node to device memory.
    // TODO: Avoid copying arrays (here and at the end) if both are on the host
    // BTNG: This is probably a deep copy.
    // xfer_node["cp_index"]           -view> cpIndexes  -devicecopy> cp_idx   -view> query_inds     -devicecopy> cpIndexes
    // xfer_node["cp_rank"]            -view> cpRanks    -devicecopy> cp_ranks -view> query_ranks    -devicecopy> cpRanks
    // xfer_node["closest_point"]      -view> closestPts -devicecopy> cp_pos   -view> query_pos      -devicecopy> closestPts
    // xfer_node["debug/min_distance"] -view> minDist    -devicecopy> cp_dist  -view> query_min_dist -devicecopy> minDist
    auto cp_idx = is_first
      ? axom::Array<axom::IndexType>(npts, npts, m_allocatorID)
      : axom::Array<axom::IndexType>(cpIndexes, m_allocatorID);
    auto cp_ranks = is_first
      ? axom::Array<axom::IndexType>(npts, npts, m_allocatorID)
      : axom::Array<axom::IndexType>(cpRanks, m_allocatorID);

    /// PROBLEM: The striding does not appear to be retained by conduit relay
    ///          We might need to transform it? or to use a single array w/ pointers into it?
    auto cp_pos = is_first ? axom::Array<PointType>(npts, npts, m_allocatorID)
                           : axom::Array<PointType>(closestPts, m_allocatorID);

    // DEBUG
    const bool has_min_distance = xfer_node.has_path("debug/min_distance");
    auto minDist = has_min_distance
      ? ArrayView_from_Node<double>(xfer_node.fetch_existing("debug/min_distance"), npts)
      : ArrayView<double>();

    auto cp_dist = has_min_distance
      ? (is_first ? axom::Array<double>(npts, npts, m_allocatorID)
                  : axom::Array<double>(minDist, m_allocatorID))
      : axom::Array<double>(0, 0, m_allocatorID);
    auto query_min_dist = cp_dist.view();
    // END DEBUG

    if(is_first)
    {
      cp_ranks.fill(-1);
      cp_idx.fill(-1);
    }
    auto query_inds = cp_idx.view();
    auto query_ranks = cp_ranks.view();
    auto query_pos = cp_pos.view();

    /// Create an ArrayView in ExecSpace that is compatible with queryPts
    PointArray execPoints(queryPts, m_allocatorID);
    auto query_pts = execPoints.view();

    // Get a device-useable iterator
    auto it = bvh->getTraverser();
    const int rank = m_rank;

          // AXOM_PERF_MARK_SECTION(
            // "ComputeClosestPoints",
      axom::for_all<ExecSpace>(
        npts,
        AXOM_LAMBDA(int32 idx) mutable {
          assert((query_ranks[idx] == -1) == (query_inds[idx] == -1)); // DBG
          PointType qpt = query_pts[idx];

          MinCandidate curr_min {};
          if(query_ranks[idx] >= 0)  // i.e. we've already found a candidate closest
          {
            curr_min.minSqDist = squared_distance(qpt, query_pos[idx]);
            curr_min.minElem = query_inds[idx];
            curr_min.minRank = query_ranks[idx];
          }

          auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
            const int candidate_idx = leaf_nodes[current_node];
            const PointType candidate_pt = m_points[candidate_idx];
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
            return sqDist <= curr_min.minSqDist &&
              sqDist <= m_sqDistanceThreshold;
          };

          // Traverse the tree, searching for the point with minimum distance.
          it.traverse_tree(qpt, searchMinDist, traversePredicate);

          // If modified, update the fields that changed
          if(curr_min.minRank == rank)
          {
            // SLIC_INFO(
            //   fmt::format("[On rank {}] Updating dist for particle {} ::"
            //               "{{ cp: {}->{}; "
            //               "dist: {}->{}; "
            //               "rank: {}->{}; "
            //               "elem: {}->{}}}",
            //               rank,
            //               idx,
            //               query_pos[idx],             //cp
            //               m_points[curr_min.minElem],
            //               query_min_dist[idx],        // dist
            //               sqrt(curr_min.minSqDist),
            //               query_ranks[idx],           // rank
            //               curr_min.minRank,
            //               query_inds[idx],            // index
            //               curr_min.minElem));

            query_inds[idx] = curr_min.minElem;
            query_ranks[idx] = curr_min.minRank;
            query_pos[idx] = m_points[curr_min.minElem];

            //DEBUG
            if(has_min_distance)
            {
              query_min_dist[idx] = sqrt(curr_min.minSqDist);
            }
          }
          // else
          // {
          //   SLIC_INFO(
          //     fmt::format("[On rank {}] Didn't update dist for particle {} ::"
          //                 "{{ cp: {}; "
          //                 "dist: {}; "
          //                 "rank: {}; "
          //                 "elem: {}}}",
          //                 rank,
          //                 idx,
          //                 query_pos[idx],
          //                 query_min_dist[idx],
          //                 query_ranks[idx],
          //                 query_inds[idx]));
          // }
        }
                               );
      // );

    // BTNG: copy args are (dst,src,numbytes)
    axom::copy(cpIndexes.data(),
               query_inds.data(),
               cpIndexes.size() * sizeof(axom::IndexType));
    axom::copy(cpRanks.data(),
               query_ranks.data(),
               cpRanks.size() * sizeof(axom::IndexType));
    axom::copy(closestPts.data(),
               query_pos.data(),
               closestPts.size() * sizeof(PointType));

    // DEBUG
    if(has_min_distance)
    {
      axom::copy(minDist.data(),
                 query_min_dist.data(),
                 minDist.size() * sizeof(double));
    }
  }

private:
  RuntimePolicy m_runtimePolicy;
  bool m_isVerbose {false};
  double m_sqDistanceThreshold;
  int m_allocatorID;
  int m_rank;
  int m_nranks;

  PointArray m_points;
  BoxArray m_objectPartitionBbs;

  std::unique_ptr<SeqBVHTree> m_bvh_seq;

#ifdef _AXOM_DCP_USE_OPENMP
  std::unique_ptr<OmpBVHTree> m_bvh_omp;
#endif

#ifdef _AXOM_DCP_USE_CUDA
  std::unique_ptr<CudaBVHTree> m_bvh_cuda;
#endif
};

}  // namespace internal

/**
 * \brief Encapsulated the Distributed closest point query for a collection of query points
 * over an "object mesh"
 *
 * The object mesh and the query mesh are provided as conduit nodes using the mesh blueprint schema.
 * Each of these are distributed over the problem's mpi ranks. This class orchestrates passing
 * the query points to all ranks whose object meshes might contain a closest point.
 *
 * \note This class assumes that the object mesh and the query mesh are decomposed over the same
 * number of mpi ranks.
 * \warning This class will need to support cases where some ranks have zero object points.
 * This is not currently supported.
 *
 * \note The class currently supports object meshes that are comprised of a collection of points.
 * In the future, we'd like to consider more general object meshes, e.g. triangle meshes.
 *
 * \note The class currently supports object meshes and query meshes with a single domain per MPI rank.
 * We intend to add support for multiple computational domains on each rank.
 *
 * To use this class, first set some parameters, such as the runtime execution policy,
 * then pass in the object mesh and build a spatial index over this mesh.
 * Finally, compute the closest points in the object mesh to each point in a query mesh
 * using the \a computeClosestPoint() function.
 *
 * \note The implementation currently assumes that the coordinates for the positions and vector field
 * data are interleaved (i.e. xyzxyzxyz....). We will relax this assumption in the future to support both
 * interleaved and strided data.
 */
class DistributedClosestPoint
{
public:
  using RuntimePolicy = DistributedClosestPointRuntimePolicy;

public:
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
    }

    return false;
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
   * \param [in] mesh_node Conduit node for the object mesh
   * \param [in] coordset The name of the coordset for the object mesh's coordinates
   *
   * \pre \a mesh_node must follow the mesh blueprint convention
   * \pre Dimension of the mesh must be 2D or 3D
   */
  void setObjectMesh(const conduit::Node& mesh_node, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(this->isValidBlueprint(mesh_node));

    // Extract the dimension and number of points from the coordinate values group
    auto valuesPath = fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(mesh_node.has_path(valuesPath));
    auto& values = mesh_node[valuesPath];

    const int dim = internal::extractDimension(values);
    setDimension(dim);

    allocateQueryInstance();

    const int N = internal::extractSize(values);

    // dispatch to implementation class over dimension
    switch(m_dimension)
    {
    case 2:
      m_dcp_2->importObjectPoints(values, N);
      break;
    case 3:
      m_dcp_3->importObjectPoints(values, N);
      break;
    }
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
  void computeClosestPoints(conduit::Node& query_node,
                            const std::string& cooordset)
  {
    SLIC_ASSERT_MSG(m_objectMeshCreated,
                    "Must call 'setObjectMesh' before calling generateBVHTree");

    SLIC_ASSERT(this->isValidBlueprint(query_node));

    // dispatch to implementation class over dimension
    switch(m_dimension)
    {
    case 2:
      m_dcp_2->setSquaredDistanceThreshold(m_sqDistanceThreshold);
      m_dcp_2->computeClosestPoints(query_node, cooordset);
      break;
    case 3:
      m_dcp_3->setSquaredDistanceThreshold(m_sqDistanceThreshold);
      m_dcp_3->computeClosestPoints(query_node, cooordset);
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
      m_dcp_2 = std::unique_ptr<internal::DistributedClosestPointImpl<2>>(
        new internal::DistributedClosestPointImpl<2>(m_runtimePolicy,
                                                     m_isVerbose));
      m_objectMeshCreated = true;
      break;
    case 3:
      m_dcp_3 = std::unique_ptr<internal::DistributedClosestPointImpl<3>>(
        new internal::DistributedClosestPointImpl<3>(m_runtimePolicy,
                                                     m_isVerbose));
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
    if(!conduit::blueprint::mpi::verify("mesh", mesh_node, info, MPI_COMM_WORLD))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }

    return success;
  }

private:
  RuntimePolicy m_runtimePolicy {RuntimePolicy::seq};
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

#endif  //  QUEST_DISTRIBUTED_CLOSEST_POINT_H_
