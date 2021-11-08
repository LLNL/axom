// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SIGNED_DISTANCE_HPP_
#define AXOM_QUEST_SIGNED_DISTANCE_HPP_

// axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/utilities/Utilities.hpp"

// primal includes
#include "axom/spin/BVH.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/utils/ZipPoint.hpp"

// mint includes
#include "axom/mint/config.hpp"
#include "axom/mint/mesh/Field.hpp"
#include "axom/mint/mesh/FieldData.hpp"
#include "axom/mint/mesh/FieldVariable.hpp"
#include "axom/mint/mesh/Mesh.hpp"

// C/C++ includes
#include <cmath>  // for std::sqrt()

namespace axom
{
namespace quest
{
namespace detail
{
/*!
 * \brief Helper class for handling unstructured surface meshes on the GPU.
 */
struct UcdMeshData
{
  int shape_type;
  mint::CellType single_cell_type;
  const mint::CellType* cell_types;
  const IndexType* cells_to_nodes;
  IndexType nodes_per_cell;
  const IndexType* cell_node_offsets;

  /*!
   * \brief Returns the type of the cell at the given index.
   */
  AXOM_HOST_DEVICE mint::CellType getCellType(IndexType cellId) const;

  /*!
   * \brief Returns the node IDs of a given cell index.
   *
   * \param [in] cellId the cell index to query
   * \param [out] nnodes the number of nodes in the given cell
   *
   * \return pointer to the node IDs in the underlying cell node array.
   */
  AXOM_HOST_DEVICE const IndexType* getCellNodeIDs(IndexType cellId,
                                                   int& nnodes) const;
};

/*!
 * \brief Creates a POD UcdMeshData object from a polymorphic mint::Mesh
 *  pointer.
 *
 * \param [in] surfaceMesh the mint::Mesh pointer with the surface mesh
 * \param [out] outSurfData the UcdMeshData object with underlying mesh data
 *
 * \return true if the mesh was an unstructured mesh, false otherwise
 */
bool SD_GetUcdMeshData(const mint::Mesh* surfaceMesh, UcdMeshData& outSurfData);

/// Enum for different 'location' types returned by primal::closest_point()
enum class ClosestPointLocType
{
  uninitialized = -1,
  vertex = 0,
  edge = 1,
  face = 2
};

/// Converts from \a loc code returned by primal::closest_point() to a \a ClosestPointLocType enum
inline ClosestPointLocType getClosestPointLocType(int loc)
{
  SLIC_ASSERT_MSG(loc >= -3,
                  "Invalid closest point location type: "
                    << loc << ". See documentation for primal::closest_point().");
  switch(loc)
  {
  case -3:  // intentional fall-through
  case -2:
  case -1:
    return ClosestPointLocType::edge;
  case 0:  // intentional fall-through
  case 1:
  case 2:
    return ClosestPointLocType::vertex;
  default:
    return ClosestPointLocType::face;
  }
}

/// A \a ClosestPointLocType is shared for a vertex or edge, unshared otherwise
inline bool isClosestPointTypeShared(ClosestPointLocType cpt)
{
  return cpt == ClosestPointLocType::edge || cpt == ClosestPointLocType::vertex;
}

}  // end namespace detail

template <int NDIMS, typename ExecSpace = axom::SEQ_EXEC>
class SignedDistance
{
public:
  using PointType = axom::primal::Point<double, NDIMS>;
  using VectorType = axom::primal::Vector<double, NDIMS>;
  using TriangleType = axom::primal::Triangle<double, NDIMS>;
  using BoxType = axom::primal::BoundingBox<double, NDIMS>;
  using ZipPoint = axom::primal::ZipIndexable<PointType>;
  using BVHTreeType = axom::spin::BVH<NDIMS, ExecSpace>;

private:
  struct MinCandidate
  {
    /// Squared distance to query point
    double minSqDist {numerics::floating_point_limits<double>::max()};
    /// Closest point on element
    PointType minPt {};
    /// Type of the closest point
    detail::ClosestPointLocType minType {
      detail::ClosestPointLocType::uninitialized};
    /// Index within mesh of closest element
    int minElem;
    /// The actual cell element
    TriangleType minTri;

    /// The normal, if the closest point is on an edge
    VectorType sumNormals {};
    /// The normal, if the closest point is a node
    VectorType sumNormalsAngWt {};
  };

public:
  /*!
   * \brief Creates a SignedDistance instance for queries on the given mesh.
   * \param [in] surfaceMesh user-supplied surface mesh.
   * \param [in] isWatertight indicates if the surface mesh is closed.
   * \param [in] computeSign indicates if distance queries should compute signs (optional).
   * \param [in] allocatorID the allocator to create the underlying BVH with (optional).
   *
   * \note computeSign defaults to \a true when not specified.
   *
   * \note The given surface mesh must be allocated in a memory space
   *  compatible with the execution space specified in the SignedDistance
   *  instantiation.
   *
   * \note The given allocatorID must be compatible with the execution space
   *  specified in the SignedDistance instantiation. By default, a compatible
   *  allocator ID used if allocatorID is not specified.
   *
   * \pre surfaceMesh != nullptr
   */
  SignedDistance(const mint::Mesh* surfaceMesh,
                 bool isWatertight,
                 bool computeSign = true,
                 int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Reinitializes a SignedDistance instance with a new surface mesh.
   *
   * \param [in] surfaceMesh user-supplied surface mesh.
   * \param [in] allocatorID the allocator to create the underlying BVH with (optional).
   *
   * \note The given surface mesh must be allocated in a memory space
   *  compatible with the execution space specified in the SignedDistance
   *  instantiation.
   *
   * \note The given allocatorID must be compatible with the execution space
   *  specified in the SignedDistance instantiation. By default, a compatible
   *  allocator ID used if allocatorID is not specified.
   *
   * \return status true if reinitialization was successful.
   * \pre surfaceMesh != nullptr
   */
  bool setMesh(const mint::Mesh* surfaceMesh,
               int allocatorID = axom::execution_space<ExecSpace>::allocatorID());

  /*!
   * \brief Computes the distance of the given point to the input surface mesh.
   *
   * \param [in] x x-coordinate of the query point
   * \param [in] y y-coordinate of the query point
   * \param [in] z z-coordinate of the query point
   *
   * \note When the input is not a closed surface mesh, the assumption is that
   *  the surface mesh divides the computational mesh domain into two regions.
   *  Hence, the surface mesh has to span the entire domain of interest, e.g.,
   *  the computational mesh at which the signed distance field is evaluated,
   *  along some plane.
   *
   * \warning The sign of the distance from a given query point is determined by
   *  a pseudo-normal which is computed at the closest point on the surface
   *  mesh. For a non-watertight mesh, the sign of the distance is not defined
   *  everywhere. Specifically, the sign is ambiguous for all points for which
   *  a normal projection onto the surface does not exist.
   *
   * \return minDist minimum signed distance of the query point to the surface.
   */
  double computeDistance(double x, double y, double z = 0.0)
  {
    PointType pt = PointType::make_point(x, y, z);
    return (computeDistance(pt));
  }

  /*!
   * \brief Computes the distance of the given point to the surface mesh.
   * \param [in] queryPnt user-supplied point.
   *
   * \note When the input is not a closed surface mesh, the assumption is that
   *  the surface mesh divides the computational mesh domain into two regions.
   *  Hence, the surface mesh has to span the entire domain of interest, e.g.,
   *  the computational mesh at which the signed distance field is evaluated,
   *  along some plane.
   *
   * \warning The sign of the distance from a given query point is determined by
   *  a pseudo-normal which is computed at the closest point on the surface
   *  mesh. For a non-watertight mesh, the sign of the distance is not defined
   *  everywhere. Specifically, the sign is ambiguous for all points for which
   *  a normal projection onto the surface does not exist.
   *
   * \return minDist the signed minimum distance to the surface mesh.
   */
  double computeDistance(const PointType& queryPnt) const;

  /*!
   * \brief Computes the distances of a set of points to the surface mesh.
   * \param [in] npts number of points to query
   * \param [in] queryPnt user-supplied point indexable type. This can be a
   *  pointer-to-array, or a ZipIndexable<PointType>.
   * \param [out] outSgnDist array to fill with corresponding signed distances
   *  for query points
   * \param [out] outClosestPts array to fill with closest points on the mesh.
   *  Optional.
   *
   * \note When the input is not a closed surface mesh, the assumption is that
   *  the surface mesh divides the computational mesh domain into two regions.
   *  Hence, the surface mesh has to span the entire domain of interest, e.g.,
   *  the computational mesh at which the signed distance field is evaluated,
   *  along some plane.
   *
   * \warning The sign of the distance from a given query point is determined by
   *  a pseudo-normal which is computed at the closest point on the surface
   *  mesh. For a non-watertight mesh, the sign of the distance is not defined
   *  everywhere. Specifically, the sign is ambiguous for all points for which
   *  a normal projection onto the surface does not exist.
   *
   * \return minDist the signed minimum distance to the surface mesh.
   *
   * \pre outSgnDist != nullptr
   */
  template <typename PointIndexable>
  void computeDistances(int npts,
                        PointIndexable queryPts,
                        double* outSgnDist,
                        PointType* outClosestPts = nullptr) const;

  /*!
   * \brief Returns a const reference to the underlying bucket tree.
   * \return ptr pointer to the underlying bucket tree
   */
  const BVHTreeType& getBVHTree() const { return m_bvh; }

private:
  /*!
   * \brief Computes the bounding box of the given cell on the surface mesh.
   * \param [in] icell the index of the cell on the surface mesh.
   * \param [in] mesh the surface mesh data
   * \param [in] meshPts the surface mesh point coordinate data
   * \return box bounding box of the cell.
   * \pre icell >= 0 && icell < m_surfaceMesh->getNumberOfCells()
   */
  AXOM_HOST_DEVICE static BoxType getCellBoundingBox(axom::IndexType icell,
                                                     const detail::UcdMeshData& mesh,
                                                     ZipPoint meshPts);

  /*!
   * \brief Checks a given candidate surface element against a query point
   *  and updates the minimum-distance candidate data if the element is closer.
   *
   * \param [in] qpt query point to check against surface element
   * \param [in,out] currMin the minimum-distance candidate data to update
   * \param [in] cellId the surface element ID to evaluate
   * \param [in] mesh the surface mesh data
   * \param [in] meshPts the surface mesh point coordinate data
   * \param [in] computeSign if true, will compute normals for the minimum
   *  candidate to use in determining sign
   */
  AXOM_HOST_DEVICE static void checkCandidate(const PointType& qpt,
                                              MinCandidate& currMin,
                                              IndexType cellId,
                                              const detail::UcdMeshData& mesh,
                                              ZipPoint meshPts,
                                              bool computeSign);

  /*!
   * \brief Computes the sign of the given query point given the closest point data
   *
   * \param [in] qpt query point to check against surface element
   * \param [in] currMin the minimum-distance surface element data
   *
   * \return sgn 1.0 if outside, -1.0 if inside
   */
  AXOM_HOST_DEVICE static double computeSign(const PointType& qpt,
                                             const MinCandidate& currMin);

private:
  bool m_isInputWatertight;        /*!< indicates if input is watertight     */
  bool m_computeSign;              /*!< indicates if queries compute sign    */
  const mint::Mesh* m_surfaceMesh; /*!< User-supplied surface mesh.          */
  BoxType m_boxDomain;             /*!< bounding box containing surface mesh */
  BVHTreeType m_bvh;               /*!< Spatial acceleration data-structure. */

  DISABLE_COPY_AND_ASSIGNMENT(SignedDistance);
};

}  // end namespace quest
}  // end namespace axom

//------------------------------------------------------------------------------
//           SignedDistance Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace quest
{
//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
SignedDistance<NDIMS, ExecSpace>::SignedDistance(const mint::Mesh* surfaceMesh,
                                                 bool isWatertight,
                                                 bool computeSign,
                                                 int allocatorID)
  : m_isInputWatertight(isWatertight)
  , m_computeSign(computeSign)
{
  // Sanity checks
  SLIC_ASSERT(surfaceMesh != nullptr);

  bool bvh_constructed = setMesh(surfaceMesh, allocatorID);
  SLIC_ASSERT(bvh_constructed);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
bool SignedDistance<NDIMS, ExecSpace>::setMesh(const mint::Mesh* surfaceMesh,
                                               int allocatorID)
{
  AXOM_PERF_MARK_FUNCTION("SignedDistance::setMesh");
  SLIC_ASSERT(surfaceMesh != nullptr);

  m_surfaceMesh = surfaceMesh;
  const axom::IndexType ncells = m_surfaceMesh->getNumberOfCells();
  const axom::IndexType nnodes = m_surfaceMesh->getNumberOfNodes();

  // Get device-usable mesh data
  const double* xs = m_surfaceMesh->getCoordinateArray(0);
  const double* ys = m_surfaceMesh->getCoordinateArray(1);
  const double* zs = nullptr;
  if(NDIMS == 3)
  {
    zs = m_surfaceMesh->getCoordinateArray(2);
  }

  ZipPoint surfPts {{xs, ys, zs}};

  detail::UcdMeshData surfaceData;
  bool mesh_valid = detail::SD_GetUcdMeshData(m_surfaceMesh, surfaceData);
  AXOM_UNUSED_VAR(mesh_valid);
  SLIC_CHECK_MSG(mesh_valid, "Input mesh is not an unstructured surface mesh");

  // compute bounding box of surface mesh
  // NOTE: this should be changed to an oriented bounding box in the future.
#ifdef AXOM_USE_RAJA
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

  double minInit = numerics::floating_point_limits<double>::max();
  double maxInit = numerics::floating_point_limits<double>::lowest();
  RAJA::ReduceMin<reduce_policy, double> xmin(minInit), ymin(minInit),
    zmin(minInit);

  RAJA::ReduceMax<reduce_policy, double> xmax(maxInit), ymax(maxInit),
    zmax(maxInit);

  for_all<ExecSpace>(
    nnodes,
    AXOM_LAMBDA(axom::IndexType inode) {
      xmin.min(xs[inode]);
      xmax.max(xs[inode]);

      ymin.min(ys[inode]);
      ymax.max(ys[inode]);

      if(NDIMS == 3)
      {
        zmin.min(zs[inode]);
        zmax.max(zs[inode]);
      }
    });
  PointType boxMin {xmin.get(), ymin.get(), zmin.get()};
  PointType boxMax {xmax.get(), ymax.get(), zmax.get()};
  m_boxDomain = BoxType {boxMin, boxMax};
#else
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    m_boxDomain.addPoint(surfPts[inode]);
  }
#endif

  // Initialize BVH with the surface elements.

  BoxType* boxes = axom::allocate<BoxType>(ncells, allocatorID);
  for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(axom::IndexType icell) {
      boxes[icell] = getCellBoundingBox(icell, surfaceData, surfPts);
    });

  // Build bounding volume hierarchy
  m_bvh.setAllocatorID(allocatorID);
  int result = m_bvh.initialize(boxes, ncells);

  axom::deallocate(boxes);
  return (result == spin::BVH_BUILD_OK);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline double SignedDistance<NDIMS, ExecSpace>::computeDistance(
  const PointType& pt) const
{
  PointType closest_pt;
  double dist;
  this->computeDistances(1, &pt, &dist, &closest_pt);
  return (dist);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
template <typename PointIndexable>
inline void SignedDistance<NDIMS, ExecSpace>::computeDistances(
  int npts,
  PointIndexable queryPts,
  double* outSgnDist,
  PointType* outClosestPts) const
{
  SLIC_ASSERT(npts > 0);
  SLIC_ASSERT(m_surfaceMesh != nullptr);
  SLIC_ASSERT(outSgnDist != nullptr);

  // Get a device-useable iterator
  auto it = m_bvh.getTraverser();

  // Get mesh data
  const double* xs = m_surfaceMesh->getCoordinateArray(0);
  const double* ys = m_surfaceMesh->getCoordinateArray(1);
  const double* zs = nullptr;
  if(NDIMS == 3)
  {
    zs = m_surfaceMesh->getCoordinateArray(2);
  }

  ZipPoint surf_pts {{xs, ys, zs}};

  const bool watertightInput = m_isInputWatertight;
  const BoxType boxDomain = m_boxDomain;
  const bool computeSigns = m_computeSign;

  detail::UcdMeshData surfaceData;
  bool result = detail::SD_GetUcdMeshData(m_surfaceMesh, surfaceData);
  AXOM_UNUSED_VAR(result);
  SLIC_CHECK_MSG(result, "Input mesh is not an unstructured surface mesh");

  AXOM_PERF_MARK_SECTION(
    "ComputeDistances",
    for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(int32 idx) {
        PointType qpt = queryPts[idx];

        MinCandidate curr_min {};

        auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
          int candidate_idx = leaf_nodes[current_node];

          checkCandidate(qpt,
                         curr_min,
                         candidate_idx,
                         surfaceData,
                         surf_pts,
                         computeSigns);
        };

        auto traversePredicate = [&](const PointType& p,
                                     const BoxType& bb) -> bool {
          return axom::primal::squared_distance(p, bb) <= curr_min.minSqDist;
        };

        // Traverse the tree, searching for the point with minimum distance.
        it.traverse_tree(qpt, searchMinDist, traversePredicate);

        double sgn = 1.0;
        if(computeSigns)
        {
          // STEP 0: if point is outside the bounding box of the surface mesh, then
          // it is outside, just return 1.0
          if(!(watertightInput && !boxDomain.contains(curr_min.minPt)))
          {
            sgn = computeSign(qpt, curr_min);
          }
        }

        outSgnDist[idx] = sqrt(curr_min.minSqDist) * sgn;
        if(outClosestPts)
        {
          outClosestPts[idx] = curr_min.minPt;
        }
      }););
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline axom::primal::BoundingBox<double, NDIMS>
SignedDistance<NDIMS, ExecSpace>::getCellBoundingBox(axom::IndexType icell,
                                                     const detail::UcdMeshData& mesh,
                                                     ZipPoint meshPts)
{
  // Get the cell type, for now we support linear triangle,quad in 3-D and
  // line segments in 2-D.
  const mint::CellType cellType = mesh.getCellType(icell);
  SLIC_ASSERT(cellType == mint::TRIANGLE || cellType == mint::QUAD ||
              cellType == mint::SEGMENT);

  // Get the cell node IDs that make up the cell
  int nnodes;
  const axom::IndexType* cellIds = mesh.getCellNodeIDs(icell, nnodes);

  // compute the cell's bounding box
  BoxType bb;

  for(int i = 0; i < nnodes; ++i)
  {
    bb.addPoint(meshPts[cellIds[i]]);
  }  // END for all cell nodes

  return (bb);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline void SignedDistance<NDIMS, ExecSpace>::checkCandidate(
  const PointType& qpt,
  MinCandidate& currMin,
  IndexType cellId,
  const detail::UcdMeshData& mesh,
  ZipPoint meshPts,
  bool computeNormal)
{
  int nnodes;
  const IndexType* nodes = mesh.getCellNodeIDs(cellId, nnodes);
  SLIC_ASSERT(nnodes <= 4);

  TriangleType surface_elems[2];

  int num_candidates = 1;
  surface_elems[0] =
    TriangleType {meshPts[nodes[0]], meshPts[nodes[1]], meshPts[nodes[2]]};

  if(nnodes == 4)
  {
    num_candidates = 2;
    surface_elems[1] =
      TriangleType {meshPts[nodes[0]], meshPts[nodes[2]], meshPts[nodes[3]]};
  }

  using axom::primal::closest_point;
  using axom::primal::squared_distance;
  using axom::utilities::isNearlyEqual;
  using detail::getClosestPointLocType;
  using detail::isClosestPointTypeShared;
  constexpr double EPS = 1e-8;

  for(int ei = 0; ei < num_candidates; ei++)
  {
    int candidate_loc;
    PointType candidate_pt =
      closest_point(qpt, surface_elems[ei], &candidate_loc, EPS);
    double sq_dist = squared_distance(qpt, candidate_pt);

    // Check the type of intersection we found
    const auto cpt_type = getClosestPointLocType(candidate_loc);

    // Determine if the closest point is on an edge or vertex
    const bool is_cpt_shared = isClosestPointTypeShared(cpt_type);

    bool shouldUpdateNormals = false;

    if(sq_dist < currMin.minSqDist)
    {
      // Clear the sum of normals if:
      const bool shouldClearVecs =
        // we're not computing normals
        !computeNormal ||
        // or we're not in a shared configuration
        !is_cpt_shared ||
        // or, if previous closest point type was different than current
        (currMin.minType != cpt_type) ||
        // finally, if there was a previous shared point -- check if approximately same as current
        !isNearlyEqual(squared_distance(candidate_pt, currMin.minPt), 0., EPS);

      currMin.minSqDist = sq_dist;
      currMin.minPt = candidate_pt;
      currMin.minType = cpt_type;
      currMin.minElem = cellId;
      currMin.minTri = surface_elems[ei];

      if(shouldClearVecs)
      {
        currMin.sumNormals = VectorType {};
        currMin.sumNormalsAngWt = VectorType {};
      }

      shouldUpdateNormals = (computeNormal && is_cpt_shared);
    }
    else
    {
      shouldUpdateNormals = computeNormal && is_cpt_shared &&
        (currMin.minType == cpt_type) &&
        isNearlyEqual(squared_distance(candidate_pt, currMin.minPt), 0., EPS);
    }

    if(shouldUpdateNormals)
    {
      VectorType norm = surface_elems[ei].normal();

      switch(cpt_type)
      {
      case detail::ClosestPointLocType::edge:
        // Candidate closest point is on an edge - add the normal of a
        // potentially-adjacent face
        currMin.sumNormals += norm;
        break;
      case detail::ClosestPointLocType::face:
        if(!surface_elems[ei].degenerate())
        {
          // Candidate closest point is on a vertex - add the angle-weighted
          // normal of a face potentially sharing a vertex
          double alpha = surface_elems[ei].angle(candidate_loc);
          currMin.sumNormalsAngWt += (norm.unitVector() * alpha);
        }
        else
        {
          SLIC_WARNING("Triangle for closest point was degenerate");
        }
        break;
      default:
        break;  // no-op
      }
    }
  }
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline double SignedDistance<NDIMS, ExecSpace>::computeSign(
  const PointType& qpt,
  const MinCandidate& currMin)
{
  double sgn = 1.0;
  // STEP 1: Select the pseudo-normal N at the closest point to calculate the sign.
  // There are effectively 3 cases based on the location of the closest point.
  VectorType N;
  switch(currMin.minType)
  {
  case detail::ClosestPointLocType::face:
    // CASE 1: closest point is on the face of the surface element
    N = currMin.minTri.normal();
    break;
  case detail::ClosestPointLocType::edge:
    // CASE 2: closest point is on an edge, use sum of normals of equidistant faces
    // TODO: Sometimes, the traversal fails to find the opposite face, so only
    // a single face's normal is accumulated here. The proper solution would be
    // to precompute edge pseudo-normals during construction time, but that
    // would also require generating cell-to-face connectivity for the surface mesh.
    N = currMin.sumNormals;
    break;
  case detail::ClosestPointLocType::vertex:
    // CASE 3: closest point is on a node, use angle-weighted pseudo-normal
    N = currMin.sumNormalsAngWt;
    break;
  default:
    SLIC_WARNING("Did not find a valid closest point for query point: " << qpt);
    break;
  }

  // STEP 2: Given the pseudo-normal, N, and the vector r from the closest point
  // to the query point, compute the sign by checking the sign of their dot product.
  VectorType r(currMin.minPt, qpt);
  double dotprod = r.dot(N);
  sgn = (dotprod >= 0.0) ? 1.0 : -1.0;
  return sgn;
}

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SIGNED_DISTANCE_HPP_
