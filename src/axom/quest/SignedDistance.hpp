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
struct UcdMeshData
{
  int shape_type;
  const IndexType* cells_to_nodes;
  IndexType nodes_per_cell;
  const IndexType* cell_node_offsets;

  AXOM_HOST_DEVICE int getCellNodeIDs(IndexType cellId, IndexType* outNodes) const;
};

bool SD_GetUcdMeshData(const mint::Mesh* surfaceMesh, UcdMeshData& outSurfData);

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
    double minSqDist =
      numerics::floating_point_limits<double>::max();  // Squared distance to query point
    PointType minPt;      // Closest point on element
    int minLoc;           // Location of closest point on element
    int minElem;          // Closest element index in mesh
    TriangleType minTri;  // The actual cell element
  };

public:
  /*!
   * \brief Creates a SignedDistance instance for queries on the given mesh.
   * \param [in] surfaceMesh user-supplied surface mesh.
   * \param [in] isWatertight indicates if the surface mesh is closed.
   * \param [in] maxObjects max number of objects for spatial decomposition.
   * \param [in] maxLevels max levels for spatial decomposition.
   * \param [in] computeSign indicates if distance queries should compute signs (optional).
   * \note computeSign defaults to \a true when not specified.
   * \pre surfaceMesh != nullptr
   */
  SignedDistance(const mint::Mesh* surfaceMesh,
                 bool isWatertight,
                 int maxObjects,
                 int maxLevels,
                 bool computeSign = true);

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
   * \return box bounding box of the cell.
   * \pre m_surfaceMesh != nullptr
   * \pre icell >= 0 && icell < m_surfaceMesh->getNumberOfCells()
   */
  BoxType getCellBoundingBox(axom::IndexType icell);

  /*!
   * \brief Checks a given candidate surface element against a query point
   *  and updates the minimum-distance candidate data if the element is closer.
   *
   * \param [in] qpt query point to check against surface element
   * \param [in,out] currMin the minimum-distance candidate data to update
   * \param [in] cellId the surface element ID to evaluate
   * \param [in] mesh the surface mesh data
   * \param [in] meshPts the surface mesh point coordinate data
   */
  AXOM_HOST_DEVICE static void checkCandidate(const PointType& qpt,
                                              MinCandidate& currMin,
                                              IndexType cellId,
                                              const detail::UcdMeshData& mesh,
                                              ZipPoint meshPts);

  /*!
   * \brief Computes the sign of the given query point given the closest point
   *  data.
   *
   * \param [in] qpt query point to check against surface element
   * \param [in] currMin the minimum-distance surface element data
   * \param [in] mesh the surface mesh data
   * \param [in] meshPts the surface mesh point coordinate data
   *
   * \return sgn 1.0 if outside, -1.0 if inside
   */
  AXOM_HOST_DEVICE static double computeSign(const PointType& qpt,
                                             const MinCandidate& currMin,
                                             const detail::UcdMeshData& mesh,
                                             ZipPoint meshPts);

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
                                                 int maxObjects,
                                                 int maxLevels,
                                                 bool computeSign)
  : m_isInputWatertight(isWatertight)
  , m_computeSign(computeSign)
{
  // Sanity checks
  SLIC_ASSERT(surfaceMesh != nullptr);
  SLIC_ASSERT(maxLevels >= 1);

  m_surfaceMesh = surfaceMesh;
  const axom::IndexType ncells = m_surfaceMesh->getNumberOfCells();
  const axom::IndexType nnodes = m_surfaceMesh->getNumberOfNodes();

  // compute bounding box of surface mesh
  // NOTE: this should be changed to an oriented bounding box in the future.
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    PointType pt;
    for(int i = 0; i < NDIMS; ++i)
    {
      pt[i] = m_surfaceMesh->getCoordinateArray(i)[inode];
    }
    m_boxDomain.addPoint(pt);
  }

  // Initialize BVH with the surface elements.

  BoxType* boxes = axom::allocate<BoxType>(ncells);
  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    boxes[icell] = this->getCellBoundingBox(icell);
  }  // END for all cells

  // Build bounding volume hierarchy
  int result = m_bvh.initialize(boxes, ncells);
  SLIC_ASSERT(result == spin::BVH_BUILD_OK);
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
  auto it = m_bvh.getIterator();

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

          checkCandidate(qpt, curr_min, candidate_idx, surfaceData, surf_pts);
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
            sgn = computeSign(qpt, curr_min, surfaceData, surf_pts);
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
SignedDistance<NDIMS, ExecSpace>::getCellBoundingBox(axom::IndexType icell)
{
  // Sanity checks
  SLIC_ASSERT(m_surfaceMesh != nullptr);
  SLIC_ASSERT(icell >= 0 && icell < m_surfaceMesh->getNumberOfCells());

  // Get the cell type, for now we support linear triangle,quad in 3-D and
  // line segments in 2-D.
  const mint::CellType cellType = m_surfaceMesh->getCellType(icell);
  SLIC_ASSERT(cellType == mint::TRIANGLE || cellType == mint::QUAD ||
              cellType == mint::SEGMENT);
  const int nnodes = axom::mint::getCellInfo(cellType).num_nodes;

  // Get the cell node IDs that make up the cell
  axom::IndexType* cellIds = new axom::IndexType[nnodes];
  m_surfaceMesh->getCellNodeIDs(icell, cellIds);

  // compute the cell's bounding box
  BoxType bb;
  PointType pt;

  for(int i = 0; i < nnodes; ++i)
  {
    m_surfaceMesh->getNode(cellIds[i], pt.data());
    bb.addPoint(pt);
  }  // END for all cell nodes

  // clean up all dynamically allocated memory
  delete[] cellIds;

  return (bb);
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline void SignedDistance<NDIMS, ExecSpace>::checkCandidate(
  const PointType& qpt,
  MinCandidate& currMin,
  IndexType cellId,
  const detail::UcdMeshData& mesh,
  ZipPoint meshPts)
{
  IndexType nodes[4];
  IndexType nnodes = mesh.getCellNodeIDs(cellId, nodes);
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

  for(int ei = 0; ei < num_candidates; ei++)
  {
    int candidate_loc;
    PointType candidate_pt =
      axom::primal::closest_point(qpt, surface_elems[ei], &candidate_loc);
    double sq_dist = axom::primal::squared_distance(qpt, candidate_pt);
    if(sq_dist < currMin.minSqDist)
    {
      currMin.minSqDist = sq_dist;
      currMin.minPt = candidate_pt;
      currMin.minLoc = candidate_loc;
      currMin.minElem = cellId;
      currMin.minTri = surface_elems[ei];
    }
  }
}

//------------------------------------------------------------------------------
template <int NDIMS, typename ExecSpace>
inline double SignedDistance<NDIMS, ExecSpace>::computeSign(
  const PointType& qpt,
  const MinCandidate& currMin,
  const detail::UcdMeshData& mesh,
  ZipPoint meshPts)
{
  // TODO

  double sgn = 1.0;
  // CASE 1: closest point is on the face of the surface element
  VectorType N = currMin.minTri.normal();
  VectorType r(currMin.minPt, qpt);
  double dotprod = r.dot(N);
  sgn = (dotprod >= 0.0) ? 1.0 : -1.0;
  return sgn;
}

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SIGNED_DISTANCE_HPP_
