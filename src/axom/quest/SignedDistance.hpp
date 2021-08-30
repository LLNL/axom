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
template <int NDIMS>
class SignedDistance
{
public:
  using PointType = axom::primal::Point<double, NDIMS>;
  using VectorType = axom::primal::Vector<double, NDIMS>;
  using TriangleType = axom::primal::Triangle<double, NDIMS>;
  using BoxType = axom::primal::BoundingBox<double, NDIMS>;
  using BVHTreeType = axom::spin::BVH<NDIMS>;

  using ExecSpace = axom::SEQ_EXEC;

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
   * \brief Destructor.
   */
  ~SignedDistance();

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

  void computeDistances(int npts,
                        const PointType* queryPts,
                        double* out_sdist,
                        PointType* out_closestPts = nullptr) const;

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
template <int NDIMS>
SignedDistance<NDIMS>::SignedDistance(const mint::Mesh* surfaceMesh,
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
  SLIC_ASSERT(result == BVH_BUILD_OK);
}

//------------------------------------------------------------------------------
template <int NDIMS>
SignedDistance<NDIMS>::~SignedDistance()
{ }

//------------------------------------------------------------------------------
template <int NDIMS>
inline double SignedDistance<NDIMS>::computeDistance(const PointType& pt) const
{
  PointType closest_pt;
  double dist;
  this->computeDistances(1, &pt, &dist, &closest_pt);
  return (dist);
}

//------------------------------------------------------------------------------
template <int NDIMS>
inline void SignedDistance<NDIMS>::computeDistances(int npts,
                                                    const PointType* queryPts,
                                                    double* out_sdist,
                                                    PointType* out_closestPts) const
{
  SLIC_ASSERT(m_surfaceMesh != nullptr);
  SLIC_ASSERT(queryPts != nullptr);
  SLIC_ASSERT(out_sdist != nullptr);

  using ZipPoint = primal::ZipIndexable<PointType>;

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

  AXOM_PERF_MARK_SECTION(
    "ComputeDistances",
    for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(int32 idx) {
        PointType qpt = queryPts[idx];

        double minSqDist = numerics::floating_point_limits<double>::max();
        PointType minPt;
        int minLoc;
        int minElem;
        TriangleType minTri;

        auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
          int candidate_idx = leaf_nodes[current_node];

          IndexType nodes[4];
          IndexType nnodes = m_surfaceMesh->getCellNodeIDs(candidate_idx, nodes);

          TriangleType surface_elems[2];

          int num_candidates = 1;
          surface_elems[0] = TriangleType {surf_pts[nodes[0]],
                                           surf_pts[nodes[1]],
                                           surf_pts[nodes[2]]};

          if(nnodes == 4)
          {
            num_candidates = 2;
            surface_elems[1] = TriangleType {surf_pts[nodes[0]],
                                             surf_pts[nodes[2]],
                                             surf_pts[nodes[3]]};
          }

          for(int ei = 0; ei < num_candidates; ei++)
          {
            int candidate_loc;
            PointType candidate_pt =
              axom::primal::closest_point(qpt, surface_elems[ei], &candidate_loc);
            double sq_dist = axom::primal::squared_distance(qpt, candidate_pt);
            if(sq_dist < minSqDist)
            {
              minSqDist = sq_dist;
              minPt = candidate_pt;
              minLoc = candidate_loc;
              minElem = candidate_idx;
              minTri = surface_elems[ei];
            }
          }
        };

        auto traversePredicate = [&](const PointType& p,
                                     const BoxType& bb) -> bool {
          return axom::primal::squared_distance(p, bb) <= minSqDist;
        };

        // Traverse the tree, searching for the point with minimum distance.
        it.traverse_tree(qpt, searchMinDist, traversePredicate);

        double sgn = 1.0;
        if(computeSigns)
        {
          //TODO

          // STEP 0: if point is outside the bounding box of the surface mesh, then
          // it is outside, just return 1.0
          if(!(watertightInput && !boxDomain.contains(minPt)))
          {
            // CASE 1: closest point is on the face of the surface element
            VectorType N = minTri.normal();
            VectorType r(minPt, qpt);
            double dotprod = r.dot(N);
            sgn = (dotprod >= 0.0) ? 1.0 : -1.0;
          }
        }

        out_sdist[idx] = sqrt(minSqDist) * sgn;
        if(out_closestPts)
        {
          out_closestPts[idx] = minPt;
        }
      }););
}

//------------------------------------------------------------------------------
template <int NDIMS>
inline axom::primal::BoundingBox<double, NDIMS>
SignedDistance<NDIMS>::getCellBoundingBox(axom::IndexType icell)
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

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SIGNED_DISTANCE_HPP_
