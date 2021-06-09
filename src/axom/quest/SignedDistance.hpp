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
#include "axom/spin/BVHTree.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/Vector.hpp"

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
  using BVHTreeType = axom::spin::BVHTree<int, NDIMS>;

private:
  /// @{
  /// \name Internal Datatype Definitions

  struct cpt_data
  {
    PointType closest_point;
    int candidate_index;
    int cpt_location;

    int nelems;
    std::vector<TriangleType> surface_elements;
    std::vector<axom::IndexType> element_ids;
    std::vector<PointType> closest_pts;
    std::vector<int> cpt_locs;
  };

  /// @}

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

  /*!
   * \brief Computes the distance of the given point to the surface mesh.
   * \note This is an overloaded method that also returns the BVH buckets and
   *  corresponding triangles used to calculate the signed distance.
   *
   * \param [in]  queryPnt user-supplied point.
   * \param [out] bvh_buckets the buckets of the BVH used to satisfy the query.
   * \param [out] triangles the triangles of the BVH used to satisfy the query.
   * \param [out] my_triangles the triangle used to compute the pseudo-normal.
   *
   * \note The variables 'triangles'/'my_triangles' are relevant in debug mode.
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
   * \return minDist the minimum signed distance to the surface mesh.
   */
  double computeDistance(const PointType& queryPnt,
                         std::vector<int>& bvh_buckets,
                         std::vector<axom::IndexType>& triangles,
                         std::vector<axom::IndexType>& my_triangles,
                         PointType& closest_pt) const;

  /*!
   * \brief Returns a const reference to the underlying bucket tree.
   * \return ptr pointer to the underlying bucket tree
   * \post ptr != nullptr
   */
  const BVHTreeType* getBVHTree() const { return m_bvhTree; }

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
   * \brief Searches through the list of candidates surface elements and
   *  computes the minimum squared distance and closest point to a surface
   *  element.
   *
   * \param [in]  pt the query point.
   * \param [in]  candidates list of IDs to the candidate surface elements
   * \param [in]  nelems total number of candidates
   * \param [out] cpt_data closest point data
   *
   * \note The list of candidates is provided by getCandidateSurfaceElements
   *
   * \see SignedDistance::getCandidateSurfaceElements()
   * \see SignedDistance::computeDistance()
   *
   * \pre candidates != nullptr
   * \pre surface_elements != nullptr
   * \pre elementIds != nullptr
   * \pre closest_pts != nullptr
   * \pre clocs != nullptr
   *
   * \post index >= 0 && index < nelems
   *
   * \return dist minimum squared distance to the surface.
   */
  double getMinSqDistance(const PointType& pt,
                          const int* candidates,
                          int nelems,
                          cpt_data* cpt) const;

  /*!
   * \brief Computes the sign of the given query point given the closest point
   *  and surface element on the surface mesh and neighboring surface elements.
   *
   * \param [in] pt the query point
   * \param [in] cpt pointer to the closest point data for this point query.
   * \param [out] my_elements list of element used to calculate pseudo-normal
   *
   * \return sign the calculated sign, 1.0 if outside, -1.0 if inside
   *
   * \pre elementIds != nullptr
   * \pre surface_elements != nullptr
   * \pre closest_pts != nullptr
   * \pre clocs != nullptr
   */
  double computeSign(const PointType& pt,
                     const cpt_data* cpt,
                     std::vector<axom::IndexType>& my_elements) const;

  /*!
   * \brief Returns a sorted list of the candidate surface elements.
   *
   * \param [in] pt the query point.
   * \param [in] bins array of the bins to check.
   * \param [in] nbins size of the bins array.
   * \param [out] surface_elements buffer to store candidate surface elements
   * \param [out] indx indirection array to access surface elements
   * \param [in] nelems total number of elements
   *
   * \note The candidate surface elements are sorted in ascending order, based
   *  on the distance of the query point and the bounding box of the surface
   *  element.
   *
   *  \pre bins != nullptr
   *  \pre m_surfaceMesh != nullptr
   *  \pre surface_elements != nullptr
   *  \pre indx != nullptr
   */
  void getCandidateSurfaceElements(const PointType& pt,
                                   const int* bins,
                                   int nbins,
                                   std::vector<int>& candidates) const;

  /*!
   * \brief Returns the maximum distance between a given point and box.
   * \param [in] b user-supplied axis-aligned bounding box.
   * \param [in] pt user-supplied point.
   * \return d maximum distance from a point to a box.
   */
  double getMaxSqDistance(const BoxType& b, const PointType& pt) const;

  /*!
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use from the calling application.
   */
  SignedDistance() : m_surfaceMesh(nullptr), m_bvhTree(nullptr) { }

private:
  bool m_isInputWatertight;        /*!< indicates if input is watertight     */
  bool m_computeSign;              /*!< indicates if queries compute sign    */
  const mint::Mesh* m_surfaceMesh; /*!< User-supplied surface mesh.          */
  BoxType m_boxDomain;             /*!< bounding box containing surface mesh */
  BVHTreeType* m_bvhTree;          /*!< Spatial acceleration data-structure. */

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
namespace detail
{
class SortByDistance
{
public:
  SortByDistance(double* dist) : m_dist(dist) { }
  ~SortByDistance() { }
  bool operator()(int i, int j) const { return (m_dist[i] < m_dist[j]); }

private:
  double* m_dist;
};

}  // end namespace detail

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

  // Initialize BucketTree with the surface elements.
  m_bvhTree = new BVHTreeType(ncells, maxLevels);

  for(axom::IndexType icell = 0; icell < ncells; ++icell)
  {
    m_bvhTree->insert(this->getCellBoundingBox(icell), icell);
  }  // END for all cells

  // Build bounding volume hierarchy
  m_bvhTree->build(maxObjects);
}

//------------------------------------------------------------------------------
template <int NDIMS>
SignedDistance<NDIMS>::~SignedDistance()
{
  delete m_bvhTree;
  m_bvhTree = nullptr;
}

//------------------------------------------------------------------------------
template <int NDIMS>
inline double SignedDistance<NDIMS>::computeDistance(const PointType& pt) const
{
  std::vector<int> buckets;
  std::vector<axom::IndexType> elements;
  std::vector<axom::IndexType> my_elements;
  PointType closest_pt;
  double dist =
    this->computeDistance(pt, buckets, elements, my_elements, closest_pt);
  return (dist);
}

//------------------------------------------------------------------------------
template <int NDIMS>
inline double SignedDistance<NDIMS>::computeDistance(
  const PointType& pt,
  std::vector<int>& buckets,
  std::vector<axom::IndexType>& AXOM_DEBUG_PARAM(elementIds),
  std::vector<axom::IndexType>& my_elements,
  PointType& closest_pt) const
{
  SLIC_ASSERT(m_surfaceMesh != nullptr);
  SLIC_ASSERT(m_bvhTree != nullptr);

  // STEP 0: get list of buckets to satisfy point query
  m_bvhTree->find(pt, buckets);

  // STEP 1: get candidate surface elements
  int nbuckets = static_cast<int>(buckets.size());
  std::vector<int> candidates;
  this->getCandidateSurfaceElements(pt, &buckets[0], nbuckets, candidates);

  const int nelems = static_cast<int>(candidates.size());

  // STEP 2: process surface elements and compute minimum distance and
  // corresponding closest point.
  cpt_data cpt;
  double minSqDist = this->getMinSqDistance(pt, &candidates[0], nelems, &cpt);
  closest_pt = cpt.closest_point;
#ifdef AXOM_DEBUG
  elementIds = cpt.element_ids;
#endif

  // STEP 3: compute sign
  double sign = m_computeSign ? this->computeSign(pt, &cpt, my_elements) : 1.;

  // STEP 4: return computed signed distance
  return (sign * std::sqrt(minSqDist));
}

//------------------------------------------------------------------------------
template <int NDIMS>
double SignedDistance<NDIMS>::computeSign(
  const PointType& pt,
  const cpt_data* cpt,
  std::vector<axom::IndexType>& AXOM_DEBUG_PARAM(my_elements)) const
{
  // Sanity checks
  SLIC_ASSERT(cpt != nullptr);

  // STEP 0: if point is outside the bounding box of the surface mesh, then
  // it is outside, just return 1.0
  if(m_isInputWatertight && !m_boxDomain.contains(pt))
  {
    /* short-circuit, pt outside bounding box of the surface mesh */
    return 1.0;
  }

  // STEP 1: Otherwise, calculate pseudo-normal N, at the closest point to
  // calculate the sign. There are effectively 3 cases based on the location
  // of the closest point

  VectorType N;  // the pseudo-normal, computed below

  const int cpt_loc = cpt->cpt_location;
  const int index = cpt->candidate_index;
  const int nelems = cpt->nelems;

  if(cpt_loc >= TriangleType::NUM_TRI_VERTS)
  {
    // CASE 1: closest point is on the face of the surface element
    N = cpt->surface_elements[index].normal();
#ifdef AXOM_DEBUG
    my_elements.push_back(cpt->element_ids[index]);
#endif
  }
  else if(cpt_loc < 0)
  {
    // CASE 2: closest point is on an edge, sum normals of adjacent facets
    for(int i = 0; i < nelems; ++i)
    {
      double dist =
        axom::primal::squared_distance(cpt->closest_point, cpt->closest_pts[i]);

      if(axom::utilities::isNearlyEqual(dist, 0.0))
      {
        N += cpt->surface_elements[i].normal();
#ifdef AXOM_DEBUG
        my_elements.push_back(cpt->element_ids[i]);
#endif
      }  // END if

    }  // END for
  }
  else
  {
    // CASE 3: closest point is on a node, use angle weighted pseudo-normal
    for(int i = 0; i < nelems; ++i)
    {
      double dist =
        axom::primal::squared_distance(cpt->closest_point, cpt->closest_pts[i]);

      if(axom::utilities::isNearlyEqual(dist, 0.0))
      {
        double alpha = cpt->surface_elements[i].angle(cpt->cpt_locs[i]);
        N += (cpt->surface_elements[i].normal().unitVector() * alpha);
#ifdef AXOM_DEBUG
        my_elements.push_back(cpt->element_ids[i]);
#endif
      }  // END if

    }  // END for
  }

  // STEP 2: Given the pseudo-normal, N, and the vector r from the closest point
  // to the query point, compute the sign by checking the sign of their dot
  // product.
  VectorType r(cpt->closest_point, pt);
  double dotprod = r.dot(N);
  double sign = (dotprod >= 0.0) ? 1.0 : -1.0;
  SLIC_ASSERT(sign == -1.0 || sign == 1.0);

  return sign;
}

//------------------------------------------------------------------------------
template <int NDIMS>
double SignedDistance<NDIMS>::getMinSqDistance(const PointType& pt,
                                               const int* candidates,
                                               int nelems,
                                               cpt_data* cpt) const
{
  SLIC_ASSERT(candidates != nullptr);
  SLIC_ASSERT(cpt != nullptr);

  cpt->nelems = nelems;
  cpt->surface_elements.resize(nelems);
  cpt->closest_pts.resize(nelems);
  cpt->cpt_locs.resize(nelems);
  cpt->element_ids.resize(nelems);

  PointType* closest_pts = &(cpt->closest_pts)[0];
  int* cpt_locs = &(cpt->cpt_locs)[0];

  double minSqDist = std::numeric_limits<double>::max();

  for(int i = 0; i < nelems; ++i)
  {
    const int cellIdx = candidates[i];
    cpt->element_ids[i] = cellIdx;

    axom::IndexType cellIds[3];
    m_surfaceMesh->getCellNodeIDs(cellIdx, cellIds);

    TriangleType surface_element;
    m_surfaceMesh->getNode(cellIds[0], surface_element[0].data());
    m_surfaceMesh->getNode(cellIds[1], surface_element[1].data());
    m_surfaceMesh->getNode(cellIds[2], surface_element[2].data());
    cpt->surface_elements[i] = surface_element;

    closest_pts[i] =
      axom::primal::closest_point(pt, surface_element, &cpt_locs[i]);

    double sqDist = axom::primal::squared_distance(pt, closest_pts[i]);
    if(sqDist < minSqDist)
    {
      minSqDist = sqDist;
      cpt->closest_point = closest_pts[i];
      cpt->cpt_location = cpt_locs[i];
      cpt->candidate_index = i;
    }

  }  // END for all elements

  return minSqDist;
}

//------------------------------------------------------------------------------
template <int NDIMS>
double SignedDistance<NDIMS>::getMaxSqDistance(const BoxType& b,
                                               const PointType& pt) const
{
  std::vector<PointType> pnts;
  BoxType::getPoints(b, pnts);

  const int npoints = pnts.size();
  SLIC_ASSERT(npoints == 4 || npoints == 8);

  double dist = std::numeric_limits<double>::min();
  for(int i = 0; i < npoints; ++i)
  {
    dist = std::max(dist, axom::primal::squared_distance(pt, pnts[i]));
  }  // END for all points

  return dist;
}

//------------------------------------------------------------------------------
template <int NDIMS>
void SignedDistance<NDIMS>::getCandidateSurfaceElements(
  const PointType& pt,
  const int* bins,
  int nbins,
  std::vector<int>& candidates) const
{
  // Sanity checks
  SLIC_ASSERT(m_surfaceMesh != nullptr);
  SLIC_ASSERT(bins != nullptr);

  // STEP 0: count total number of surface elements
  int nelems = 0;
  for(int ibin = 0; ibin < nbins; ++ibin)
  {
    const int bucketIdx = bins[ibin];
    nelems += m_bvhTree->getBucketNumObjects(bucketIdx);
  }

  // STEP 0: allocate array to store the distance of the bounding box of each
  // surface element to the query point.
  int* surface_elements = new int[nelems];
  double* dist = new double[nelems];
  int* objectIds = new int[nelems];
  int* indx = new int[nelems];

  // STEP 1: get flat array of the surface elements and compute distances
  int icount(static_cast<int>(0));
  for(int ibin = 0; ibin < nbins; ++ibin)
  {
    const int binIdx = bins[ibin];
    const int nobjects = m_bvhTree->getBucketNumObjects(binIdx);
    const int* objList = m_bvhTree->getBucketObjectArray(binIdx);

    for(int iobject = 0; iobject < nobjects; ++iobject)
    {
      const int objectId = objList[iobject];
      const int cellIdx = m_bvhTree->getObjectData(objectId);
      BoxType bbox = m_bvhTree->getObjectBox(objectId);

      //        dist[ icount ]           = this->getMaxSqDistance( bbox, pt );
      dist[icount] = axom::primal::squared_distance(pt, bbox);
      objectIds[icount] = objectId;
      surface_elements[icount] = cellIdx;
      indx[icount] = icount;
      ++icount;

    }  // END for all objects within the bin

  }  // END for all bins

  SLIC_ASSERT(icount == nelems);

  // STEP 2: sort the indices based on distance
  std::sort(indx, indx + nelems, detail::SortByDistance(dist));

  // STEP 3: filter out candidates
  candidates.push_back(surface_elements[indx[0]]);
  const BoxType& bbox = m_bvhTree->getObjectBox(objectIds[indx[0]]);
  double maxDist = this->getMaxSqDistance(bbox, pt);

  for(int i = 1; i < nelems; ++i)
  {
    const int idx = indx[i];
    if(dist[idx] < maxDist)
    {
      candidates.push_back(surface_elements[idx]);
    }
  }

  // STEP 4: delete temporary distance array
  delete[] dist;
  delete[] objectIds;
  delete[] indx;
  delete[] surface_elements;
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
