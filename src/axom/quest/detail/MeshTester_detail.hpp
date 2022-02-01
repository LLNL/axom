// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESH_TESTER_DETAIL_HPP_
#define AXOM_QUEST_MESH_TESTER_DETAIL_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/mint.hpp"

// RAJA includes
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
  #include "axom/mint/execution/internal/structured_exec.hpp"
#endif

// Acceleration data structure includes
#include "axom/spin/BVH.hpp"
#include "axom/spin/ImplicitGrid.hpp"

namespace axom
{
namespace quest
{
namespace detail
{
using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
using Triangle3 = primal::Triangle<double, 3>;
using SpatialBoundingBox = primal::BoundingBox<double, 3>;
using UniformGrid3 = spin::UniformGrid<int, 3>;
using Point3 = primal::Point<double, 3>;

inline detail::Triangle3 getMeshTriangle(axom::IndexType i,
                                         detail::UMesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getNumberOfCellNodes(i) == 3);

  detail::Triangle3 tri;

  const axom::IndexType* triCell = surface_mesh->getCellNodeIDs(i);

  const double* x = surface_mesh->getCoordinateArray(mint::X_COORDINATE);
  const double* y = surface_mesh->getCoordinateArray(mint::Y_COORDINATE);
  const double* z = surface_mesh->getCoordinateArray(mint::Z_COORDINATE);

  for(int n = 0; n < 3; ++n)
  {
    const axom::IndexType nodeIdx = triCell[n];
    tri[n][0] = x[nodeIdx];
    tri[n][1] = y[nodeIdx];
    tri[n][2] = z[nodeIdx];
  }

  return tri;
}
enum class AccelType
{
  ImplicitGrid,
  BVH
};

template <AccelType accel, typename ExecSpace, typename FloatType>
struct CandidateFinder;

template <typename ExecSpace, typename FloatType>
struct CandidateFinderBase;

/*!
 * \class CandidateFinderBase
 *
 * \brief Base class to handle common operations for mesh testing.
 */
template <typename ExecSpace, typename FloatType>
struct CandidateFinderBase
{
  using BoxType = typename primal::BoundingBox<FloatType, 3>;
  using PointType = typename primal::Point<FloatType, 3>;
#ifdef AXOM_USE_UMPIRE
  static constexpr bool ExecOnDevice =
    axom::execution_space<ExecSpace>::onDevice();
  static constexpr MemorySpace Space =
    ExecOnDevice ? axom::MemorySpace::Device : axom::MemorySpace::Host;
  static constexpr MemorySpace HostSpace = axom::MemorySpace::Host;
#else
  static constexpr MemorySpace Space = axom::MemorySpace::Dynamic;
  static constexpr MemorySpace HostSpace = axom::MemorySpace::Dynamic;
#endif

  /*!
   * \brief Initializes the CandidateFinder query object with a surface mesh.
   *
   * \param [in] surface_mesh The triangular surface mesh to query.
   * \param [in] intersectionThreshold The tolerance threshold to use for
   *  triangle intersection tests.
   */
  CandidateFinderBase(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* surface_mesh,
                      double intersectionThreshold)
    : m_surfaceMesh(surface_mesh)
    , m_intersectionThreshold(intersectionThreshold)
  { }

  void initialize();

  /*!
   * \brief Runs a query to find pairs of triangle cell indices intersecting
   *  in the surface mesh, as well as degenerate triangles in the mesh.
   *
   * \param [in] firstIndex The first indices of intersecting pairs
   * \param [in] secondIndex The second indices of intersecting pairs
   * \param [in] degenerateIndices Indices of degenerate mesh triangles
   */
  void findTriMeshIntersections(axom::Array<IndexType>& firstIndex,
                                axom::Array<IndexType>& secondIndex,
                                axom::Array<IndexType>& degenerateIndices);

protected:
  /*!
   * \brief Returns the candidate intersection pairs by using a spatial query
   *  data structure.
   *
   * \param [out] offsets Offsets into the candidates array for each mesh cell
   * \param [out] counts  The number of candidates for each mesh cell
   * \return candidates The flat array of candidate indices
   *
   * \note Upon completion, the triangular mesh element at index i has:
   *  * counts[ i ] candidates
   *  * Candidate intersections with indices stored in the range:
   *    [offsets[i], offsets[i] + counts[i])
   * \note This should be implemented in derived CandidateFinder template
   *  specializations for each supported acceleration data structure.
   */
  virtual axom::ArrayView<IndexType, 1, Space> getCandidates(
    axom::Array<IndexType, 1, Space>& offsets,
    axom::Array<IndexType, 1, Space>& counts) = 0;

  mint::UnstructuredMesh<mint::SINGLE_SHAPE>* m_surfaceMesh;
  double m_intersectionThreshold;
  int m_ncells;
  axom::Array<detail::Triangle3, 1, Space> m_tris;
  axom::Array<BoxType, 1, Space> m_aabbs;
  axom::Array<int, 1, Space> m_degenerate;
};

template <typename ExecSpace, typename FloatType>
void CandidateFinderBase<ExecSpace, FloatType>::initialize()
{
  const int ncells = m_surfaceMesh->getNumberOfCells();

  m_tris.resize(ncells);
  axom::ArrayView<detail::Triangle3, 1, Space> p_tris = m_tris;

  // Marks each cell/triangle as degenerate (1) or not (0)
  m_degenerate.resize(ncells);
  axom::ArrayView<int, 1, Space> p_degenerate = m_degenerate;

  // Each access-aligned bounding box represented by 2 (x,y,z) points
  m_aabbs.resize(ncells);
  axom::ArrayView<BoxType, 1, Space> p_aabbs = m_aabbs;

  // Initialize the bounding box for each Triangle and marks
  // if the Triangle is degenerate.
  mint::for_all_cells<ExecSpace, mint::xargs::coords>(
    m_surfaceMesh,
    AXOM_LAMBDA(IndexType cellIdx,
                numerics::Matrix<double> & coords,
                const IndexType* nodeIds) {
      AXOM_UNUSED_VAR(nodeIds);

      detail::Triangle3 tri;

      for(IndexType inode = 0; inode < 3; ++inode)
      {
        const double* node = coords.getColumn(inode);
        tri[inode][0] = node[mint::X_COORDINATE];
        tri[inode][1] = node[mint::Y_COORDINATE];
        tri[inode][2] = node[mint::Z_COORDINATE];
      }  // END for all cells nodes

      p_degenerate[cellIdx] = (tri.degenerate() ? 1 : 0);

      p_tris[cellIdx] = tri;

      p_aabbs[cellIdx] = compute_bounding_box(tri);
    });
}

template <typename ExecSpace, typename FloatType>
void CandidateFinderBase<ExecSpace, FloatType>::findTriMeshIntersections(
  axom::Array<IndexType>& firstIndex,
  axom::Array<IndexType>& secondIndex,
  axom::Array<IndexType>& degenerateIndices)
{
  const int ncells = m_surfaceMesh->getNumberOfCells();

  using IndexArray = axom::Array<IndexType, 1, Space>;
  using IndexView = axom::ArrayView<IndexType, 1, Space>;

  using HostIndexArray = axom::Array<IndexType, 1, HostSpace>;
#ifdef AXOM_USE_RAJA
  using atomic_pol = typename axom::execution_space<ExecSpace>::atomic_policy;
#endif

  // Get CSR arrays for candidate data
  IndexArray offsets, counts;
  IndexView candidates = getCandidates(offsets, counts);

  IndexArray indices(candidates.size());
  IndexArray validCandidates(candidates.size());

  IndexView p_indices = indices;
  IndexView p_validCandidates = validCandidates;

  IndexType numCandidates;
  {
    IndexArray numValidCandidates(1);
    numValidCandidates.fill(0);
    IndexView p_numValidCandidates = numValidCandidates;

    IndexView p_offsets = offsets;
    IndexView p_counts = counts;

    // Initialize triangle indices and valid candidates
    for_all<ExecSpace>(
      ncells,
      AXOM_LAMBDA(IndexType i) {
        for(int j = 0; j < p_counts[i]; j++)
        {
          if(i < candidates[p_offsets[i] + j])
          {
#ifdef AXOM_USE_RAJA
            auto idx = RAJA::atomicAdd<atomic_pol>(&p_numValidCandidates[0], 1);
#else
            auto idx = p_numValidCandidates[0]++;
#endif
            p_indices[idx] = i;
            p_validCandidates[idx] = candidates[p_offsets[i] + j];
          }
        }
      });

    axom::copy(&numCandidates, numValidCandidates.data(), sizeof(IndexType));
  }

  IndexArray firstIsectPair(candidates.size());
  IndexArray secondIsectPair(candidates.size());
  IndexType isectCounter;
  {
    IndexArray numIsectPairs(1);
    IndexView p_numIsectPairs = numIsectPairs;

    IndexView p_firstIsectPair = firstIsectPair;
    IndexView p_secondIsectPair = secondIsectPair;
    axom::ArrayView<detail::Triangle3, 1, Space> p_tris = m_tris;

    double intersectionThreshold = m_intersectionThreshold;

    // Perform triangle-triangle tests
    for_all<ExecSpace>(
      numCandidates,
      AXOM_LAMBDA(IndexType i) {
        int index = p_indices[i];
        int candidate = p_validCandidates[i];
        if(primal::intersect(p_tris[index],
                             p_tris[candidate],
                             false,
                             intersectionThreshold))
        {
#ifdef AXOM_USE_RAJA
          auto idx = RAJA::atomicAdd<atomic_pol>(&p_numIsectPairs[0], 1);
#else
          auto idx = p_numIsectPairs[0];
          p_numIsectPairs[0]++;
#endif
          p_firstIsectPair[idx] = index;
          p_secondIsectPair[idx] = candidate;
        }
      });
    axom::copy(&isectCounter, numIsectPairs.data(), sizeof(IndexType));
  }

  firstIsectPair.resize(isectCounter);
  secondIsectPair.resize(isectCounter);

  {
    // copy results to output
    firstIndex = firstIsectPair;
    secondIndex = secondIsectPair;
    HostIndexArray host_degenerate = m_degenerate;
    for(int i = 0; i < host_degenerate.size(); i++)
    {
      if(host_degenerate[i] == 1)
      {
        degenerateIndices.push_back(i);
      }
    }
  }
}

/*!
 * \brief Specialization of CandidateFinder using a Bounding Volume Hierarchy
 *  (BVH) to perform broad-phase collision detection.
 */
template <typename ExecSpace, typename FloatType>
struct CandidateFinder<AccelType::BVH, ExecSpace, FloatType>
  : public CandidateFinderBase<ExecSpace, FloatType>
{
  using BaseClass = CandidateFinderBase<ExecSpace, FloatType>;
  using BaseClass::CandidateFinderBase;
  using BaseClass::HostSpace;
  using BaseClass::Space;

  virtual axom::ArrayView<IndexType, 1, Space> getCandidates(
    axom::Array<IndexType, 1, Space>& offsets,
    axom::Array<IndexType, 1, Space>& counts) override
  {
    int allocatorId = axom::detail::getAllocatorID<Space>();
    spin::BVH<3, ExecSpace, FloatType> bvh;
    bvh.setAllocatorID(allocatorId);
    bvh.initialize(this->m_aabbs.view(), this->m_aabbs.size());

    offsets.resize(this->m_aabbs.size());
    counts.resize(this->m_aabbs.size());

    // Search for intersecting bounding boxes of triangles
    IndexType* candidatesData = nullptr;
    IndexType ncandidates = bvh.findBoundingBoxes(offsets.data(),
                                                  counts.data(),
                                                  candidatesData,
                                                  this->m_aabbs.size(),
                                                  this->m_aabbs.view());

    m_currCandidates.reset(candidatesData);
    return axom::ArrayView<IndexType, 1, Space>(candidatesData, ncandidates);
  }

  // Helper functor used by unique_ptr for deleting the result array returned
  // from the BVH query
  struct CandidateDeleter
  {
    void operator()(IndexType* ptr) const { axom::deallocate(ptr); }
  };

  std::unique_ptr<IndexType, CandidateDeleter> m_currCandidates {
    nullptr,
    CandidateDeleter {}};
};

/*!
 * \brief Specialization of CandidateFinder using an Implicit Grid data
 *  structure to perform broad-phase collision detection.
 */
template <typename ExecSpace, typename FloatType>
struct CandidateFinder<AccelType::ImplicitGrid, ExecSpace, FloatType>
  : public CandidateFinderBase<ExecSpace, FloatType>
{
  using BaseClass = CandidateFinderBase<ExecSpace, FloatType>;
  using BaseClass::CandidateFinderBase;
  using BaseClass::HostSpace;
  using BaseClass::Space;
  using typename BaseClass::BoxType;
  using typename BaseClass::PointType;

  void initialize(int spatialIndexResolution)
  {
    BaseClass::initialize();
#ifdef AXOM_USE_RAJA
    using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
    RAJA::ReduceMin<reduce_pol, double> xmin(DBL_MAX), ymin(DBL_MAX),
      zmin(DBL_MAX);
    RAJA::ReduceMax<reduce_pol, double> xmax(DBL_MIN), ymax(DBL_MIN),
      zmax(DBL_MIN);

    // Get the global bounding box.
    mint::for_all_nodes<ExecSpace, mint::xargs::xyz>(
      this->m_surfaceMesh,
      AXOM_LAMBDA(IndexType, double x, double y, double z) {
        xmin.min(x);
        xmax.max(x);
        ymin.min(y);
        ymax.max(y);
        zmin.min(z);
        zmax.max(z);
      });

    m_globalBox = BoxType(PointType {xmin.get(), ymin.get(), zmin.get()},
                          PointType {xmax.get(), ymax.get(), zmax.get()});
#else
    BoxType global_box;

    // Get the global bounding box.
    mint::for_all_nodes<ExecSpace, mint::xargs::xyz>(
      surface_mesh,
      [=, &global_box](IndexType, double x, double y, double z) {
        global_box.addPoint(PointType {x, y, z});
      });

    // Slightly scale the box
    global_box.scale(1.0001);
    m_globalBox = global_box;
#endif

    // find the specified resolution.  If we're passed a number less than one,
    // use the cube root of the number of triangles.
    if(spatialIndexResolution < 1)
    {
      spatialIndexResolution = (int)(1 + std::pow(this->m_aabbs.size(), 1 / 3.));
    }
    m_resolutions = axom::primal::Point<int, 3>(spatialIndexResolution);
  }

  virtual axom::ArrayView<IndexType, 1, Space> getCandidates(
    axom::Array<IndexType, 1, Space>& offsets,
    axom::Array<IndexType, 1, Space>& counts) override
  {
    int allocatorId = axom::detail::getAllocatorID<Space>();
    axom::spin::ImplicitGrid<3, ExecSpace, IndexType> gridIndex(
      m_globalBox,
      &m_resolutions,
      this->m_aabbs.size(),
      allocatorId);
    gridIndex.insert(this->m_aabbs.size(), this->m_aabbs.data());
    axom::Array<IndexType> offsetsTmp, countsTmp;
    gridIndex.getCandidatesAsArray(this->m_aabbs.size(),
                                   this->m_aabbs.data(),
                                   offsetsTmp,
                                   countsTmp,
                                   m_currCandidates);

    offsets = offsetsTmp;
    counts = countsTmp;
    return m_currCandidates;
  }

  BoxType m_globalBox;
  axom::primal::Point<int, 3> m_resolutions;
  axom::Array<IndexType> m_currCandidates;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom
#endif
