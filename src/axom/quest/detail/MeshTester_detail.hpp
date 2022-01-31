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

  CandidateFinderBase(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* surface_mesh,
                      double intersectionThreshold);

  void findTriMeshIntersections(axom::Array<IndexType>& firstIndex,
                                axom::Array<IndexType>& secondIndex,
                                axom::Array<IndexType>& degenerateIndices);

  virtual void getCandidates(axom::Array<IndexType, 1, Space>& offsets,
                             axom::Array<IndexType, 1, Space>& counts,
                             axom::Array<IndexType, 1, Space>& candidates) = 0;

  mint::UnstructuredMesh<mint::SINGLE_SHAPE>* m_surfaceMesh;
  double m_intersectionThreshold;
  int m_ncells;
  axom::Array<detail::Triangle3, 1, Space> m_tris;
  axom::Array<BoxType, 1, Space> m_aabbs;
  axom::Array<int, 1, Space> m_degenerate;
};

template <typename ExecSpace, typename FloatType>
CandidateFinderBase<ExecSpace, FloatType>::CandidateFinderBase(
  mint::UnstructuredMesh<mint::SINGLE_SHAPE>* surface_mesh,
  double intersectionThreshold)
  : m_surfaceMesh(surface_mesh)
  , m_intersectionThreshold(intersectionThreshold)
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
    surface_mesh,
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
  IndexArray offsets, counts, candidates;
  getCandidates(offsets, counts, candidates);

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
    IndexView p_candidates = candidates;

    // Initialize triangle indices and valid candidates
    for_all<ExecSpace>(
      ncells,
      AXOM_LAMBDA(IndexType i) {
        for(int j = 0; j < p_counts[i]; j++)
        {
          if(i < p_candidates[p_offsets[i] + j])
          {
#ifdef AXOM_USE_RAJA
            auto idx = RAJA::atomicAdd<atomic_pol>(&p_numValidCandidates[0], 1);
#else
            auto idx = p_numValidCandidates[0]++;
#endif
            p_indices[idx] = i;
            p_validCandidates[idx] = p_candidates[p_offsets[i] + j];
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

template <typename ExecSpace, typename FloatType>
struct CandidateFinder<AccelType::BVH, ExecSpace, FloatType>
  : public CandidateFinderBase<ExecSpace, FloatType>
{
  using BaseClass = CandidateFinderBase<ExecSpace, FloatType>;
  using BaseClass::CandidateFinderBase;
  using BaseClass::HostSpace;
  using BaseClass::Space;

  virtual void getCandidates(axom::Array<IndexType, 1, Space>& offsets,
                             axom::Array<IndexType, 1, Space>& counts,
                             axom::Array<IndexType, 1, Space>& candidates) override
  {
    int allocatorId = axom::detail::getAllocatorID<Space>();
    spin::BVH<3, ExecSpace, FloatType> bvh;
    bvh.setAllocatorID(allocatorId);
    bvh.initialize(this->m_aabbs.view(), this->m_aabbs.size());

    offsets.resize(this->m_aabbs.size());
    counts.resize(this->m_aabbs.size());

    // Search for intersecting bounding boxes of triangles
    IndexType* candidatesData = nullptr;
    bvh.findBoundingBoxes(offsets.data(),
                          counts.data(),
                          candidatesData,
                          this->m_aabbs.size(),
                          this->m_aabbs.view());

    IndexType ncandidates;
    {
      axom::Array<IndexType, 1, Space> ncandidates_buf(1);
      axom::ArrayView<IndexType, 1, Space> p_ncandidates = ncandidates_buf;
      axom::ArrayView<IndexType, 1, Space> p_offsets = offsets;
      axom::ArrayView<IndexType, 1, Space> p_counts = counts;
      IndexType lastIdx = offsets.size() - 1;
      for_all<ExecSpace>(
        1,
        AXOM_LAMBDA(IndexType) {
          p_ncandidates[0] = p_offsets[lastIdx] + p_counts[lastIdx];
        });
      axom::copy(&ncandidates, ncandidates_buf.data(), sizeof(IndexType));
    }
    axom::ArrayView<IndexType, 1, Space> candidateBuf(candidatesData,
                                                      ncandidates);
    candidates = candidateBuf;
    axom::deallocate(candidatesData);
  }
};

template <typename ExecSpace, typename FloatType>
struct CandidateFinder<AccelType::ImplicitGrid, ExecSpace, FloatType>
  : public CandidateFinderBase<ExecSpace, FloatType>
{
  using BaseClass = CandidateFinderBase<ExecSpace, FloatType>;
  using BaseClass::HostSpace;
  using BaseClass::Space;
  using typename BaseClass::BoxType;
  using typename BaseClass::PointType;

  CandidateFinder(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* surface_mesh,
                  int spatialIndexResolution,
                  double intersectionThreshold)
    : BaseClass(surface_mesh, intersectionThreshold)
  {
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

  virtual void getCandidates(axom::Array<IndexType, 1, Space>& offsets,
                             axom::Array<IndexType, 1, Space>& counts,
                             axom::Array<IndexType, 1, Space>& candidates) override
  {
    int allocatorId = axom::detail::getAllocatorID<Space>();
    axom::spin::ImplicitGrid<3, ExecSpace, IndexType> gridIndex(
      m_globalBox,
      &m_resolutions,
      this->m_aabbs.size(),
      allocatorId);
    gridIndex.insert(this->m_aabbs.size(), this->m_aabbs.data());
    axom::Array<IndexType> offsetsTmp, countsTmp, candidatesTmp;
    gridIndex.getCandidatesAsArray(this->m_aabbs.size(),
                                   this->m_aabbs.data(),
                                   offsetsTmp,
                                   countsTmp,
                                   candidatesTmp);

    offsets = offsetsTmp;
    counts = countsTmp;
    candidates = candidatesTmp;
  }

  BoxType m_globalBox;
  axom::primal::Point<int, 3> m_resolutions;
};

}  // namespace detail
}  // namespace quest
}  // namespace axom
#endif
