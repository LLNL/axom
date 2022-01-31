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

namespace axom
{
namespace quest
{
namespace detail
{

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

}  // namespace detail
}  // namespace quest
}  // namespace axom
#endif
