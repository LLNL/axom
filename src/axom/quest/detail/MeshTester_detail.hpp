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
#ifdef AXOM_USE_UMPIRE
  static constexpr bool ExecOnDevice =
    axom::execution_space<ExecSpace>::onDevice();
  static constexpr MemorySpace Space =
    ExecOnDevice ? axom::MemorySpace::Device : axom::MemorySpace::Host;
#else
  static constexpr MemorySpace Space = axom::MemorySpace::Dynamic;
#endif

  CandidateFinderBase(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* surface_mesh,
                      double intersectionThreshold);

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

}  // namespace detail
}  // namespace quest
}  // namespace axom
#endif
