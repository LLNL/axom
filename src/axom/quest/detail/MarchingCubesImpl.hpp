// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/marching_cubes_lookup.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "conduit_blueprint.hpp"
#include "axom/fmt.hpp"

#ifndef __WHERE
  #define __STRINGIZE(x) __STRINGIZE2(x)
  #define __STRINGIZE2(x) #x
  //!@brief String literal for code location
  #define __WHERE \
    __FILE__ ":" __STRINGIZE(__LINE__) "(" + std::string(__func__) + ") "
#endif

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
template <typename T, int DIM, typename U>
static void add_to_StackArray(axom::StackArray<T, DIM>& a, U b)
{
  for(int d = 0; d < DIM; ++d) a[d] += b;
}

//!@brief Reverse the order of a StackArray.
template <typename T, int DIM>
static void reverse(axom::StackArray<T, DIM>& a)
{
  for(int d = 0; d < DIM / 2; ++d)
  {
    std::swap(a[d], a[DIM - 1 - d]);
  }
  return;
}

/*!
  @brief Computations for MarchingCubesSingleDomain

  Spatial dimension templating is here, to keep out of higher level
  classes MarchCubes and MarchingCubesSingleDomain.

  Usage:
  @beginverbatim
    MarchingCubesImpl<2, ExecSpace> impl;
    impl.initialize(&domain, coordsetPath, fcnPath, maskPath);
    impl.set_contour_value(contourVal);
    impl.mark_crossings();
    impl.scan_crossings();
    impl.compute_surface();
    impl.populate_surface_mesh(mesh, cellIdField);
  @endverbatim
*/
template <int DIM, typename ExecSpace>
struct MarchingCubesImpl : public MarchingCubesSingleDomain::ImplBase
{
  using Point = axom::primal::Point<double, DIM>;
  using MdimIdx = axom::StackArray<axom::IndexType, DIM>;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  using ReducePolicy = typename execution_space<ExecSpace>::reduce_policy;
  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured domain
    @param coordsetPath Where coordinates are in dom
    @param fcnPath Where nodal function is in dom
    @param maskPath Where cell mask function is in dom

    Set up views to domain data and allocate other data to work on the
    given domain.
  */
  AXOM_HOST void initialize(const conduit::Node& dom,
                            const std::string& coordsetPath,
                            const std::string& fcnPath,
                            const std::string& maskPath) override
  {
    clear();

    // Data sizes
    const conduit::Node& dimsNode =
      dom.fetch_existing("topologies/mesh/elements/dims");
    for(int d = 0; d < DIM; ++d)
    {
      m_bShape[d] = dimsNode[d].as_int();
    }
    m_cShape = m_bShape;
    reverse(m_cShape);
    // This should work but breaks gcc11 on 64-bit linux:
    // m_pShape = m_cShape + 1;
    m_pShape = m_cShape;
    add_to_StackArray(m_pShape, 1);

    // Domain's node coordinates
    {
      const conduit::Node& coordValues =
        dom.fetch_existing(coordsetPath + "/values");
      const bool isInterleaved =
        conduit::blueprint::mcarray::is_interleaved(coordValues);
      const int coordSp = isInterleaved ? DIM : 1;
      for(int d = 0; d < DIM; ++d)
      {
        const double* coordsPtr = coordValues[d].as_double_ptr();
        m_coordsViews[d] =
          axom::ArrayView<const double, DIM>(coordsPtr, m_pShape, coordSp);
      }
    }

    // Nodal function
    {
      auto& fcnValues = dom.fetch_existing(fcnPath + "/values");
      const double* fcnPtr = fcnValues.as_double_ptr();
      m_fcnView = axom::ArrayView<const double, DIM>(fcnPtr, m_pShape);
    }

    // Mask
    {
      const int* maskPtr = nullptr;
      if(!maskPath.empty())
      {
        auto& maskValues = dom.fetch_existing(maskPath + "/values");
        maskPtr = maskValues.as_int_ptr();
      }
      if(maskPtr)
      {
        m_maskView = axom::ArrayView<const int, DIM>(maskPtr, m_cShape);
      }
    }

    m_caseIds =
      axom::Array<std::uint16_t, DIM, axom::MemorySpace::Dynamic>(m_cShape);
  }

  /*!
    @brief Implementation of virtual mark_crossings.

    Virtual methods cannot be templated, so this implementation
    delegates to methods templated on DIM.
  */
  void mark_crossings() override { mark_crossings_dim(); }

// RAJA_INDEX_VALUE_T doesn't work on rzansel.  Maybe old RAJA?
#define USE_RAJA_INDEX_VALUE_TYPES 0

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type mark_crossings_dim()
  {
#if defined(AXOM_USE_RAJA)
  #if USE_RAJA_INDEX_VALUE_TYPES == 1
    RAJA_INDEX_VALUE_T(IIDX, axom::IndexType, "IIDX");
    RAJA_INDEX_VALUE_T(JIDX, axom::IndexType, "JIDX");
    RAJA::TypedRangeSegment<JIDX> jRange(0, m_cShape[0]);
    RAJA::TypedRangeSegment<IIDX> iRange(0, m_cShape[1]);
  #else
    RAJA::RangeSegment jRange(0, m_cShape[0]);
    RAJA::RangeSegment iRange(0, m_cShape[1]);
  #endif
    using EXEC_POL = RAJA::KernelPolicy<
      // TODO: Why does LoopPolicy here not compile on rzansel with cuda?  I have to use RAJA::seq_exec
      RAJA::statement::
        For<1, LoopPolicy, RAJA::statement::For<0, LoopPolicy, RAJA::statement::Lambda<0>>>>;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange),
  #if USE_RAJA_INDEX_VALUE_TYPES == 1
      // TODO: Why does AXOM_LAMBDA here not compile on rzansel with cuda?
      AXOM_LAMBDA(IIDX it, JIDX jt)
  #else
      AXOM_LAMBDA(int i, int j)
  #endif
      {
  #if USE_RAJA_INDEX_VALUE_TYPES == 1
        auto& i = *it;
        auto& j = *jt;
  #endif
        const bool skipZone = !m_maskView.empty() && bool(m_maskView(j, i));
        if(!skipZone)
        {
          // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {m_fcnView(j    , i    ),
             m_fcnView(j    , i + 1),
             m_fcnView(j + 1, i + 1),
             m_fcnView(j + 1, i    )};
          // clang-format on

          auto crossingCase = compute_crossing_case(nodalValues);
          m_caseIds(j, i) = crossingCase;
        }
      });
#else
    for(int j = 0; j < m_cShape[0]; ++j)
    {
      for(int i = 0; i < m_cShape[1]; ++i)
      {
        const bool skipZone = !m_maskView.empty() && bool(m_maskView(i, j));
        if(!skipZone)
        {
          // clang-format on
          double nodalValues[CELL_CORNER_COUNT] = {m_fcnView(j, i),
                                                   m_fcnView(j, i + 1),
                                                   m_fcnView(j + 1, i + 1),
                                                   m_fcnView(j + 1, i)};
          // clang-format off

          auto crossingCase = compute_crossing_case(nodalValues);
          SLIC_ASSERT(crossingCase >= 0 && crossingCase < 16);
          m_caseIds(j,i) = crossingCase;
        }
      }
    }
#endif
  }
  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type mark_crossings_dim()
  {
#if defined(AXOM_USE_RAJA)
#if USE_RAJA_INDEX_VALUE_TYPES == 1
    RAJA_INDEX_VALUE_T(IIDX, axom::IndexType, "IIDX");
    RAJA_INDEX_VALUE_T(JIDX, axom::IndexType, "JIDX");
    RAJA_INDEX_VALUE_T(KIDX, axom::IndexType, "KIDX");
    RAJA::TypedRangeSegment<KIDX> kRange(0, m_cShape[0]);
    RAJA::TypedRangeSegment<JIDX> jRange(0, m_cShape[1]);
    RAJA::TypedRangeSegment<IIDX> iRange(0, m_cShape[2]);
#else
    RAJA::RangeSegment kRange(0, m_cShape[0]);
    RAJA::RangeSegment jRange(0, m_cShape[1]);
    RAJA::RangeSegment iRange(0, m_cShape[2]);
#endif
    using EXEC_POL = RAJA::KernelPolicy<RAJA::statement::For<
      2,
      LoopPolicy,
      RAJA::statement::
        For<1, LoopPolicy, RAJA::statement::For<0, LoopPolicy, RAJA::statement::Lambda<0>>>>>;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange, kRange),
#if USE_RAJA_INDEX_VALUE_TYPES == 1
      AXOM_LAMBDA(IIDX it, JIDX jt, KIDX kt)
#else
      AXOM_LAMBDA(int i, int j, int k)
#endif
      {
#if USE_RAJA_INDEX_VALUE_TYPES == 1
        auto& i = *it;
        auto& j = *jt;
        auto& k = *kt;
#endif
        const bool skipZone = !m_maskView.empty() && bool(m_maskView(k, j, i));
        if(!skipZone)
        {
          // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {m_fcnView(k    , j    , i + 1),
             m_fcnView(k    , j + 1, i + 1),
             m_fcnView(k    , j + 1, i    ),
             m_fcnView(k    , j    , i    ),
             m_fcnView(k + 1, j    , i + 1),
             m_fcnView(k + 1, j + 1, i + 1),
             m_fcnView(k + 1, j + 1, i    ),
             m_fcnView(k + 1, j    , i    )};
          // clang-format on

          auto crossingCase = compute_crossing_case(nodalValues);
          m_caseIds(k, j, i) = crossingCase;
        }
      });
#else
    for(int k = 0; k < m_cShape[0]; ++k)
    {
      for(int j = 0; j < m_cShape[1]; ++j)
      {
        for(int i = 0; i < m_cShape[2]; ++i)
        {
          const bool skipZone = !m_maskView.empty() && bool(m_maskView(k, j, i));
          if(!skipZone)
          {
            // clang-format off
            double nodalValues[CELL_CORNER_COUNT] =
              { m_fcnView(k  , j  , i+1)
              , m_fcnView(k  , j+1, i+1)
              , m_fcnView(k  , j+1, i  )
              , m_fcnView(k  , j  , i  )
              , m_fcnView(k+1, j  , i+1)
              , m_fcnView(k+1, j+1, i+1)
              , m_fcnView(k+1, j+1, i  )
              , m_fcnView(k+1, j  , i  ) };
            // clang-format on

            auto crossingCase = compute_crossing_case(nodalValues);
            SLIC_ASSERT(crossingCase >= 0 && crossingCase < 256);
            m_caseIds(k, j, i) = crossingCase;
          }
        }
      }
    }
#endif
  }

  /*!
    @brief Populate the 1D m_crossings array, one entry for each parent cell that
    crosses the surface.  We sum up the number of surface cells from the crossings,
    allocate space, then populate it.
  */
  void scan_crossings() override
  {
    const axom::IndexType parentCellCount = m_caseIds.size();

#if defined(AXOM_USE_LAMBDA)
    RAJA::ReduceSum<ReducePolicy, axom::IndexType> vsum(0);
    RAJA::forall<LoopPolicy>(
      RAJA::RangeSegment(0, parentCellCount),
      AXOM_LAMBDA(RAJA::Index_type n) {
        vsum += bool(num_surface_cells(m_caseIds.flatIndex(n)));
      });
    m_crossingCount = static_cast<axom::IndexType>(vsum.get());
#else
    axom::IndexType vsum = 0;
    for(axom::IndexType n = 0; n < parentCellCount; ++n)
    {
      vsum += bool(num_surface_cells(m_caseIds.flatIndex(n)));
    }
    m_crossingCount = vsum;
#endif

    m_crossings.resize(m_crossingCount, {0, 0});
    axom::ArrayView<CrossingInfo> crossingsView = m_crossings.view();

    axom::Array<int> addCells(m_crossingCount,
                              m_crossingCount,
                              m_crossings.getAllocatorID());
    axom::ArrayView<int> addCellsView = addCells.view();

    axom::IndexType* crossingId = axom::allocate<axom::IndexType>(
      1,
      axom::execution_space<ExecSpace>::allocatorID());
    *crossingId = 0;
    /*
      This loop can't be parallelized, even with devices, due
      (*crossingId)'s dependence on the if outcome.  Must use SEQ_EXEC,
      regardless of ExecSpace.

      NOTE: ExecSpace works anyway.  Why?
    */
    axom::for_all<ExecSpace>(
      0,
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType n) {
        auto caseId = m_caseIds.flatIndex(n);
        auto ccc = num_surface_cells(caseId);
        if(ccc != 0)
        {
          addCellsView[*crossingId] = ccc;
          m_crossings[*crossingId].caseNum = caseId;
          m_crossings[*crossingId].parentCellNum = n;
          ++(*crossingId);
        }
      });
    assert(*crossingId == m_crossingCount);
    axom::deallocate(crossingId);

    axom::Array<axom::IndexType> prefixSum(m_crossingCount,
                                           m_crossingCount,
                                           m_crossings.getAllocatorID());
    axom::ArrayView<axom::IndexType> prefixSumView = prefixSum.view();

    auto copyFirstSurfaceCellId = AXOM_LAMBDA(axom::IndexType n)
    {
      crossingsView[n].firstSurfaceCellId = prefixSumView[n];
    };
#if defined(AXOM_USE_RAJA)
    RAJA::exclusive_scan<LoopPolicy>(
      RAJA::make_span(addCellsView.data(), m_crossingCount),
      RAJA::make_span(prefixSumView.data(), m_crossingCount),
      RAJA::operators::plus<axom::IndexType> {});
    RAJA::forall<LoopPolicy>(RAJA::RangeSegment(0, m_crossingCount),
                             copyFirstSurfaceCellId);
#else
    if(m_crossingCount > 0)
    {
      prefixSumView[0] = 0;
      for(axom::IndexType i = 1; i < m_crossingCount; ++i)
      {
        prefixSumView[i] = prefixSumView[i - 1] + addCellsView[i - 1];
      }
      for(axom::IndexType i = 0; i < m_crossingCount; ++i)
      {
        copyFirstSurfaceCellId(i);
      }
    }
#endif

    m_surfaceCellCount = m_crossings.empty()
      ? 0
      : m_crossings.back().firstSurfaceCellId +
        num_surface_cells(m_crossings.back().caseNum);
  }

  void compute_surface() override
  {
    /*
      Reserve surface mesh data space so we can add data without
      reallocation.
    */
    const axom::IndexType surfaceNodeCount = DIM * m_surfaceCellCount;

    m_surfaceCoords.resize(surfaceNodeCount);
    m_surfaceCellCorners.resize(m_surfaceCellCount);
    m_surfaceCellParents.resize(m_surfaceCellCount);

    axom::for_all<ExecSpace>(
      0,
      m_crossingCount,
      AXOM_LAMBDA(axom::IndexType iCrossing) {
        const auto& crossingInfo = m_crossings[iCrossing];
        const IndexType crossingCellCount =
          num_surface_cells(crossingInfo.caseNum);
        SLIC_ASSERT(crossingCellCount > 0);

        // Parent cell data for interpolating new node coordinates.
        Point cornerCoords[CELL_CORNER_COUNT];
        double cornerValues[CELL_CORNER_COUNT];
        get_corner_coords_and_values(crossingInfo.parentCellNum,
                                     cornerCoords,
                                     cornerValues);

        /*
          Create the new cell and its DIM nodes.  New node are on
          parent cell edges where they intersects the isosurface.
          linear_interp for the exact coordinates.

          TODO: The varying crossingCellCount value may inhibit device
          performance.  Try grouping m_crossings items that have the
          same values for crossingCellCount.
        */
        for(int iCell = 0; iCell < crossingCellCount; ++iCell)
        {
          IndexType surfaceCellId = crossingInfo.firstSurfaceCellId + iCell;
          m_surfaceCellParents[surfaceCellId] = crossingInfo.parentCellNum;
          for(int d = 0; d < DIM; ++d)
          {
            IndexType surfaceNodeId = surfaceCellId * DIM + d;
            m_surfaceCellCorners[surfaceCellId][d] = surfaceNodeId;

            const int edge = cases_table(crossingInfo.caseNum, iCell * DIM + d);
            linear_interp(edge,
                          cornerCoords,
                          cornerValues,
                          m_surfaceCoords[surfaceNodeId]);
          }
        }
      });
  }

  // These 4 functions provides access to the look-up table
  // whether on host or device.  Is there a more elegant way
  // to put static 1D and 2D arrays on host and device?  BTNG.
  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type num_surface_cells(
    int iCase) const
  {
#define _MC_LOOKUP_NUM_SEGMENTS
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_SEGMENTS
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return num_segments[iCase];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type cases_table(
    int iCase,
    int iEdge) const
  {
#define _MC_LOOKUP_CASES2D
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_CASES2D
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return cases2D[iCase][iEdge];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type num_surface_cells(
    int iCase) const
  {
#define _MC_LOOKUP_NUM_TRIANGLES
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_TRIANGLES
    SLIC_ASSERT(iCase >= 0 && iCase < 256);
    return num_triangles[iCase];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type cases_table(
    int iCase,
    int iEdge) const
  {
#define _MC_LOOKUP_CASES3D
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_CASES3D
    SLIC_ASSERT(iCase >= 0 && iCase < 256);
    return cases3D[iCase][iEdge];
  }

  /*!
    @brief Compute multidimensional index from flat cell index
    in domain data.
  */
  AXOM_HOST_DEVICE axom::StackArray<axom::IndexType, DIM> multidim_cell_index(
    axom::IndexType flatId) const
  {
    axom::IndexType strides[DIM] = {1};
    for(int d = 1; d < DIM; ++d) strides[d] = strides[d - 1] * m_bShape[d - 1];

    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d = DIM - 1; d >= 0; --d)
    {
      rval[d] = flatId / strides[d];
      flatId -= rval[d] * strides[d];
    }
    return rval;
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type
  get_corner_coords_and_values(IndexType cellNum,
                               Point cornerCoords[],
                               double cornerValues[])
  {
    const auto c = multidim_cell_index(cellNum);
    const auto& i = c[0];
    const auto& j = c[1];

    const auto& x = m_coordsViews[0];
    const auto& y = m_coordsViews[1];

    // clang-format off
    cornerCoords[0] = { x(j  ,i  ), y(j  ,i  ) };
    cornerCoords[1] = { x(j  ,i+1), y(j  ,i+1) };
    cornerCoords[2] = { x(j+1,i+1), y(j+1,i+1) };
    cornerCoords[3] = { x(j+1,i  ), y(j+1,i  ) };

    cornerValues[0] = m_fcnView(j  ,i  );
    cornerValues[1] = m_fcnView(j  ,i+1);
    cornerValues[2] = m_fcnView(j+1,i+1);
    cornerValues[3] = m_fcnView(j+1,i  );
    // clang-format on
  }
  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type
  get_corner_coords_and_values(IndexType cellNum,
                               Point cornerCoords[],
                               double cornerValues[])
  {
    const auto c = multidim_cell_index(cellNum);
    const auto& i = c[0];
    const auto& j = c[1];
    const auto& k = c[2];

    const auto& x = m_coordsViews[0];
    const auto& y = m_coordsViews[1];
    const auto& z = m_coordsViews[2];

    // clang-format off
    cornerCoords[0] = { x(k  ,j  ,i+1), y(k  ,j  ,i+1), z(k  ,j  ,i+1) };
    cornerCoords[1] = { x(k  ,j+1,i+1), y(k  ,j+1,i+1), z(k  ,j+1,i+1) };
    cornerCoords[2] = { x(k  ,j+1,i  ), y(k  ,j+1,i  ), z(k  ,j+1,i  ) };
    cornerCoords[3] = { x(k  ,j  ,i  ), y(k  ,j  ,i  ), z(k  ,j  ,i  ) };
    cornerCoords[4] = { x(k+1,j  ,i+1), y(k+1,j  ,i+1), z(k+1,j  ,i+1) };
    cornerCoords[5] = { x(k+1,j+1,i+1), y(k+1,j+1,i+1), z(k+1,j+1,i+1) };
    cornerCoords[6] = { x(k+1,j+1,i  ), y(k+1,j+1,i  ), z(k+1,j+1,i  ) };
    cornerCoords[7] = { x(k+1,j  ,i  ), y(k+1,j  ,i  ), z(k+1,j  ,i  ) };

    cornerValues[0] = m_fcnView(k  ,j  ,i+1);
    cornerValues[1] = m_fcnView(k  ,j+1,i+1);
    cornerValues[2] = m_fcnView(k  ,j+1,i  );
    cornerValues[3] = m_fcnView(k  ,j  ,i  );
    cornerValues[4] = m_fcnView(k+1,j  ,i+1);
    cornerValues[5] = m_fcnView(k+1,j+1,i+1);
    cornerValues[6] = m_fcnView(k+1,j+1,i  );
    cornerValues[7] = m_fcnView(k+1,j  ,i  );
    // clang-format on
  }

  void populate_surface_mesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
    const std::string& cellIdField) const override
  {
    if(!cellIdField.empty() &&
       !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
    {
      mesh.createField<axom::IndexType>(cellIdField,
                                        axom::mint::CELL_CENTERED,
                                        DIM);
    }

    const axom::IndexType addedCellCount = m_surfaceCellCorners.size();
    const axom::IndexType addedNodeCount = m_surfaceCoords.size();
    if(addedCellCount != 0)
    {
      const axom::IndexType priorCellCount = mesh.getNumberOfCells();
      const axom::IndexType priorNodeCount = mesh.getNumberOfNodes();
      mesh.reserveNodes(priorNodeCount + addedNodeCount);
      mesh.reserveCells(priorCellCount + addedCellCount);

      mesh.appendNodes((double*)m_surfaceCoords.data(), m_surfaceCoords.size());
      for(int n = 0; n < addedCellCount; ++n)
      {
        MdimIdx cornerIds = m_surfaceCellCorners[n];
        add_to_StackArray(cornerIds, priorNodeCount);
        mesh.appendCell(cornerIds);
      }
      axom::IndexType numComponents = -1;
      axom::IndexType* dstPtr =
        mesh.getFieldPtr<axom::IndexType>(cellIdField,
                                          axom::mint::CELL_CENTERED,
                                          numComponents);
      SLIC_ASSERT(numComponents == DIM);
      axom::ArrayView<axom::StackArray<axom::IndexType, DIM>> dstView(
        (axom::StackArray<axom::IndexType, DIM>*)dstPtr,
        priorCellCount + addedCellCount);
      for(axom::IndexType i = 0; i < addedCellCount; ++i)
      {
        dstView[priorCellCount + i] =
          multidim_cell_index(m_surfaceCellParents[i]);
      }
    }
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type linear_interp(
    int edgeIdx,
    const Point cornerCoords[4],
    const double nodeValues[4],
    Point& crossingPt)
  {
    // STEP 0: get the edge node indices
    // 2 nodes define the edge.  n1 and n2 are the indices of
    // the nodes w.r.t. the square or cubic zone.  There is a
    // agreed-on ordering of these indices in the arrays xx, yy,
    // zz, nodeValues, crossingPt.
    int n1 = edgeIdx;
    int n2 = (edgeIdx == 3) ? 0 : edgeIdx + 1;

    // STEP 1: get the fields and coordinates from the two points
    const double f1 = nodeValues[n1];
    const double f2 = nodeValues[n2];

    const Point& p1 = cornerCoords[n1];
    const Point& p2 = cornerCoords[n2];

    // STEP 2: check whether the interpolated point is at one of the two corners.
    if(axom::utilities::isNearlyEqual(m_contourVal, f1) ||
       axom::utilities::isNearlyEqual(f1, f2))
    {
      crossingPt = p1;  // memcpy(crossingPt, p1, DIM * sizeof(double));
      return;
    }

    if(axom::utilities::isNearlyEqual(m_contourVal, f2))
    {
      crossingPt = p2;  // memcpy(crossingPt, p2, DIM * sizeof(double));
      return;
    }

    // STEP 3: point is in between the edge points, interpolate its position
    constexpr double ptiny = 1.0e-80;
    const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
    const double w = (m_contourVal - f1) / df;
    for(int d = 0; d < DIM; ++d)
    {
      crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
    }
  }
  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type linear_interp(
    int edgeIdx,
    const Point cornerCoords[8],
    const double nodeValues[8],
    Point& crossingPt)
  {
    // STEP 0: get the edge node indices
    // 2 nodes define the edge.  n1 and n2 are the indices of
    // the nodes w.r.t. the square or cubic zone.  There is a
    // agreed-on ordering of these indices in the arrays
    // cornerCoords, nodeValues, hex_edge_table.
    const int hex_edge_table[] = {
      0, 1, 1, 2, 2, 3, 3, 0,  // base
      4, 5, 5, 6, 6, 7, 7, 4,  // top
      0, 4, 1, 5, 2, 6, 3, 7   // vertical
    };

    int n1 = hex_edge_table[edgeIdx * 2];
    int n2 = hex_edge_table[edgeIdx * 2 + 1];

    // STEP 1: get the fields and coordinates from the two points
    const double f1 = nodeValues[n1];
    const double f2 = nodeValues[n2];

    const Point& p1 = cornerCoords[n1];
    const Point& p2 = cornerCoords[n2];

    // STEP 2: check whether the interpolated point is at one of the two corners.
    if(axom::utilities::isNearlyEqual(m_contourVal, f1) ||
       axom::utilities::isNearlyEqual(f1, f2))
    {
      crossingPt = p1;
      return;
    }

    if(axom::utilities::isNearlyEqual(m_contourVal, f2))
    {
      crossingPt = p2;
      return;
    }

    // STEP 3: point is not at corner; interpolate its position
    constexpr double ptiny = 1.0e-80;
    const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
    const double w = (m_contourVal - f1) / df;
    for(int d = 0; d < DIM; ++d)
    {
      crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
    }
  }

  void set_contour_value(double contourVal) override
  {
    m_contourVal = contourVal;
  }

  //!@brief Compute the case index into cases2D or cases3D.
  AXOM_HOST_DEVICE int compute_crossing_case(const double* f)
  {
    int index = 0;
    for(int n = 0; n < CELL_CORNER_COUNT; ++n)
    {
      if(f[n] >= m_contourVal)
      {
        const int bit = (1 << n);
        index |= bit;
      }
    }
    return index;
  }

  //!@brief Clear data so you can rerun with a different contour value.
  void clear()
  {
    m_surfaceCoords.clear();
    m_surfaceCellCorners.clear();
    m_surfaceCellParents.clear();
    m_crossingCount = 0;
    m_surfaceCellCount = 0;
  }

  /*
    Internal data arrays use memory compatible with memory for ExecSpace.
  */
  MarchingCubesImpl()
    : m_crossings(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_surfaceCoords(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_surfaceCellCorners(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_surfaceCellParents(0, 0, execution_space<ExecSpace>::allocatorID())
  { }

  /*!
    @brief Info for a parent cell intersecting the surface.
  */
  struct CrossingInfo
  {
    CrossingInfo(axom::IndexType parentCellNum_, std::uint16_t caseNum_)
      : parentCellNum(parentCellNum_)
      , caseNum(caseNum_)
      , firstSurfaceCellId(std::numeric_limits<axom::IndexType>::max())
    { }
    axom::IndexType parentCellNum;       //!< @brief Flat index of parent cell.
    std::uint16_t caseNum;               //!< @brief Index in cases2D or cases3D
    axom::IndexType firstSurfaceCellId;  //!< @brief First index for generated cells.
  };
  axom::Array<CrossingInfo> m_crossings;

  MdimIdx m_bShape;  //!< @brief Blueprint cell data shape.
  MdimIdx m_cShape;  //!< @brief Cell-centered array shape for ArrayViews.
  MdimIdx m_pShape;  //!< @brief Node-centered array shape for ArrayViews.

  // Views of parent domain data.
  axom::ArrayView<const double, DIM> m_coordsViews[DIM];
  axom::ArrayView<const double, DIM> m_fcnView;
  axom::ArrayView<const int, DIM> m_maskView;

  //!@brief Crossing case for each computational mesh cell.
  // TODO: Put this in correct memory space.
  axom::Array<std::uint16_t, DIM, axom::MemorySpace::Dynamic> m_caseIds;

  //!@brief Number of parent cells crossing the contour surface.
  axom::IndexType m_crossingCount = 0;
  //!@brief Number of surface cells from crossings.
  axom::IndexType m_surfaceCellCount = 0;
  axom::IndexType get_surface_cell_count() const override
  {
    return m_surfaceCellCount;
  }

  //!@brief Number of corners (nodes) on each cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;

  //!@name Internal representation of generated surface mesh.
  //@{
  //!@brief Coordinates of generated surface nodes.
  axom::Array<Point> m_surfaceCoords;

  //!@brief Corners (index into m_surfaceCoords) of generated surface cells.
  axom::Array<MdimIdx> m_surfaceCellCorners;

  //!@brief Computational cell (flat index) crossing the surface cell.
  axom::Array<IndexType> m_surfaceCellParents;
  //@}

  double m_contourVal = 0.0;
};

}  // end namespace marching_cubes
}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
