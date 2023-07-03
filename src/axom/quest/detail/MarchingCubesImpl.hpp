// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/marching_cubes_lookup.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/constants.hpp"
#include "axom/mint/execution/internal/structured_exec.hpp"
#include "conduit_blueprint.hpp"
#include "axom/fmt.hpp"

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

  ExecSpace is the general execution space, like axom::SEQ_EXEC and
  axom::CUDA_EXEC<256>.  SequentialExecSpace is used for loops that
  cannot be parallelized but must access data allocated for ExecSpace.
  Use something like axom::SEQ_EXEC or axom::CUDA_EXEC<1>.
*/
template <int DIM, typename ExecSpace, typename SequentialExecSpace>
class MarchingCubesImpl : public MarchingCubesSingleDomain::ImplBase
{
public:
  using Point = axom::primal::Point<double, DIM>;
  using MIdx = axom::StackArray<axom::IndexType, DIM>;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  using ReducePolicy = typename execution_space<ExecSpace>::reduce_policy;
  using SequentialLoopPolicy =
    typename execution_space<SequentialExecSpace>::loop_policy;
  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;
  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured mesh domain
    @param coordsetPath Where coordinates are in dom
    @param fcnPath Where nodal function is in dom
    @param maskPath Where cell mask function is in dom

    Set up views to domain data and allocate other data to work on the
    given domain.

    The above data from the domain MUST be in a memory space
    compatible with ExecSpace.
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

    m_bStrides[0] = 1;
    for(int d = 1; d < DIM; ++d)
      m_bStrides[d] = m_bStrides[d - 1] * m_bShape[d - 1];

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
          axom::ArrayView<const double, DIM, MemorySpace>(coordsPtr,
                                                          m_pShape,
                                                          coordSp);
      }
    }

    // Nodal function
    {
      auto& fcnValues = dom.fetch_existing(fcnPath + "/values");
      const double* fcnPtr = fcnValues.as_double_ptr();
      m_fcnView =
        axom::ArrayView<const double, DIM, MemorySpace>(fcnPtr, m_pShape);
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
        m_maskView =
          axom::ArrayView<const int, DIM, MemorySpace>(maskPtr, m_cShape);
      }
    }

    m_caseIds = axom::Array<std::uint16_t, DIM, MemorySpace>(m_cShape);
  }

  //!@brief Set a value to find the contour for.
  void setContourValue(double contourVal) override
  {
    m_contourVal = contourVal;
  }

  /*!
    @brief Implementation of virtual markCrossings.

    Virtual methods cannot be templated, so this implementation
    delegates to a name templated on DIM.
  */
  void markCrossings() override { markCrossings_dim(); }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type markCrossings_dim()
  {
    m_caseIds.fill(0);

    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal);

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment jRange(0, m_cShape[0]);
    RAJA::RangeSegment iRange(0, m_cShape[1]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop2d_policy;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        mcu.computeCaseId(i, j);
      });
#else
    for(int j = 0; j < m_cShape[0]; ++j)
    {
      for(int i = 0; i < m_cShape[1]; ++i)
      {
        mcu.computeCaseId(i, j);
      }
    }
#endif
  }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type markCrossings_dim()
  {
    m_caseIds.fill(0);

    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal);

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment kRange(0, m_cShape[0]);
    RAJA::RangeSegment jRange(0, m_cShape[1]);
    RAJA::RangeSegment iRange(0, m_cShape[2]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop3d_policy;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange, kRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        mcu.computeCaseId(i, j, k);
      });
#else
    for(int k = 0; k < m_cShape[0]; ++k)
    {
      for(int j = 0; j < m_cShape[1]; ++j)
      {
        for(int i = 0; i < m_cShape[2]; ++i)
        {
          mcu.computeCaseId(i, j, k);
        }
      }
    }
#endif
  }

  /*!
    @brief Implementation used by MarchingCubesImpl::markCrossings_dim()
    containing just the objects needed for that part, to be made available
    on devices.
  */
  struct MarkCrossings_Util
  {
    axom::ArrayView<std::uint16_t, DIM, MemorySpace> caseIdsView;
    axom::ArrayView<const double, DIM, MemorySpace> fcnView;
    axom::ArrayView<const int, DIM, MemorySpace> maskView;
    double contourVal;
    MarkCrossings_Util(axom::Array<std::uint16_t, DIM, MemorySpace>& caseIds,
                       axom::ArrayView<const double, DIM, MemorySpace>& fcnView_,
                       axom::ArrayView<const int, DIM, MemorySpace>& maskView_,
                       double contourVal_)
      : caseIdsView(caseIds.view())
      , fcnView(fcnView_)
      , maskView(maskView_)
      , contourVal(contourVal_)
    { }

    //!@brief Compute the case index into cases2D or cases3D.
    AXOM_HOST_DEVICE inline int computeCrossingCase(const double* f) const
    {
      int index = 0;
      for(int n = 0; n < CELL_CORNER_COUNT; ++n)
      {
        if(f[n] >= contourVal)
        {
          const int bit = (1 << n);
          index |= bit;
        }
      }
      return index;
    }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2>::type
    computeCaseId(axom::IndexType i, axom::IndexType j) const
    {
      const bool skipZone = !maskView.empty() && bool(maskView(j, i));
      if(!skipZone)
      {
        // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {fcnView(j    , i    ),
             fcnView(j    , i + 1),
             fcnView(j + 1, i + 1),
             fcnView(j + 1, i    )};
        // clang-format on
        caseIdsView(j, i) = computeCrossingCase(nodalValues);
      }
    }

    //!@brief Populate m_caseIds with crossing indices.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3>::type
    computeCaseId(axom::IndexType i, axom::IndexType j, axom::IndexType k) const
    {
      const bool skipZone = !maskView.empty() && bool(maskView(k, j, i));
      if(!skipZone)
      {
        // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {fcnView(k    , j    , i + 1),
             fcnView(k    , j + 1, i + 1),
             fcnView(k    , j + 1, i    ),
             fcnView(k    , j    , i    ),
             fcnView(k + 1, j    , i + 1),
             fcnView(k + 1, j + 1, i + 1),
             fcnView(k + 1, j + 1, i    ),
             fcnView(k + 1, j    , i    )};
        // clang-format on
        caseIdsView(k, j, i) = computeCrossingCase(nodalValues);
      }
    }
  };  // MarkCrossings_Util

  /*!
    @brief Populate the 1D m_crossings array, one entry for each
    parent cell that crosses the contour.

    We sum up the number of contour surface cells from the crossings,
    allocate space, then populate it.
  */
  void scanCrossings() override
  {
    const axom::IndexType parentCellCount = m_caseIds.size();
    auto caseIdsView = m_caseIds.view();
#if defined(AXOM_USE_RAJA)
    RAJA::ReduceSum<ReducePolicy, axom::IndexType> vsum(0);
    RAJA::forall<LoopPolicy>(
      RAJA::RangeSegment(0, parentCellCount),
      AXOM_LAMBDA(RAJA::Index_type n) {
        vsum += bool(num_contour_cells(caseIdsView.flatIndex(n)));
      });
    m_crossingCount = static_cast<axom::IndexType>(vsum.get());
#else
    axom::IndexType vsum = 0;
    for(axom::IndexType n = 0; n < parentCellCount; ++n)
    {
      vsum += bool(num_contour_cells(caseIdsView.flatIndex(n)));
    }
    m_crossingCount = vsum;
#endif

    m_crossings.resize(m_crossingCount, {0, 0});
    axom::ArrayView<CrossingInfo, 1, MemorySpace> crossingsView =
      m_crossings.view();

    axom::Array<int, 1, MemorySpace> addCells(m_crossingCount, m_crossingCount);
    const axom::ArrayView<int, 1, MemorySpace> addCellsView = addCells.view();

    axom::IndexType* crossingId = axom::allocate<axom::IndexType>(
      1,
      axom::detail::getAllocatorID<MemorySpace>());

    auto loopBody = AXOM_LAMBDA(axom::IndexType n)
    {
      auto caseId = caseIdsView.flatIndex(n);
      auto ccc = num_contour_cells(caseId);
      if(ccc != 0)
      {
        addCellsView[*crossingId] = ccc;
        crossingsView[*crossingId].caseNum = caseId;
        crossingsView[*crossingId].parentCellNum = n;
        ++(*crossingId);
      }
    };

#if defined(AXOM_USE_RAJA)
    /*
      The m_crossings filling loop isn't data-parallel and shouldn't
      be parallelized.  This contrived RAJA::forall forces it to run
      sequentially.
    */
    RAJA::forall<SequentialLoopPolicy>(
      RAJA::RangeSegment(0, 1),
      [=] AXOM_HOST_DEVICE(int /* i */) {
        *crossingId = 0;
        for(axom::IndexType n = 0; n < parentCellCount; ++n)
        {
          loopBody(n);
        }
        // assert(*crossingId == m_crossingCount);
      });
#else
    *crossingId = 0;
    for(axom::IndexType n = 0; n < parentCellCount; ++n)
    {
      loopBody(n);
    }
    SLIC_ASSERT(*crossingId == m_crossingCount);
#endif

    axom::deallocate(crossingId);

    axom::Array<axom::IndexType, 1, MemorySpace> prefixSum(m_crossingCount,
                                                           m_crossingCount);
    const axom::ArrayView<axom::IndexType, 1, MemorySpace> prefixSumView =
      prefixSum.view();

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

    // Data from the last crossing tells us how many contour cells there are.
    if(m_crossings.empty())
    {
      m_contourCellCount = 0;
    }
    else
    {
      CrossingInfo back;
      axom::copy(&back,
                 m_crossings.data() + m_crossings.size() - 1,
                 sizeof(CrossingInfo));
      m_contourCellCount =
        back.firstSurfaceCellId + num_contour_cells(back.caseNum);
    }
  }

  /*!
    @brief Implementation used by MarchingCubesImpl::computeContour().
    containing just the objects needed for that part, to be made available
    on devices.
  */
  struct ComputeContour_Util
  {
    double contourVal;
    MIdx bStrides;
    axom::ArrayView<const double, DIM, MemorySpace> fcnView;
    axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM> coordsViews;
    ComputeContour_Util(
      double contourVal_,
      const MIdx& bStrides_,
      const axom::ArrayView<const double, DIM, MemorySpace>& fcnView_,
      const axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>
        coordsViews_)
      : contourVal(contourVal_)
      , bStrides(bStrides_)
      , fcnView(fcnView_)
      , coordsViews(coordsViews_)
    { }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type
    get_corner_coords_and_values(IndexType cellNum,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];

      const auto c = multidim_cell_index(cellNum);
      const auto& i = c[0];
      const auto& j = c[1];

      // clang-format off
      cornerCoords[0] = { x(j    , i    ), y(j    , i    ) };
      cornerCoords[1] = { x(j    , i + 1), y(j    , i + 1) };
      cornerCoords[2] = { x(j + 1, i + 1), y(j + 1, i + 1) };
      cornerCoords[3] = { x(j + 1, i    ), y(j + 1, i    ) };

      cornerValues[0] = fcnView(j    , i    );
      cornerValues[1] = fcnView(j    , i + 1);
      cornerValues[2] = fcnView(j + 1, i + 1);
      cornerValues[3] = fcnView(j + 1, i    );
      // clang-format on
    }
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type
    get_corner_coords_and_values(IndexType cellNum,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];
      const auto& z = coordsViews[2];

      const auto c = multidim_cell_index(cellNum);
      const auto& i = c[0];
      const auto& j = c[1];
      const auto& k = c[2];

      // clang-format off
      cornerCoords[0] = { x(k  , j  , i+1), y(k  , j  , i+1), z(k  , j  , i+1) };
      cornerCoords[1] = { x(k  , j+1, i+1), y(k  , j+1, i+1), z(k  , j+1, i+1) };
      cornerCoords[2] = { x(k  , j+1, i  ), y(k  , j+1, i  ), z(k  , j+1, i  ) };
      cornerCoords[3] = { x(k  , j  , i  ), y(k  , j  , i  ), z(k  , j  , i  ) };
      cornerCoords[4] = { x(k+1, j  , i+1), y(k+1, j  , i+1), z(k+1, j  , i+1) };
      cornerCoords[5] = { x(k+1, j+1, i+1), y(k+1, j+1, i+1), z(k+1, j+1, i+1) };
      cornerCoords[6] = { x(k+1, j+1, i  ), y(k+1, j+1, i  ), z(k+1, j+1, i  ) };
      cornerCoords[7] = { x(k+1, j  , i  ), y(k+1, j  , i  ), z(k+1, j  , i  ) };

      cornerValues[0] = fcnView(k  , j  , i+1);
      cornerValues[1] = fcnView(k  , j+1, i+1);
      cornerValues[2] = fcnView(k  , j+1, i  );
      cornerValues[3] = fcnView(k  , j  , i  );
      cornerValues[4] = fcnView(k+1, j  , i+1);
      cornerValues[5] = fcnView(k+1, j+1, i+1);
      cornerValues[6] = fcnView(k+1, j+1, i  );
      cornerValues[7] = fcnView(k+1, j  , i  );
      // clang-format on
    }

    //!@brief Interpolate for the contour location crossing a parent edge.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type linear_interp(
      int edgeIdx,
      const Point cornerCoords[4],
      const double nodeValues[4],
      Point& crossingPt) const
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
      if(axom::utilities::isNearlyEqual(contourVal, f1) ||
         axom::utilities::isNearlyEqual(f1, f2))
      {
        crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt = p2;
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      constexpr double ptiny = axom::primal::PRIMAL_TINY;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (contourVal - f1) / df;
      for(int d = 0; d < DIM; ++d)
      {
        crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
      }
    }

    //!@brief Interpolate for the contour location crossing a parent edge.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type linear_interp(
      int edgeIdx,
      const Point cornerCoords[8],
      const double nodeValues[8],
      Point& crossingPt) const
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
      if(axom::utilities::isNearlyEqual(contourVal, f1) ||
         axom::utilities::isNearlyEqual(f1, f2))
      {
        crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt = p2;
        return;
      }

      // STEP 3: point is not at corner; interpolate its position
      constexpr double ptiny = axom::primal::PRIMAL_TINY;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (contourVal - f1) / df;
      for(int d = 0; d < DIM; ++d)
      {
        crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
      }
    }

    AXOM_HOST_DEVICE inline axom::StackArray<axom::IndexType, DIM>
    multidim_cell_index(axom::IndexType flatId) const
    {
      axom::StackArray<axom::IndexType, DIM> rval;
      for(int d = DIM - 1; d >= 0; --d)
      {
        rval[d] = flatId / bStrides[d];
        flatId -= rval[d] * bStrides[d];
      }
      return rval;
    }
  };  // ComputeContour_Util

  void computeContour() override
  {
    auto crossingsView = m_crossings.view();

    /*
      Reserve contour mesh data space so we can add data without
      reallocation.
    */
    const axom::IndexType contourNodeCount = DIM * m_contourCellCount;
    m_contourNodeCoords.resize(contourNodeCount);
    m_contourCellCorners.resize(m_contourCellCount);
    m_contourCellParents.resize(m_contourCellCount);

    auto nodeCoordsView = m_contourNodeCoords.view();
    auto cellCornersView = m_contourCellCorners.view();
    auto cellParentsView = m_contourCellParents.view();

    ComputeContour_Util ccu(m_contourVal, m_bStrides, m_fcnView, m_coordsViews);

    auto loopBody = AXOM_LAMBDA(axom::IndexType iCrossing)
    {
      const auto& crossingInfo = crossingsView[iCrossing];
      const IndexType crossingCellCount = num_contour_cells(crossingInfo.caseNum);
      SLIC_ASSERT(crossingCellCount > 0);

      // Parent cell data for interpolating new node coordinates.
      Point cornerCoords[CELL_CORNER_COUNT];
      double cornerValues[CELL_CORNER_COUNT];
      ccu.get_corner_coords_and_values(crossingInfo.parentCellNum,
                                       cornerCoords,
                                       cornerValues);

      /*
        Create the new cell and its DIM nodes.  New nodes are on
        parent cell edges where the edge intersects the isocontour.
        linear_interp for the exact coordinates.

        TODO: The varying crossingCellCount value may inhibit device
        performance.  Try grouping m_crossings items that have the
        same values for crossingCellCount.
      */
      for(int iCell = 0; iCell < crossingCellCount; ++iCell)
      {
        IndexType contourCellId = crossingInfo.firstSurfaceCellId + iCell;
        cellParentsView[contourCellId] = crossingInfo.parentCellNum;
        for(int d = 0; d < DIM; ++d)
        {
          IndexType contourNodeId = contourCellId * DIM + d;
          cellCornersView[contourCellId][d] = contourNodeId;

          const int edge = cases_table(crossingInfo.caseNum, iCell * DIM + d);
          ccu.linear_interp(edge,
                            cornerCoords,
                            cornerValues,
                            nodeCoordsView[contourNodeId]);
        }
      }
    };
    axom::for_all<ExecSpace>(0, m_crossingCount, loopBody);
  }

  // These 4 functions provide access to the look-up table
  // whether on host or device.  Is there a more elegant way
  // to put static 1D and 2D arrays on both host and device?  BTNG.

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2, int>::type
  num_contour_cells(int iCase) const
  {
#define _MC_LOOKUP_NUM_SEGMENTS
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_SEGMENTS
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return num_segments[iCase];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 2, int>::type
  cases_table(int iCase, int iEdge) const
  {
#define _MC_LOOKUP_CASES2D
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_CASES2D
    SLIC_ASSERT(iCase >= 0 && iCase < 16);
    return cases2D[iCase][iEdge];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3, int>::type
  num_contour_cells(int iCase) const
  {
#define _MC_LOOKUP_NUM_TRIANGLES
#include "marching_cubes_lookup.hpp"
#undef _MC_LOOKUP_NUM_TRIANGLES
    SLIC_ASSERT(iCase >= 0 && iCase < 256);
    return num_triangles[iCase];
  }

  template <int TDIM = DIM>
  AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3, int>::type
  cases_table(int iCase, int iEdge) const
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
  AXOM_HOST_DEVICE inline axom::StackArray<axom::IndexType, DIM>
  multidim_cell_index(axom::IndexType flatId) const
  {
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d = DIM - 1; d >= 0; --d)
    {
      rval[d] = flatId / m_bStrides[d];
      flatId -= rval[d] * m_bStrides[d];
    }
    return rval;
  }

  /*!
    @brief Output contour mesh to a mint::UnstructuredMesh object.

    mint uses host memory.  If internal memory isn't on the host,
    make a temporary copy of it on the host.
  */
  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
    const std::string& cellIdField) const override
  {
    auto internalAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto hostAllocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    if(internalAllocatorID == hostAllocatorID)
    {
      populateContourMesh(mesh,
                          cellIdField,
                          m_contourNodeCoords,
                          m_contourCellCorners,
                          m_contourCellParents);
    }
    else
    {
      axom::Array<Point, 1, MemorySpace::Dynamic> contourNodeCoords(
        m_contourNodeCoords,
        hostAllocatorID);
      axom::Array<MIdx, 1, MemorySpace::Dynamic> contourCellCorners(
        m_contourCellCorners,
        hostAllocatorID);
      axom::Array<IndexType, 1, MemorySpace::Dynamic> contourCellParents(
        m_contourCellParents,
        hostAllocatorID);

      populateContourMesh(mesh,
                          cellIdField,
                          contourNodeCoords,
                          contourCellCorners,
                          contourCellParents);
    }
  }

  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
    const std::string& cellIdField,
    axom::Array<Point, 1, MemorySpace::Dynamic> contourNodeCoords,
    axom::Array<MIdx, 1, MemorySpace::Dynamic> contourCellCorners,
    axom::Array<IndexType, 1, MemorySpace::Dynamic> contourCellParents) const
  {
    if(!cellIdField.empty() &&
       !mesh.hasField(cellIdField, axom::mint::CELL_CENTERED))
    {
      mesh.createField<axom::IndexType>(cellIdField,
                                        axom::mint::CELL_CENTERED,
                                        DIM);
    }

    const axom::IndexType addedCellCount = contourCellCorners.size();
    const axom::IndexType addedNodeCount = contourNodeCoords.size();
    if(addedCellCount != 0)
    {
      const axom::IndexType priorCellCount = mesh.getNumberOfCells();
      const axom::IndexType priorNodeCount = mesh.getNumberOfNodes();
      mesh.reserveNodes(priorNodeCount + addedNodeCount);
      mesh.reserveCells(priorCellCount + addedCellCount);

      mesh.appendNodes((double*)contourNodeCoords.data(),
                       contourNodeCoords.size());
      for(int n = 0; n < addedCellCount; ++n)
      {
        MIdx cornerIds = contourCellCorners[n];
        add_to_StackArray(cornerIds, priorNodeCount);
        mesh.appendCell(cornerIds);
      }
      axom::IndexType numComponents = -1;
      axom::IndexType* dstPtr =
        mesh.getFieldPtr<axom::IndexType>(cellIdField,
                                          axom::mint::CELL_CENTERED,
                                          numComponents);
      SLIC_ASSERT(numComponents == DIM);
      // Bump corner indices by priorCellCount to avoid indices
      // used by other parents domains.
      axom::ArrayView<axom::StackArray<axom::IndexType, DIM>> dstView(
        (axom::StackArray<axom::IndexType, DIM>*)dstPtr,
        priorCellCount + addedCellCount);
      for(axom::IndexType i = 0; i < addedCellCount; ++i)
      {
        dstView[priorCellCount + i] = multidim_cell_index(contourCellParents[i]);
      }
    }
  }

  //!@brief Compute the case index into cases2D or cases3D.
  AXOM_HOST_DEVICE inline int compute_crossing_case(const double* f) const
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
    m_contourNodeCoords.clear();
    m_contourCellCorners.clear();
    m_contourCellParents.clear();
    m_crossingCount = 0;
    m_contourCellCount = 0;
  }

  /*!
    @brief Constructor.
  */
  MarchingCubesImpl()
    : m_allocatorID(execution_space<ExecSpace>::allocatorID())
    , m_crossings(0, 0)
    , m_contourNodeCoords(0, 0)
    , m_contourCellCorners(0, 0)
    , m_contourCellParents(0, 0)
  { }

  /*!
    @brief Info for a parent cell intersecting the contour surface.
  */
  struct CrossingInfo
  {
    CrossingInfo() { }
    CrossingInfo(axom::IndexType parentCellNum_, std::uint16_t caseNum_)
      : parentCellNum(parentCellNum_)
      , caseNum(caseNum_)
      , firstSurfaceCellId(std::numeric_limits<axom::IndexType>::max())
    { }
    axom::IndexType parentCellNum;       //!< @brief Flat index of parent cell.
    std::uint16_t caseNum;               //!< @brief Index in cases2D or cases3D
    axom::IndexType firstSurfaceCellId;  //!< @brief First index for generated cells.
  };

private:
  const int m_allocatorID;  //!< @brief ExecSpace-based allocator ID for all internal data
  MIdx m_bShape;            //!< @brief Blueprint cell data shape.
  MIdx m_cShape;    //!< @brief Cell-centered array shape for ArrayViews.
  MIdx m_pShape;    //!< @brief Node-centered array shape for ArrayViews.
  MIdx m_bStrides;  //!< @brief Strides for m_bShape arrays.

  // Views of parent domain data.
  // DIM coordinate components, each on a DIM-dimensional mesh.
  using CoordViews =
    axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>;
  CoordViews m_coordsViews;
  axom::ArrayView<const double, DIM, MemorySpace> m_fcnView;
  axom::ArrayView<const int, DIM, MemorySpace> m_maskView;

  //!@brief Crossing case for each computational mesh cell.
  axom::Array<std::uint16_t, DIM, MemorySpace> m_caseIds;

  //!@brief Info on every parent cell that crosses the contour surface.
  axom::Array<CrossingInfo, 1, MemorySpace> m_crossings;

  //!@brief Number of parent cells crossing the contour surface.
  axom::IndexType m_crossingCount = 0;
  //!@brief Number of contour surface cells from crossings.
  axom::IndexType m_contourCellCount = 0;
  axom::IndexType getContourCellCount() const override
  {
    return m_contourCellCount;
  }

  //!@brief Number of corners (nodes) on each parent cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;

  //!@name Internal representation of generated contour mesh.
  //@{
  //!@brief Coordinates of generated surface mesh nodes.
  axom::Array<Point, 1, MemorySpace> m_contourNodeCoords;

  //!@brief Corners (index into m_contourNodeCoords) of generated contour cells.
  axom::Array<MIdx, 1, MemorySpace> m_contourCellCorners;

  //!@brief Flat index of computational cell crossing the contour cell.
  axom::Array<IndexType, 1, MemorySpace> m_contourCellParents;
  //@}

  double m_contourVal = 0.0;
};

}  // end namespace marching_cubes
}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
