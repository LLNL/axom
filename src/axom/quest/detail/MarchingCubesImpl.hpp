// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifndef AXOM_USE_CONDUIT
  #error "MarchingCubesImpl.hpp requires conduit"
#endif
#include "conduit_blueprint.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/quest/ArrayIndexer.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/constants.hpp"
#include "axom/mint/execution/internal/structured_exec.hpp"
#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{
namespace detail
{
namespace marching_cubes
{
/*!
  @brief Computations for MarchingCubesSingleDomain

  Spatial dimension templating is here, to keep out of higher level
  classes MarchCubes and MarchingCubesSingleDomain.

  ExecSpace is the general execution space, like axom::SEQ_EXEC and
  axom::CUDA_EXEC<256>.

  See MarchingCubesImpl for the difference between that class and
  MarchingCubesImpl.
*/
template <int DIM, typename ExecSpace, typename SequentialExecSpace>
class MarchingCubesImpl : public MarchingCubesSingleDomain::ImplBase
{
public:
  using Point = axom::primal::Point<double, DIM>;
  using MIdx = axom::StackArray<axom::IndexType, DIM>;
  using FacetIdType = int;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  using ReducePolicy = typename execution_space<ExecSpace>::reduce_policy;
  using SequentialLoopPolicy =
    typename execution_space<SequentialExecSpace>::loop_policy;
  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;

  AXOM_HOST MarchingCubesImpl(int allocatorID)
    : m_allocatorID(allocatorID)
    , m_caseIds(emptyShape(), m_allocatorID)
    , m_crossingParentIds(0, 0, m_allocatorID)
    , m_facetIncrs(0, 0, m_allocatorID)
    , m_firstFacetIds(0, 0, m_allocatorID)
  { }

  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured mesh domain
    @param topologyName Name of mesh topology (see blueprint
           mesh documentation)
    @param maskFieldName Name of integer cell mask function is in dom

    Set up views to domain data and allocate other data to work on the
    given domain.

    The above data from the domain MUST be in a memory space
    compatible with ExecSpace.
  */
  AXOM_HOST void initialize(const conduit::Node& dom,
                            const std::string& topologyName,
                            const std::string& maskFieldName) override
  {
    SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(dom.fetch_existing(
                  axom::fmt::format("topologies/{}", topologyName))) == DIM);

    // clear();

    m_mvu = axom::quest::MeshViewUtil<DIM, MemorySpace>(dom, topologyName);

    m_bShape = m_mvu.getCellShape();
    m_coordsViews = m_mvu.getConstCoordsViews(false);
    if(!maskFieldName.empty())
    {
      m_maskView = m_mvu.template getConstFieldView<int>(maskFieldName, false);
    }

    /*
      TODO: To get good cache performance, we should make m_caseIds
      row-major if fcn is that way, and vice versa.  However, Array
      only support column-major, so we're stuck with that for now.
    */
    // m_caseIds.resize(m_bShape, 0); // This unexpectedly fails.
    m_caseIds =
      axom::Array<std::uint16_t, DIM, MemorySpace>(m_bShape, m_allocatorID);
    m_caseIds.fill(0);
  }

  /*!
    @brief Set the scale field name
    @param fcnFieldName Name of nodal function is in dom
  */
  void setFunctionField(const std::string& fcnFieldName) override
  {
    m_fcnView = m_mvu.template getConstFieldView<double>(fcnFieldName, false);
  }

  void setContourValue(double contourVal) override
  {
    m_contourVal = contourVal;
  }

  /*!
    @brief Implementation of virtual markCrossings.

    Virtual methods cannot be templated, so this implementation
    delegates to an implementation templated on DIM.
  */
  void markCrossings() override { markCrossings_dim(); }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type markCrossings_dim()
  {
    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal);

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment jRange(0, m_bShape[1]);
    RAJA::RangeSegment iRange(0, m_bShape[0]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop2d_policy;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        mcu.computeCaseId(i, j);
      });
#else
    for(int j = 0; j < m_bShape[1]; ++j)
    {
      for(int i = 0; i < m_bShape[0]; ++i)
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
    MarkCrossings_Util mcu(m_caseIds, m_fcnView, m_maskView, m_contourVal);

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment kRange(0, m_bShape[2]);
    RAJA::RangeSegment jRange(0, m_bShape[1]);
    RAJA::RangeSegment iRange(0, m_bShape[0]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop3d_policy;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange, kRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        mcu.computeCaseId(i, j, k);
      });
#else
    for(int k = 0; k < m_bShape[2]; ++k)
    {
      for(int j = 0; j < m_bShape[1]; ++j)
      {
        for(int i = 0; i < m_bShape[0]; ++i)
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
      const bool useZone = maskView.empty() || bool(maskView(i, j));
      if(useZone)
      {
        // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {fcnView(i    , j    ),
             fcnView(i + 1, j    ),
             fcnView(i + 1, j + 1),
             fcnView(i    , j + 1)};
        // clang-format on
        caseIdsView(i, j) = computeCrossingCase(nodalValues);
      }
    }

    //!@brief Populate m_caseIds with crossing indices.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE inline typename std::enable_if<TDIM == 3>::type
    computeCaseId(axom::IndexType i, axom::IndexType j, axom::IndexType k) const
    {
      const bool useZone = maskView.empty() || bool(maskView(i, j, k));
      if(useZone)
      {
        // clang-format off
          double nodalValues[CELL_CORNER_COUNT] =
            {fcnView(i + 1, j    , k    ),
             fcnView(i + 1, j + 1, k    ),
             fcnView(i    , j + 1, k    ),
             fcnView(i    , j    , k    ),
             fcnView(i + 1, j    , k + 1),
             fcnView(i + 1, j + 1, k + 1),
             fcnView(i    , j + 1, k + 1),
             fcnView(i    , j    , k + 1)};
        // clang-format on
        caseIdsView(i, j, k) = computeCrossingCase(nodalValues);
      }
    }
  };  // MarkCrossings_Util

  void scanCrossings() override
  {
    constexpr MarchingCubesDataParallelism autoPolicy =
      std::is_same<ExecSpace, axom::SEQ_EXEC>::value
      ? MarchingCubesDataParallelism::hybridParallel
      :
#ifdef AXOM_USE_OPENMP
      std::is_same<ExecSpace, axom::OMP_EXEC>::value
        ? MarchingCubesDataParallelism::hybridParallel
        :
#endif
        MarchingCubesDataParallelism::fullParallel;

    if(m_dataParallelism ==
         axom::quest::MarchingCubesDataParallelism::hybridParallel ||
       (m_dataParallelism == axom::quest::MarchingCubesDataParallelism::byPolicy &&
        autoPolicy == MarchingCubesDataParallelism::hybridParallel))
    {
      scanCrossings_hybridParallel();
    }
    else
    {
      scanCrossings_fullParallel();
    }
  }

  void scanCrossings_fullParallel()
  {
#if defined(AXOM_USE_RAJA)
  #ifdef __INTEL_LLVM_COMPILER
    // Intel oneAPI compiler segfaults with OpenMP RAJA scan
    using ScanPolicy =
      typename axom::execution_space<axom::SEQ_EXEC>::loop_policy;
  #else
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;
  #endif
#endif
    const axom::IndexType parentCellCount = m_caseIds.size();
    auto caseIdsView = m_caseIds.view();

    //
    // Initialize crossingFlags to 0 or 1, prefix-sum it, and count the
    // crossings.
    //
    // All 1D array data alligns with m_caseId, which is column-major
    // (the default for Array class), leading to *column-major* parent
    // cell ids, regardless of the ordering of the input mesh data.
    //

    axom::Array<axom::IndexType, 1, MemorySpace> crossingFlags(parentCellCount);
    auto crossingFlagsView = crossingFlags.view();
    axom::for_all<ExecSpace>(
      0,
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType parentCellId) {
        auto numContourCells =
          num_contour_cells(caseIdsView.flatIndex(parentCellId));
        crossingFlagsView[parentCellId] = bool(numContourCells);
      });

    axom::Array<axom::IndexType, 1, MemorySpace> scannedFlags(1 + parentCellCount);
    auto scannedFlagsView = scannedFlags.view();
    RAJA::inclusive_scan<ScanPolicy>(
      RAJA::make_span(crossingFlags.data(), parentCellCount),
      RAJA::make_span(scannedFlags.data() + 1, parentCellCount),
      RAJA::operators::plus<axom::IndexType> {});

    axom::copy(&m_crossingCount,
               scannedFlags.data() + scannedFlags.size() - 1,
               sizeof(axom::IndexType));

    //
    // Generate crossing-cells index list and corresponding facet counts.
    //

    m_crossingParentIds.resize(m_crossingCount);
    m_facetIncrs.resize(m_crossingCount);
    m_firstFacetIds.resize(1 + m_crossingCount);

    auto crossingParentIdsView = m_crossingParentIds.view();
    auto facetIncrsView = m_facetIncrs.view();

    auto loopBody = AXOM_LAMBDA(axom::IndexType parentCellId)
    {
      if(scannedFlagsView[parentCellId] != scannedFlagsView[1 + parentCellId])
      {
        auto crossingId = scannedFlagsView[parentCellId];
        auto facetIncr = num_contour_cells(caseIdsView.flatIndex(parentCellId));
        crossingParentIdsView[crossingId] = parentCellId;
        facetIncrsView[crossingId] = facetIncr;
      }
    };
    axom::for_all<ExecSpace>(0, parentCellCount, loopBody);

    //
    // Prefix-sum the facets counts to get first facet id for each crossing
    // and the total number of facets.
    //

    RAJA::inclusive_scan<ScanPolicy>(
      RAJA::make_span(m_facetIncrs.data(), m_crossingCount),
      RAJA::make_span(m_firstFacetIds.data() + 1, m_crossingCount),
      RAJA::operators::plus<axom::IndexType> {});

    axom::copy(&m_facetCount,
               m_firstFacetIds.data() + m_firstFacetIds.size() - 1,
               sizeof(axom::IndexType));
    m_firstFacetIds.resize(m_crossingCount);
  }

  void scanCrossings_hybridParallel()
  {
    //
    // Compute number of crossings in m_caseIds
    //
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

    //
    // Allocate space for crossing info
    //
    m_crossingParentIds.resize(m_crossingCount);
    m_facetIncrs.resize(m_crossingCount);
    m_firstFacetIds.resize(1 + m_crossingCount);

    auto crossingParentIdsView = m_crossingParentIds.view();
    auto facetIncrsView = m_facetIncrs.view();

    axom::IndexType* crossingId = axom::allocate<axom::IndexType>(
      1,
      axom::detail::getAllocatorID<MemorySpace>());

    auto loopBody = AXOM_LAMBDA(axom::IndexType n)
    {
      auto caseId = caseIdsView.flatIndex(n);
      auto ccc = num_contour_cells(caseId);
      if(ccc != 0)
      {
        facetIncrsView[*crossingId] = ccc;
        crossingParentIdsView[*crossingId] = n;
        ++(*crossingId);
      }
    };

#if defined(AXOM_USE_RAJA)
    /*
      loopBody isn't data-parallel and shouldn't be parallelized.
      This contrived RAJA::forall forces it to run sequentially.
    */
    RAJA::forall<SequentialLoopPolicy>(
      RAJA::RangeSegment(0, 1),
      [=] AXOM_HOST_DEVICE(int /* i */) {
        *crossingId = 0;
        for(axom::IndexType n = 0; n < parentCellCount; ++n)
        {
          loopBody(n);
        }
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

    // axom::Array<axom::IndexType, 1, MemorySpace> prefixSum(m_crossingCount, m_crossingCount);
    const auto firstFacetIdsView = m_firstFacetIds.view();

#if defined(AXOM_USE_RAJA)
    // Intel oneAPI compiler segfaults with OpenMP RAJA scan
  #ifdef __INTEL_LLVM_COMPILER
    using ScanPolicy =
      typename axom::execution_space<axom::SEQ_EXEC>::loop_policy;
  #else
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;
  #endif
    RAJA::inclusive_scan<ScanPolicy>(
      RAJA::make_span(facetIncrsView.data(), m_crossingCount),
      RAJA::make_span(firstFacetIdsView.data() + 1, m_crossingCount),
      RAJA::operators::plus<axom::IndexType> {});
#else
    if(m_crossingCount > 0)
    {
      firstFacetIdsView[0] = 0;
      for(axom::IndexType i = 1; i < 1 + m_crossingCount; ++i)
      {
        firstFacetIdsView[i] = firstFacetIdsView[i - 1] + facetIncrsView[i - 1];
      }
    }
#endif
    axom::copy(&m_facetCount,
               m_firstFacetIds.data() + m_firstFacetIds.size() - 1,
               sizeof(axom::IndexType));
    m_firstFacetIds.resize(m_crossingCount);
  }

  void computeContour() override
  {
    const auto facetIncrsView = m_facetIncrs.view();
    const auto firstFacetIdsView = m_firstFacetIds.view();
    const auto crossingParentIdsView = m_crossingParentIds.view();
    const auto caseIdsView = m_caseIds.view();

    // Internal contour mesh data to populate
    axom::ArrayView<axom::IndexType, 2> facetNodeIdsView = m_facetNodeIds;
    axom::ArrayView<double, 2> facetNodeCoordsView = m_facetNodeCoords;
    axom::ArrayView<axom::IndexType> facetParentIdsView = m_facetParentIds;
    const axom::IndexType facetIndexOffset = m_facetIndexOffset;

    ComputeContour_Util ccu(m_contourVal,
                            m_caseIds.strides(),
                            m_fcnView,
                            m_coordsViews);

    auto gen_for_parent_cell = AXOM_LAMBDA(axom::IndexType crossingId)
    {
      auto parentCellId = crossingParentIdsView[crossingId];
      auto caseId = caseIdsView.flatIndex(parentCellId);
      Point cornerCoords[CELL_CORNER_COUNT];
      double cornerValues[CELL_CORNER_COUNT];
      ccu.get_corner_coords_and_values(parentCellId, cornerCoords, cornerValues);

      auto additionalFacets = facetIncrsView[crossingId];
      auto firstFacetId = facetIndexOffset + firstFacetIdsView[crossingId];

      for(axom::IndexType fId = 0; fId < additionalFacets; ++fId)
      {
        axom::IndexType newFacetId = firstFacetId + fId;
        axom::IndexType firstCornerId = newFacetId * DIM;

        facetParentIdsView[newFacetId] = parentCellId;

        for(axom::IndexType d = 0; d < DIM; ++d)
        {
          axom::IndexType newCornerId = firstCornerId + d;
          facetNodeIdsView[newFacetId][d] = newCornerId;

          int edge = cases_table(caseId, fId * DIM + d);
          ccu.linear_interp(edge,
                            cornerCoords,
                            cornerValues,
                            &facetNodeCoordsView(newCornerId, 0));
        }
      }
    };

    axom::for_all<ExecSpace>(0, m_crossingCount, gen_for_parent_cell);
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
    axom::ArrayIndexer<axom::IndexType, DIM> indexer;
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
      , indexer(bStrides_)
      , fcnView(fcnView_)
      , coordsViews(coordsViews_)
    { }

    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type
    get_corner_coords_and_values(IndexType parentCellId,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];

      const auto c = indexer.toMultiIndex(parentCellId);
      const auto& i = c[0];
      const auto& j = c[1];

      // clang-format off
      cornerCoords[0] = { x(i  , j  ), y(i  , j  ) };
      cornerCoords[1] = { x(i+1, j  ), y(i+1, j  ) };
      cornerCoords[2] = { x(i+1, j+1), y(i+1, j+1) };
      cornerCoords[3] = { x(i  , j+1), y(i  , j+1) };

      cornerValues[0] = fcnView(i  , j  );
      cornerValues[1] = fcnView(i+1, j  );
      cornerValues[2] = fcnView(i+1, j+1);
      cornerValues[3] = fcnView(i  , j+1);
      // clang-format on
    }
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3>::type
    get_corner_coords_and_values(IndexType parentCellId,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];
      const auto& z = coordsViews[2];

      const auto c = indexer.toMultiIndex(parentCellId);
      const auto& i = c[0];
      const auto& j = c[1];
      const auto& k = c[2];

      // clang-format off
      cornerCoords[0] = { x(i+1, j  , k  ), y(i+1, j  , k  ), z(i+1, j  , k  ) };
      cornerCoords[1] = { x(i+1, j+1, k  ), y(i+1, j+1, k  ), z(i+1, j+1, k  ) };
      cornerCoords[2] = { x(i  , j+1, k  ), y(i  , j+1, k  ), z(i  , j+1, k  ) };
      cornerCoords[3] = { x(i  , j  , k  ), y(i  , j  , k  ), z(i  , j  , k  ) };
      cornerCoords[4] = { x(i+1, j  , k+1), y(i+1, j  , k+1), z(i+1, j  , k+1) };
      cornerCoords[5] = { x(i+1, j+1, k+1), y(i+1, j+1, k+1), z(i+1, j+1, k+1) };
      cornerCoords[6] = { x(i  , j+1, k+1), y(i  , j+1, k+1), z(i  , j+1, k+1) };
      cornerCoords[7] = { x(i  , j  , k+1), y(i  , j  , k+1), z(i  , j  , k+1) };

      cornerValues[0] = fcnView(i+1, j  , k  );
      cornerValues[1] = fcnView(i+1, j+1, k  );
      cornerValues[2] = fcnView(i  , j+1, k  );
      cornerValues[3] = fcnView(i  , j  , k  );
      cornerValues[4] = fcnView(i+1, j  , k+1);
      cornerValues[5] = fcnView(i+1, j+1, k+1);
      cornerValues[6] = fcnView(i  , j+1, k+1);
      cornerValues[7] = fcnView(i  , j  , k+1);
      // clang-format on
    }

    //!@brief Interpolate for the contour location crossing a parent edge.
    template <int TDIM = DIM>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2>::type linear_interp(
      int edgeIdx,
      const Point cornerCoords[4],
      const double nodeValues[4],
      double* /* Point& */ crossingPt) const
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
        crossingPt[0] = p1[0];
        crossingPt[1] = p1[1];  // crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt[0] = p2[0];
        crossingPt[1] = p2[1];  // crossingPt = p2;
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
      double* /* Point& */ crossingPt) const
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
        crossingPt[0] = p1[0];
        crossingPt[1] = p1[1];
        crossingPt[2] = p1[2];  // crossingPt = p1;
        return;
      }

      if(axom::utilities::isNearlyEqual(contourVal, f2))
      {
        crossingPt[0] = p2[0];
        crossingPt[1] = p2[1];
        crossingPt[2] = p2[2];  // crossingPt = p2;
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
  };  // ComputeContour_Util

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

  /*!
    @brief Constructor.
  */
  MarchingCubesImpl() { }

private:
  int m_allocatorID;

  axom::quest::MeshViewUtil<DIM, MemorySpace> m_mvu;
  MIdx m_bShape;  //!< @brief Blueprint cell data shape.

  // Views of parent domain data.
  // DIM coordinate components, each on a DIM-dimensional mesh.
  using CoordViews =
    axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>;
  CoordViews m_coordsViews;
  axom::ArrayView<const double, DIM, MemorySpace> m_fcnView;
  axom::ArrayView<const int, DIM, MemorySpace> m_maskView;

  //!@brief Crossing case for each computational mesh cell.
  axom::Array<std::uint16_t, DIM, MemorySpace> m_caseIds;

  //!@brief Number of parent cells crossing the contour surface.
  axom::IndexType m_crossingCount = 0;

  //!@brief Number of contour surface cells from all crossings.
  axom::IndexType m_facetCount = 0;
  axom::IndexType getContourCellCount() const override { return m_facetCount; }

  //!@brief Parent cell id (flat index into m_caseIds) for each crossing.
  axom::Array<axom::IndexType, 1, MemorySpace> m_crossingParentIds;

  //!@brief Number of surface mesh facets added by each crossing.
  axom::Array<FacetIdType, 1, MemorySpace> m_facetIncrs;

  //!@brief First index of facets for each crossing.
  axom::Array<axom::IndexType, 1, MemorySpace> m_firstFacetIds;

  //!@brief Number of corners (nodes) on each parent cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;

  double m_contourVal = 0.0;

  axom::StackArray<axom::IndexType, DIM> emptyShape()
  {
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      rval[d] = 0;
    }
    return rval;
  }
};

/*!
  @brief Allocate a MarchingCubesImpl object, template-specialized
  for caller-specified runtime policy and physical dimension.
*/
static std::unique_ptr<axom::quest::MarchingCubesSingleDomain::ImplBase>
newMarchingCubesImpl(MarchingCubes::RuntimePolicy runtimePolicy,
                     int allocatorID,
                     int dim)
{
  using ImplBase = axom::quest::MarchingCubesSingleDomain::ImplBase;

  SLIC_ASSERT(dim >= 2 && dim <= 3);
  std::unique_ptr<ImplBase> impl;
  if(runtimePolicy == MarchingCubes::RuntimePolicy::seq)
  {
    impl = dim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::SEQ_EXEC, axom::SEQ_EXEC>(allocatorID))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::SEQ_EXEC, axom::SEQ_EXEC>(allocatorID));
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(runtimePolicy == MarchingCubes::RuntimePolicy::omp)
  {
    impl = dim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::OMP_EXEC, axom::SEQ_EXEC>(allocatorID))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::OMP_EXEC, axom::SEQ_EXEC>(allocatorID));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(runtimePolicy == MarchingCubes::RuntimePolicy::cuda)
  {
    impl = dim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>(
            allocatorID))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::CUDA_EXEC<256>, axom::CUDA_EXEC<1>>(
            allocatorID));
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(runtimePolicy == MarchingCubes::RuntimePolicy::hip)
  {
    impl = dim == 2
      ? std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<2, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>(
            allocatorID))
      : std::unique_ptr<ImplBase>(
          new MarchingCubesImpl<3, axom::HIP_EXEC<256>, axom::HIP_EXEC<1>>(
            allocatorID));
  }
#endif
  else
  {
    SLIC_ERROR(axom::fmt::format(
      "MarchingCubesSingleDomain has no implementation for runtime policy {}",
      runtimePolicy));
  }
  return impl;
}

}  // end namespace marching_cubes
}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
