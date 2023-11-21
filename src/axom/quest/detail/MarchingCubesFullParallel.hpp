// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT
  #include "conduit_blueprint.hpp"

  #include "axom/core/execution/execution_space.hpp"
  #include "axom/quest/ArrayIndexer.hpp"
  #include "axom/quest/detail/marching_cubes_lookup.hpp"
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
  MarchingCubesFullParallel.
*/
template <int DIM, typename ExecSpace>
class MarchingCubesFullParallel : public MarchingCubesSingleDomain::ImplBase
{
public:
  using Point = axom::primal::Point<double, DIM>;
  using MIdx = axom::StackArray<axom::IndexType, DIM>;
  using LoopPolicy = typename execution_space<ExecSpace>::loop_policy;
  using ReducePolicy = typename execution_space<ExecSpace>::reduce_policy;
  static constexpr auto MemorySpace = execution_space<ExecSpace>::memory_space;
  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured mesh domain
    @param topologyName Name of mesh topology (see blueprint
           mesh documentation)
    @param fcnFieldName Name of nodal function is in dom
    @param maskFieldName Name of integer cell mask function is in dom

    Set up views to domain data and allocate other data to work on the
    given domain.

    The above data from the domain MUST be in a memory space
    compatible with ExecSpace.
  */
  AXOM_HOST void initialize(const conduit::Node& dom,
                            const std::string& topologyName,
                            const std::string& fcnFieldName,
                            const std::string& maskFieldName = {}) override
  {
    SLIC_ASSERT(conduit::blueprint::mesh::topology::dims(dom.fetch_existing(
                  axom::fmt::format("topologies/{}", topologyName))) == DIM);

    clear();

    axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(dom, topologyName);

    m_bShape = mvu.getCellShape();
    m_coordsViews = mvu.getConstCoordsViews(false);
    m_fcnView = mvu.template getConstFieldView<double>(fcnFieldName, false);
    if(!maskFieldName.empty())
    {
      m_maskView = mvu.template getConstFieldView<int>(maskFieldName, false);
    }

    /*
      TODO: To get good cache performance, we should make m_caseIds
      row-major if fcn is that way, and vice versa.  However, Array
      only support column-major, so we're stuck with that for now.
    */
    m_caseIds = axom::Array<std::uint16_t, DIM, MemorySpace>(m_bShape);
    m_caseIds.fill(0);
  }

  //!@brief Set a value to find the contour for.
  void setContourValue(double contourVal) override
  {
    m_contourVal = contourVal;
  }

  void computeContourMesh() override
  {
    markCrossings();
    scanCrossings();
    computeContour();
  }

  /*!
    @brief Implementation of virtual markCrossings.

    Virtual methods cannot be templated, so this implementation
    delegates to a name templated on DIM.
  */
  void markCrossings() { markCrossings_dim(); }

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
    @brief Implementation used by MarchingCubesFullParallel::markCrossings_dim()
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

  /*!
    @brief Populate the 1D m_contourNodeCoords, m_contourCellCorners
    and m_contourCellParents arrays that defines the unstructured
    contour mesh.
  */
  void scanCrossings()
  {
  #ifdef __INTEL_LLVM_COMPILER
    // Intel oneAPI compiler segfaults with OpenMP RAJA scan
    using ScanPolicy =
      typename axom::execution_space<axom::SEQ_EXEC>::loop_policy;
  #else
    using ScanPolicy = typename axom::execution_space<ExecSpace>::loop_policy;
  #endif

    const axom::IndexType parentCellCount = m_caseIds.size();
    auto caseIdsView = m_caseIds.view();

    // Compute number of surface facets added by each parent cell.
    m_facetIncrs.resize(parentCellCount);
    const axom::ArrayView<int, 1, MemorySpace> facetIncrsView = m_facetIncrs.view();

  #if defined(AXOM_USE_RAJA)
    axom::for_all<ExecSpace>(
      0,
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType parentCellId) {
        facetIncrsView.flatIndex(parentCellId) =
          num_contour_cells(caseIdsView.flatIndex(parentCellId));
      });
  #else
    for(axom::IndexType pcId = 0; pcId < parentCellCount; ++pcId)
    {
      facetIncrsView.flatIndex(pcId) =
        num_contour_cells(caseIdsView.flatIndex(pcId));
    }
  #endif

    // Compute index of first facet added by each parent cell
    // (whether the cell generates any facet!).
    m_firstFacetIds.resize(parentCellCount);
    const axom::ArrayView<axom::IndexType, 1, MemorySpace> firstFacetIdsView =
      m_firstFacetIds.view();
  #if defined(AXOM_USE_RAJA)
    RAJA::exclusive_scan<ScanPolicy>(
      RAJA::make_span(facetIncrsView.data(), parentCellCount),
      RAJA::make_span(firstFacetIdsView.data(), parentCellCount),
      RAJA::operators::plus<axom::IndexType> {});

      // m_facetIncrs and m_firstFacetIds, combined with m_caseIds,
      // are all we need to compute the surface mesh.
  #else
    firstFacetIdsView[0] = 0;
    for(axom::IndexType pcId = 1; pcId < parentCellCount; ++pcId)
    {
      firstFacetIdsView[pcId] =
        firstFacetIdsView[pcId - 1] + facetIncrsView[pcId - 1];
    }
  #endif

    // Use last facet info to compute number of facets in domain.
    // In case data is on device, copy to host before computing.
    axom::IndexType firstFacetIds_back = 0;
    axom::IndexType facetIncrs_back = 0;
    axom::copy(&firstFacetIds_back,
               m_firstFacetIds.data() + m_firstFacetIds.size() - 1,
               sizeof(firstFacetIds_back));
    axom::copy(&facetIncrs_back,
               m_facetIncrs.data() + m_facetIncrs.size() - 1,
               sizeof(facetIncrs_back));
    m_facetCount = firstFacetIds_back + facetIncrs_back;

    // Allocate space for surface mesh.
    const axom::IndexType cornersCount = DIM * m_facetCount;
    m_contourCellParents.resize(m_facetCount);
    m_contourCellCorners.resize(m_facetCount);
    m_contourNodeCoords.resize(cornersCount);
  }

  void computeContour()
  {
    //
    // Fill in surface mesh data.
    //

    const axom::IndexType parentCellCount = m_caseIds.size();
    const auto facetIncrsView = m_facetIncrs.view();
    const auto firstFacetIdsView = m_firstFacetIds.view();
    const auto caseIdsView = m_caseIds.view();

    // sortedIndices are parent cell indices, sorted by number
    // of facets in them.
    axom::Array<axom::IndexType, 1, MemorySpace> sortedFacetIncrs(m_facetIncrs);
    axom::Array<axom::IndexType, 1, MemorySpace> sortedIndices(parentCellCount);
    auto sortedIndicesView = sortedIndices.view();
    axom::for_all<ExecSpace>(0, parentCellCount,
                             AXOM_LAMBDA(axom::IndexType pcId) {
                               sortedIndicesView[pcId] = pcId;
                             });
    RAJA::stable_sort_pairs<LoopPolicy>(RAJA::make_span(sortedFacetIncrs.data(), parentCellCount),
                                        RAJA::make_span(sortedIndices.data(), parentCellCount),
                                        RAJA::operators::greater<axom::IndexType>{});

    auto contourCellParentsView = m_contourCellParents.view();
    auto contourCellCornersView = m_contourCellCorners.view();
    auto contourNodeCoordsView = m_contourNodeCoords.view();

    ComputeContour_Util ccu(m_contourVal,
                            m_caseIds.strides(),
                            m_fcnView,
                            m_coordsViews);
    auto gen_for_parent_cell = AXOM_LAMBDA(axom::IndexType loopIndex)
    {
      axom::IndexType parentCellId = sortedIndicesView[loopIndex];
      Point cornerCoords[CELL_CORNER_COUNT];
      double cornerValues[CELL_CORNER_COUNT];
      ccu.get_corner_coords_and_values(parentCellId, cornerCoords, cornerValues);

      auto additionalFacets = facetIncrsView[parentCellId];
      auto firstFacetId = firstFacetIdsView[parentCellId];

      auto caseId = caseIdsView.flatIndex(parentCellId);

      for(axom::IndexType fId = 0; fId < additionalFacets; ++fId)
      {
        axom::IndexType newFacetId = firstFacetId + fId;
        axom::IndexType firstCornerId = newFacetId * DIM;

        contourCellParentsView[newFacetId] = parentCellId;

        for(axom::IndexType d = 0; d < DIM; ++d)
        {
          axom::IndexType newCornerId = firstCornerId + d;
          contourCellCornersView[newFacetId][d] = newCornerId;

          int edge = cases_table(caseId, fId * DIM + d);
          ccu.linear_interp(edge,
                            cornerCoords,
                            cornerValues,
                            contourNodeCoordsView[newCornerId]);
        }
      }
    };

  #if defined(AXOM_USE_RAJA)
    axom::for_all<ExecSpace>(0, parentCellCount, gen_for_parent_cell);
  #else
    for(axom::IndexType pcId = 1; pcId < parentCellCount; ++pcId)
    {
      gen_for_parent_cell(pcId);
    }
  #endif
  }

  /*!
    @brief Implementation used by MarchingCubesFullParallel::computeContour().
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
    get_corner_coords_and_values(IndexType cellNum,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];

      const auto c = indexer.toMultiIndex(cellNum);
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
    get_corner_coords_and_values(IndexType cellNum,
                                 Point cornerCoords[],
                                 double cornerValues[]) const
    {
      const auto& x = coordsViews[0];
      const auto& y = coordsViews[1];
      const auto& z = coordsViews[2];

      const auto c = indexer.toMultiIndex(cellNum);
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

  /*!
    @brief Output contour mesh to a mint::UnstructuredMesh object.
  */
  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
    const std::string& cellIdField) const override
  {
    auto internalAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto hostAllocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

    /*
      mint uses host memory.  If internal memory is on the host, use
      it.  Otherwise, make a temporary copy of it on the host.
    */
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

  //!@brief Output contour mesh to a mint::UnstructuredMesh object.
  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& mesh,
    const std::string& cellIdField,
    const axom::Array<Point, 1, MemorySpace::Dynamic> contourNodeCoords,
    const axom::Array<MIdx, 1, MemorySpace::Dynamic> contourCellCorners,
    const axom::Array<IndexType, 1, MemorySpace::Dynamic> contourCellParents) const
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
        // Bump corner indices by priorNodeCount to avoid indices
        // used by other parents domains.
        for(int d = 0; d < DIM; ++d)
        {
          cornerIds[d] += priorNodeCount;
        }
        mesh.appendCell(cornerIds);
      }
      axom::IndexType numComponents = -1;
      axom::IndexType* cellIdPtr =
        mesh.getFieldPtr<axom::IndexType>(cellIdField,
                                          axom::mint::CELL_CENTERED,
                                          numComponents);
      SLIC_ASSERT(numComponents == DIM);
      axom::ArrayView<axom::StackArray<axom::IndexType, DIM>> cellIdView(
        (axom::StackArray<axom::IndexType, DIM>*)cellIdPtr,
        priorCellCount + addedCellCount);
      axom::ArrayIndexer<axom::IndexType, DIM> si(m_caseIds.shape(), 'c');
      for(axom::IndexType i = 0; i < addedCellCount; ++i)
      {
        cellIdView[priorCellCount + i] = si.toMultiIndex(contourCellParents[i]);
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
    m_facetCount = 0;
  }

  /*!
    @brief Constructor.
  */
  MarchingCubesFullParallel()
    : m_contourNodeCoords(0, 0)
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
    axom::IndexType parentCellNum;  //!< @brief Flat index of parent cell in m_caseIds.
    std::uint16_t caseNum;          //!< @brief Index in cases2D or cases3D
    axom::IndexType firstSurfaceCellId;  //!< @brief First index for generated cells.
  };

private:
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

  //!@brief Number of surface mesh facets added by computational mesh cells.
  axom::Array<int, 1, MemorySpace> m_facetIncrs;

  //!@brief First index of facets in computational mesh cells.
  axom::Array<axom::IndexType, 1, MemorySpace> m_firstFacetIds;

  //!@brief Number of contour surface cells from crossings.
  axom::IndexType m_facetCount = 0;
  axom::IndexType getContourCellCount() const override { return m_facetCount; }

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

static std::unique_ptr<axom::quest::MarchingCubesSingleDomain::ImplBase>
newMarchingCubesFullParallel(MarchingCubesRuntimePolicy runtimePolicy, int dim)
{
  using ImplBase = axom::quest::MarchingCubesSingleDomain::ImplBase;
  using RuntimePolicy = axom::quest::MarchingCubesRuntimePolicy;

  SLIC_ASSERT(dim >= 2 && dim <= 3);
  std::unique_ptr<ImplBase> impl;
  if(runtimePolicy == RuntimePolicy::seq)
  {
    impl = dim == 2 ? std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<2, axom::SEQ_EXEC>)
                    : std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<3, axom::SEQ_EXEC>);
  }
  #ifdef _AXOM_MARCHINGCUBES_USE_OPENMP
  else if(runtimePolicy == RuntimePolicy::omp)
  {
    impl = dim == 2 ? std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<2, axom::OMP_EXEC>)
                    : std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<3, axom::OMP_EXEC>);
  }
  #endif
  #ifdef _AXOM_MARCHINGCUBES_USE_CUDA
  else if(runtimePolicy == RuntimePolicy::cuda)
  {
    impl = dim == 2 ? std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<2, axom::CUDA_EXEC<256>>)
                    : std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<3, axom::CUDA_EXEC<256>>);
  }
  #endif
  #ifdef _AXOM_MARCHINGCUBES_USE_HIP
  else if(runtimePolicy == RuntimePolicy::hip)
  {
    impl = dim == 2 ? std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<2, axom::HIP_EXEC<256>>)
                    : std::unique_ptr<ImplBase>(
                        new MarchingCubesFullParallel<3, axom::HIP_EXEC<256>>);
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
#endif  // AXOM_USE_CONDUIT

}  // end namespace marching_cubes
}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
