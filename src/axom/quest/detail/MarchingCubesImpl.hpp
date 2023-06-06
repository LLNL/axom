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
*/
template <int DIM, typename ExecSpace>
class MarchingCubesImpl : public MarchingCubesSingleDomain::ImplBase
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

    m_caseIds = axom::Array<std::uint16_t, DIM>(
      m_cShape,
      axom::execution_space<ExecSpace>::allocatorID());
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
    auto caseIdsView = m_caseIds.view();

    auto loopBody = AXOM_LAMBDA(axom::IndexType i, axom::IndexType j)
    {
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
        caseIdsView(j, i) = crossingCase;
      }
    };

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment jRange(0, m_cShape[0]);
    RAJA::RangeSegment iRange(0, m_cShape[1]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop2d_policy;
    RAJA::kernel<EXEC_POL>(RAJA::make_tuple(iRange, jRange), loopBody);
#else
    for(int j = 0; j < m_cShape[0]; ++j)
    {
      for(int i = 0; i < m_cShape[1]; ++i)
      {
        loopBody(i, j);
      }
    }
#endif
  }
  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type markCrossings_dim()
  {
    auto caseIdsView = m_caseIds.view();

    auto loopBody =
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k)
    {
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
        caseIdsView(k, j, i) = crossingCase;
      }
    };

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment kRange(0, m_cShape[0]);
    RAJA::RangeSegment jRange(0, m_cShape[1]);
    RAJA::RangeSegment iRange(0, m_cShape[2]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<ExecSpace>::loop3d_policy;
    RAJA::kernel<EXEC_POL>(RAJA::make_tuple(iRange, jRange, kRange), loopBody);
#else
    for(int k = 0; k < m_cShape[0]; ++k)
    {
      for(int j = 0; j < m_cShape[1]; ++j)
      {
        for(int i = 0; i < m_cShape[2]; ++i)
        {
          loopBody(i, j, k);
        }
      }
    }
#endif
  }

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

#if defined(AXOM_USE_LAMBDA)
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
    axom::ArrayView<CrossingInfo> crossingsView = m_crossings.view();

    axom::Array<int> addCells(m_crossingCount,
                              m_crossingCount,
                              m_crossings.getAllocatorID());
    axom::ArrayView<int> addCellsView = addCells.view();

    axom::IndexType* crossingId = axom::allocate<axom::IndexType>(
      1,
      axom::execution_space<ExecSpace>::allocatorID());
    *crossingId = 0;
    axom::for_all<ExecSpace>(
      0,
      parentCellCount,
      AXOM_LAMBDA(axom::IndexType n) {
        auto caseId = caseIdsView.flatIndex(n);
        auto ccc = num_contour_cells(caseId);
        if(ccc != 0)
        {
          addCellsView[*crossingId] = ccc;
          m_crossings[*crossingId].caseNum = caseId;
          m_crossings[*crossingId].parentCellNum = n;
          ++(*crossingId);
        }
      });
    SLIC_ASSERT(*crossingId == m_crossingCount);
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

    m_contourCellCount = m_crossings.empty()
      ? 0
      : m_crossings.back().firstSurfaceCellId +
        num_contour_cells(m_crossings.back().caseNum);
  }

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

    axom::for_all<ExecSpace>(
      0,
      m_crossingCount,
      AXOM_LAMBDA(axom::IndexType iCrossing) {
        const auto& crossingInfo = crossingsView[iCrossing];
        const IndexType crossingCellCount =
          num_contour_cells(crossingInfo.caseNum);
        SLIC_ASSERT(crossingCellCount > 0);

        // Parent cell data for interpolating new node coordinates.
        Point cornerCoords[CELL_CORNER_COUNT];
        double cornerValues[CELL_CORNER_COUNT];
        get_corner_coords_and_values(crossingInfo.parentCellNum,
                                     cornerCoords,
                                     cornerValues);

        /*
          Create the new cell and its DIM nodes.  New node are on
          parent cell edges where the edge intersects the isocontour.
          linear_interp for the exact coordinates.

          TODO: The varying crossingCellCount value may inhibit device
          performance.  Try grouping m_crossings items that have the
          same values for crossingCellCount.
        */
        for(int iCell = 0; iCell < crossingCellCount; ++iCell)
        {
          IndexType contourCellId = crossingInfo.firstSurfaceCellId + iCell;
          m_contourCellParents[contourCellId] = crossingInfo.parentCellNum;
          for(int d = 0; d < DIM; ++d)
          {
            IndexType contourNodeId = contourCellId * DIM + d;
            m_contourCellCorners[contourCellId][d] = contourNodeId;

            const int edge = cases_table(crossingInfo.caseNum, iCell * DIM + d);
            linear_interp(edge,
                          cornerCoords,
                          cornerValues,
                          m_contourNodeCoords[contourNodeId]);
          }
        }
      });
  }

  // These 4 functions provides access to the look-up table
  // whether on host or device.  Is there a more elegant way
  // to put static 1D and 2D arrays on both host and device?  BTNG.
  template <int TDIM = DIM>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, int>::type num_contour_cells(
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
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, int>::type num_contour_cells(
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
    axom::StackArray<axom::IndexType, DIM> rval;
    for(int d = DIM - 1; d >= 0; --d)
    {
      rval[d] = flatId / m_bStrides[d];
      flatId -= rval[d] * m_bStrides[d];
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
    cornerCoords[0] = { x(j    , i    ), y(j    , i    ) };
    cornerCoords[1] = { x(j    , i + 1), y(j    , i + 1) };
    cornerCoords[2] = { x(j + 1, i + 1), y(j + 1, i + 1) };
    cornerCoords[3] = { x(j + 1, i    ), y(j + 1, i    ) };

    cornerValues[0] = m_fcnView(j    , i    );
    cornerValues[1] = m_fcnView(j    , i + 1);
    cornerValues[2] = m_fcnView(j + 1, i + 1);
    cornerValues[3] = m_fcnView(j + 1, i    );
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
    cornerCoords[0] = { x(k  , j  , i+1), y(k  , j  , i+1), z(k  , j  , i+1) };
    cornerCoords[1] = { x(k  , j+1, i+1), y(k  , j+1, i+1), z(k  , j+1, i+1) };
    cornerCoords[2] = { x(k  , j+1, i  ), y(k  , j+1, i  ), z(k  , j+1, i  ) };
    cornerCoords[3] = { x(k  , j  , i  ), y(k  , j  , i  ), z(k  , j  , i  ) };
    cornerCoords[4] = { x(k+1, j  , i+1), y(k+1, j  , i+1), z(k+1, j  , i+1) };
    cornerCoords[5] = { x(k+1, j+1, i+1), y(k+1, j+1, i+1), z(k+1, j+1, i+1) };
    cornerCoords[6] = { x(k+1, j+1, i  ), y(k+1, j+1, i  ), z(k+1, j+1, i  ) };
    cornerCoords[7] = { x(k+1, j  , i  ), y(k+1, j  , i  ), z(k+1, j  , i  ) };

    cornerValues[0] = m_fcnView(k  , j  , i+1);
    cornerValues[1] = m_fcnView(k  , j+1, i+1);
    cornerValues[2] = m_fcnView(k  , j+1, i  );
    cornerValues[3] = m_fcnView(k  , j  , i  );
    cornerValues[4] = m_fcnView(k+1, j  , i+1);
    cornerValues[5] = m_fcnView(k+1, j+1, i+1);
    cornerValues[6] = m_fcnView(k+1, j+1, i  );
    cornerValues[7] = m_fcnView(k+1, j  , i  );
    // clang-format on
  }

  //!@brief Output contour mesh to a mint::UnstructuredMesh object.
  void populateContourMesh(
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

    const axom::IndexType addedCellCount = m_contourCellCorners.size();
    const axom::IndexType addedNodeCount = m_contourNodeCoords.size();
    if(addedCellCount != 0)
    {
      const axom::IndexType priorCellCount = mesh.getNumberOfCells();
      const axom::IndexType priorNodeCount = mesh.getNumberOfNodes();
      mesh.reserveNodes(priorNodeCount + addedNodeCount);
      mesh.reserveCells(priorCellCount + addedCellCount);

      mesh.appendNodes((double*)m_contourNodeCoords.data(),
                       m_contourNodeCoords.size());
      for(int n = 0; n < addedCellCount; ++n)
      {
        MIdx cornerIds = m_contourCellCorners[n];
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
          multidim_cell_index(m_contourCellParents[i]);
      }
    }
  }

  //!@brief Interpolate for the contour location crossing a parent edge.
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
    constexpr double ptiny = axom::primal::PRIMAL_TINY;
    const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
    const double w = (m_contourVal - f1) / df;
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
    constexpr double ptiny = axom::primal::PRIMAL_TINY;
    const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
    const double w = (m_contourVal - f1) / df;
    for(int d = 0; d < DIM; ++d)
    {
      crossingPt[d] = p1[d] + w * (p2[d] - p1[d]);
    }
  }

  //!@brief Set a value to find the contour for.
  void setContourValue(double contourVal) override
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
    : m_crossings(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_contourNodeCoords(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_contourCellCorners(0, 0, execution_space<ExecSpace>::allocatorID())
    , m_contourCellParents(0, 0, execution_space<ExecSpace>::allocatorID())
  { }

  /*!
    @brief Info for a parent cell intersecting the contour surface.
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

private:
  MIdx m_bShape;  //!< @brief Blueprint cell data shape.
  MIdx m_cShape;  //!< @brief Cell-centered array shape for ArrayViews.
  MIdx m_pShape;  //!< @brief Node-centered array shape for ArrayViews.
  axom::IndexType m_bStrides[DIM];  //!< @brief Strides for m_bShape arrays.

  // Views of parent domain data.
  axom::ArrayView<const double, DIM> m_coordsViews[DIM];
  axom::ArrayView<const double, DIM> m_fcnView;
  axom::ArrayView<const int, DIM> m_maskView;

  //!@brief Crossing case for each computational mesh cell.
  axom::Array<std::uint16_t, DIM, MemorySpace> m_caseIds;

  //!@brief Info on every parent cell that crosses the contour surface.
  axom::Array<CrossingInfo> m_crossings;

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
  axom::Array<Point> m_contourNodeCoords;

  //!@brief Corners (index into m_contourNodeCoords) of generated contour cells.
  axom::Array<MIdx> m_contourCellCorners;

  //!@brief Flat index of computational cell crossing the contour cell.
  axom::Array<IndexType> m_contourCellParents;
  //@}

  double m_contourVal = 0.0;
};

}  // end namespace marching_cubes
}  // end namespace detail
}  // end namespace quest
}  // end namespace axom
