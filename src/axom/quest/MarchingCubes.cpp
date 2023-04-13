// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/marching_cubes_lookup.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
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

//!@brief Add scalar value to every component in StackArray.
template<typename T, int DIM>
static axom::StackArray<T, DIM> operator+(const axom::StackArray<T, DIM>& left, T right)
{
  axom::StackArray<T, DIM> rval = left;
  for(int d=0; d<DIM; ++d) rval[d] += right;
  return rval;
}

//!@brief Reverse the order of a StackArray.
template<typename T, int DIM>
void reverse(axom::StackArray<T, DIM>& a)
{
  for(int d=0; d<DIM/2; ++d)
  {
    std::swap(a[d], a[DIM-1-d]);
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
    impl.generate_surface();
    impl.populate_surface_mesh(mesh, m_cellIdField);
  @endverbatim
*/
template<int DIM, typename ExecSpace>
struct MarchingCubesImpl {
  using Point = axom::primal::Point<double, DIM>;
  using MdimIdx = axom::StackArray<axom::IndexType, DIM>;
  /*!
    @brief Initialize data to a blueprint domain.
    @param dom Blueprint structured domain
    @param coordsetPath Where coordinates are in dom
    @param fcnPath Where nodal function is in dom
    @param maskPath Where cell mask function is in dom

    Set up views to external data and allocate internal data to work
    on the given domain.
  */
  AXOM_HOST void initialize(const conduit::Node& dom,
                            const std::string& coordsetPath,
                            const std::string& fcnPath,
                            const std::string& maskPath)
    {
      clear();
      m_dom = &dom;

      // Data sizes
      const conduit::Node& dimsNode =
        m_dom->fetch_existing("topologies/mesh/elements/dims");
      for(int d = 0; d < DIM; ++d)
      {
        m_bShape[d] = dimsNode[d].as_int();
      }
      m_cShape = m_bShape;
      reverse(m_cShape);
      m_pShape = m_cShape + 1;

      // Domain's node coordinates
      {
        const conduit::Node& coordValues =
          m_dom->fetch_existing(coordsetPath + "/values");
        bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
        const double* coordsPtrs[DIM];
        for(int d=0; d<DIM; ++d) coordsPtrs[d] = coordValues[d].as_double_ptr();
        const int coordSp = isInterleaved ? DIM : 1;
        for( int d=0; d<DIM; ++d )
        {
          m_coordsViews[d] = axom::ArrayView<const double, DIM>(coordsPtrs[d], m_pShape, coordSp);
        }
      }

      // Nodal function
      {
        auto& fcnValues = m_dom->fetch_existing(fcnPath + "/values");
        const double* fcnPtr = fcnValues.as_double_ptr();
        m_fcnView = axom::ArrayView<const double, DIM>(fcnPtr, m_pShape);
      }

      // Mask
      {
        const int* maskPtr = nullptr;
        if(!maskPath.empty())
        {
          auto& maskValues = m_dom->fetch_existing(maskPath + "/values");
          maskPtr = maskValues.as_int_ptr();
        }
        if (maskPtr)
        {
          m_maskView = axom::ArrayView<const int, DIM>(maskPtr, m_cShape);
        }
      }

      m_caseIds = axom::Array<std::uint16_t, DIM, axom::MemorySpace::Dynamic>(m_cShape);
    }

  //!@brief Populate m_caseIds with crossing indices.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type mark_crossings()
    {
      for(int j = 0; j < m_cShape[0]; ++j)
      {
        for(int i = 0; i < m_cShape[1]; ++i)
        {
          const bool skipZone = !m_maskView.empty() && bool(m_maskView(i, j));
          if(!skipZone)
          {
            // clang-format on
            double nodalValues[CELL_CORNER_COUNT] =
              { nodalValues[0] = m_fcnView(j  , i  )
              , nodalValues[1] = m_fcnView(j  , i+1)
              , nodalValues[2] = m_fcnView(j+1, i+1)
              , nodalValues[3] = m_fcnView(j+1, i  ) };
            // clang-format off

            auto crossingCase = compute_crossing_case(nodalValues);
            m_caseIds(j,i) = crossingCase;
          }
        }
      }
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type mark_crossings()
    {
      for(int k = 0; k < m_cShape[0]; ++k)
      {
        for(int j = 0; j < m_cShape[1]; ++j)
        {
          for(int i = 0; i < m_cShape[2]; ++i)
          {
            const bool skipZone = !m_maskView.empty() && bool(m_maskView(k, j, i));
            if(!skipZone)
            {
              double nodalValues[CELL_CORNER_COUNT];

              // clang-format off
              nodalValues[0] = m_fcnView(k    , j    , i + 1);
              nodalValues[1] = m_fcnView(k    , j + 1, i + 1);
              nodalValues[2] = m_fcnView(k    , j + 1, i    );
              nodalValues[3] = m_fcnView(k    , j    , i    );
              nodalValues[4] = m_fcnView(k + 1, j    , i + 1);
              nodalValues[5] = m_fcnView(k + 1, j + 1, i + 1);
              nodalValues[6] = m_fcnView(k + 1, j + 1, i    );
              nodalValues[7] = m_fcnView(k + 1, j    , i    );
              // clang-format on

              auto crossingCase = compute_crossing_case(nodalValues);
              m_caseIds(k,j,i) = crossingCase;
            }
          }
        }
      }
    }

  /*!
    @brief Scan m_caseIds to count number of crossing cells and determine
    some offset for numbering the generated nodes.
  */
  void scan_crossings()
    {
      // Number of cells a crossing can generate:
      const int* crossingCellCounts = DIM == 2 ?
        detail::num_segments : detail::num_triangles;

      // Reserve memory for crossing info.
      axom::IndexType crossingCount = 0;
      const auto cellCount = m_caseIds.size();
      for(int n=0; n<cellCount; ++n)
      {
        auto caseId = m_caseIds.flatIndex(n);
        auto surfaceCellCount = crossingCellCounts[caseId];
        crossingCount += bool(surfaceCellCount);
      }
      m_crossings.reserve(crossingCount);

      // Populate crossing info.
      axom::IndexType firstSurfaceCellIdx = 0;
      axom::IndexType parentCellCount = m_caseIds.size();
      for(axom::IndexType parentCellNum=0; parentCellNum<parentCellCount; ++parentCellNum)
      {
        auto caseId = m_caseIds.flatIndex(parentCellNum);
        auto surfaceCellCount = crossingCellCounts[caseId];
        if(surfaceCellCount != 0)
        {
          m_crossings.push_back({parentCellNum, caseId, firstSurfaceCellIdx});
          firstSurfaceCellIdx += surfaceCellCount;
        }
      }
      SLIC_ASSERT(m_crossings.size() == crossingCount);
      std::cout << __WHERE << m_crossings.size() << " crossings found." << std::endl;
    }

  //!@brief Compute multidimensional index from flat cell index.
  axom::StackArray<axom::IndexType, DIM> multidim_cell_index(
    axom::IndexType flatId)
    {
      axom::IndexType strides[DIM] = {1};
      for( int d=1; d<DIM; ++d) strides[d] = strides[d-1]*m_bShape[d-1];
      if(flatId >= strides[DIM-1]*m_bShape[DIM-1])
      {
        SLIC_ERROR("flatId is too big.");
      }

      axom::StackArray<axom::IndexType, DIM> rval;
      for(int d=DIM-1; d>=0; --d)
      {
        rval[d] = flatId/strides[d];
        flatId -= rval[d]*strides[d];
      }
      return rval;
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type generate_surface()
    {
      if(m_crossings.empty())
      {
        return;
      }
      // Determine needed capacity in output mesh and reserve.
      // TODO: For multidomain mesh, capacity should be summed over domains
      // and memory pre-allocated outside this class.
      axom::IndexType surfaceCellCount = m_crossings.back().firstSurfaceCellIdx
        + detail::num_segments[m_crossings.back().caseNum];
      axom::IndexType surfaceNodeCount = DIM*surfaceCellCount;

      m_surfaceNodeCoords.clear();
      m_surfaceCellCorners.clear();
      m_surfaceCellParents.clear();
      m_surfaceNodeCoords.resize(surfaceNodeCount);
      m_surfaceCellCorners.resize(surfaceCellCount);
      m_surfaceCellParents.resize(surfaceCellCount);

      IndexType surfaceCellId = 0;
      IndexType surfaceNodeId = 0;
      for(axom::IndexType iCrossing=0; iCrossing<m_crossings.size(); ++iCrossing)
      {
        const auto& crossingInfo = m_crossings[iCrossing];
        IndexType surfaceCellCount = detail::num_segments[crossingInfo.caseNum];
        surfaceCellId = crossingInfo.firstSurfaceCellIdx;
        surfaceNodeId = crossingInfo.firstSurfaceCellIdx*DIM;
        SLIC_ASSERT(surfaceCellCount > 0);

        const auto c = multidim_cell_index(crossingInfo.parentCellNum);
        const auto& i = c[0];
        const auto& j = c[1];

        // clang-format off
        Point cornerCoords[CELL_CORNER_COUNT] =
          { {m_coordsViews[0](j  ,i  ), m_coordsViews[1](j  ,i  )}
          , {m_coordsViews[0](j  ,i+1), m_coordsViews[1](j  ,i+1)}
          , {m_coordsViews[0](j+1,i+1), m_coordsViews[1](j+1,i+1)}
          , {m_coordsViews[0](j+1,i  ), m_coordsViews[1](j+1,i  )} };

        double cornerValues[CELL_CORNER_COUNT] =
          { m_fcnView(j  ,i  )
          , m_fcnView(j  ,i+1)
          , m_fcnView(j+1,i+1)
          , m_fcnView(j+1,i  ) };
        // clang-format on

        for(int iCell = 0; iCell < surfaceCellCount; ++iCell)
        {
          const int e1 = detail::cases2D[crossingInfo.caseNum][iCell * 2];
          const int e2 = detail::cases2D[crossingInfo.caseNum][iCell * 2 + 1];
          Point p; // double p[2];

          linear_interp(e1, cornerCoords, cornerValues, p);
          m_surfaceNodeCoords[surfaceNodeId][0] = p[0];
          m_surfaceNodeCoords[surfaceNodeId][1] = p[1];
          m_surfaceCellCorners[surfaceCellId][0] = surfaceNodeId;
          ++surfaceNodeId;

          linear_interp(e2, cornerCoords, cornerValues, p);
          m_surfaceNodeCoords[surfaceNodeId][0] = p[0];
          m_surfaceNodeCoords[surfaceNodeId][1] = p[1];
          m_surfaceCellCorners[surfaceCellId][1] = surfaceNodeId;
          ++surfaceNodeId;

          m_surfaceCellParents[surfaceCellId] = crossingInfo.parentCellNum;

          ++surfaceCellId;
        }  // END for all segments
      }
      assert(surfaceCellId == surfaceCellCount);
      assert(surfaceNodeId == surfaceNodeCount);
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type generate_surface()
    {
      if(m_crossings.empty())
      {
        return;
      }
      // Determine needed capacity in output mesh and reserve.
      // TODO: For multidomain mesh, capacity should be summed over domains
      // and memory pre-allocated outside this class.
      axom::IndexType surfaceCellCount = m_crossings.back().firstSurfaceCellIdx
        + detail::num_triangles[m_crossings.back().caseNum];
      axom::IndexType surfaceNodeCount = DIM*surfaceCellCount;

      m_surfaceNodeCoords.clear();
      m_surfaceCellCorners.clear();
      m_surfaceCellParents.clear();
      m_surfaceNodeCoords.resize(surfaceNodeCount);
      m_surfaceCellCorners.resize(surfaceCellCount);
      m_surfaceCellParents.resize(surfaceCellCount);

      IndexType surfaceCellId = 0;
      IndexType surfaceNodeId = 0;
      for(axom::IndexType iCrossing=0; iCrossing<m_crossings.size(); ++iCrossing)
      {
        const auto& crossingInfo = m_crossings[iCrossing];
        IndexType surfaceCellCount = detail::num_triangles[crossingInfo.caseNum];
        surfaceCellId = crossingInfo.firstSurfaceCellIdx;
        surfaceNodeId = crossingInfo.firstSurfaceCellIdx*DIM;
        SLIC_ASSERT(surfaceCellCount > 0);

        const auto c = multidim_cell_index(crossingInfo.parentCellNum);
        const auto& i = c[0];
        const auto& j = c[1];
        const auto& k = c[2];

        // clang-format off
        Point cornerCoords[CELL_CORNER_COUNT] =
          { {m_coordsViews[0](k  ,j  ,i+1), m_coordsViews[1](k  ,j  ,i+1), m_coordsViews[2](k  ,j  ,i+1)}
          , {m_coordsViews[0](k  ,j+1,i+1), m_coordsViews[1](k  ,j+1,i+1), m_coordsViews[2](k  ,j+1,i+1)}
          , {m_coordsViews[0](k  ,j+1,i  ), m_coordsViews[1](k  ,j+1,i  ), m_coordsViews[2](k  ,j+1,i  )}
          , {m_coordsViews[0](k  ,j  ,i  ), m_coordsViews[1](k  ,j  ,i  ), m_coordsViews[2](k  ,j  ,i  )}
          , {m_coordsViews[0](k+1,j  ,i+1), m_coordsViews[1](k+1,j  ,i+1), m_coordsViews[2](k+1,j  ,i+1)}
          , {m_coordsViews[0](k+1,j+1,i+1), m_coordsViews[1](k+1,j+1,i+1), m_coordsViews[2](k+1,j+1,i+1)}
          , {m_coordsViews[0](k+1,j+1,i  ), m_coordsViews[1](k+1,j+1,i  ), m_coordsViews[2](k+1,j+1,i  )}
          , {m_coordsViews[0](k+1,j  ,i  ), m_coordsViews[1](k+1,j  ,i  ), m_coordsViews[2](k+1,j  ,i  )} };

        double cornerValues[CELL_CORNER_COUNT] =
          { m_fcnView(k  ,j  ,i+1)
          , m_fcnView(k  ,j+1,i+1)
          , m_fcnView(k  ,j+1,i  )
          , m_fcnView(k  ,j  ,i  )
          , m_fcnView(k+1,j  ,i+1)
          , m_fcnView(k+1,j+1,i+1)
          , m_fcnView(k+1,j+1,i  )
          , m_fcnView(k+1,j  ,i  ) };
        // clang-format on

        for(int iCell = 0; iCell < surfaceCellCount; ++iCell)
        {
          const int e1 = detail::cases3D[crossingInfo.caseNum][iCell * 3];
          const int e2 = detail::cases3D[crossingInfo.caseNum][iCell * 3 + 1];
          const int e3 = detail::cases3D[crossingInfo.caseNum][iCell * 3 + 2];
          Point p; // double p[3];

          linear_interp(e1, cornerCoords, cornerValues, p);
          m_surfaceNodeCoords[surfaceNodeId] = p;
          m_surfaceCellCorners[surfaceCellId][0] = surfaceNodeId;
          ++surfaceNodeId;

          linear_interp(e2, cornerCoords, cornerValues, p);
          m_surfaceNodeCoords[surfaceNodeId] = p;
          m_surfaceCellCorners[surfaceCellId][1] = surfaceNodeId;
          ++surfaceNodeId;

          linear_interp(e3, cornerCoords, cornerValues, p);
          m_surfaceNodeCoords[surfaceNodeId] = p;
          m_surfaceCellCorners[surfaceCellId][2] = surfaceNodeId;
          ++surfaceNodeId;

          m_surfaceCellParents[surfaceCellId] = crossingInfo.parentCellNum;

          ++surfaceCellId;
        }  // END for all segments
      }
      assert(surfaceCellId == surfaceCellCount);
      assert(surfaceNodeId == surfaceNodeCount);
    }

  void populate_surface_mesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> *mesh,
    const std::string& cellIdField)
  {
    const axom::IndexType addedCellCount = m_surfaceCellCorners.size();
    const axom::IndexType addedNodeCount = m_surfaceNodeCoords.size();
    if(addedCellCount != 0) {
      const axom::IndexType priorCellCount = mesh->getNumberOfCells();
      const axom::IndexType priorNodeCount = mesh->getNumberOfNodes();
      mesh->reserveNodes(priorNodeCount + addedNodeCount);
      mesh->reserveCells(priorCellCount + addedCellCount);

      mesh->appendNodes((double*)m_surfaceNodeCoords.data(),
                        m_surfaceNodeCoords.size());
      for( int n=0; n<addedCellCount; ++n )
      {
        const MdimIdx
          cornerIds = m_surfaceCellCorners[n] + priorNodeCount;
        mesh->appendCell(cornerIds);
      }
      // TODO: Replace with device-capable copy.
      axom::IndexType* dst =
        mesh->getFieldPtr<axom::IndexType>(cellIdField, axom::mint::CELL_CENTERED);
      memcpy(dst + priorCellCount,
             m_surfaceCellParents.data(),
             sizeof(axom::IndexType)*addedCellCount);
    }
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type linear_interp(
    int edgeIdx,
    const Point cornerCoords[4],
    const double nodeValues[4],
    Point& xyz)
    {
      // STEP 0: get the edge node indices
      // 2 nodes define the edge.  n1 and n2 are the indices of
      // the nodes w.r.t. the square or cubic zone.  There is a
      // agreed-on ordering of these indices in the arrays xx, yy,
      // zz, nodeValues, xyz.
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
        xyz = p1; // memcpy(xyz, p1, DIM * sizeof(double));
        return;
      }

      if(axom::utilities::isNearlyEqual(m_contourVal, f2))
      {
        xyz = p2; // memcpy(xyz, p2, DIM * sizeof(double));
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      constexpr double ptiny = 1.0e-80;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (m_contourVal - f1) / df;
      for(int d = 0; d < DIM; ++d)
      {
        xyz[d] = p1[d] + w * (p2[d] - p1[d]);
      }
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type linear_interp(
    int edgeIdx,
    const Point cornerCoords[8],
    const double nodeValues[8],
    Point& xyz)
    {
      // STEP 0: get the edge node indices
      // 2 nodes define the edge.  n1 and n2 are the indices of
      // the nodes w.r.t. the square or cubic zone.  There is a
      // agreed-on ordering of these indices in the arrays xx, yy,
      // zz, nodeValues, xyz.
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
        xyz = p1; // memcpy(&xyz, p1, DIM * sizeof(double));
        return;
      }

      if(axom::utilities::isNearlyEqual(m_contourVal, f2))
      {
        xyz = p2; // memcpy(&xyz, p2, DIM * sizeof(double));
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      constexpr double ptiny = 1.0e-80;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (m_contourVal - f1) / df;
      for(int d = 0; d < DIM; ++d)
      {
        xyz[d] = p1[d] + w * (p2[d] - p1[d]);
      }
// std::cout<<__WHERE<< xyz << std::endl;
    }

  void set_contour_value(double contourVal)
    {
      m_contourVal = contourVal;
    }

  //!@brief Compute the case index into case2D or case3D.
  int compute_crossing_case(const double* f)
    {
      int index = 0;
      for(int n = 0; n < CELL_CORNER_COUNT; ++n)
      {
        if(f[n] >= m_contourVal)
        {
          const int mask = (1 << n);
          index |= mask;
        }
      }
      return (index);
    }

  //!@brief Clear data so you can rerun with a different contour value.
  void clear()
  {
    m_surfaceNodeCoords.clear();
    m_surfaceCellCorners.clear();
    m_surfaceCellParents.clear();
  }

  /*!
    @brief Info for a parent cell intersecting the surface.
  */
  struct CrossingInfo {
    CrossingInfo(axom::IndexType parentCellNum_,
                 std::uint16_t caseNum_,
                 axom::IndexType firstSurfaceCellIdx_)
      : parentCellNum(parentCellNum_)
      , caseNum(caseNum_)
      , firstSurfaceCellIdx(firstSurfaceCellIdx_) {}
    axom::IndexType parentCellNum; //!< @brief Flat index of parent cell.
    std::uint16_t caseNum;   //!< @brief Index in cases2D or cases3D
    axom::IndexType firstSurfaceCellIdx; //!< @brief First index for generated cells.
  };
  axom::Array<CrossingInfo> m_crossings;

  const conduit::Node *m_dom;

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

  //!@brief Number of corners (nodes) on each cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;

  //!@name Output surface mesh.
  //@{
  //!@brief Coordinates of generated surface nodes.
  axom::Array<Point> m_surfaceNodeCoords;

  //!@brief Corners (index into m_surfaceNodeCoords) of generated surface cells.
  axom::Array<MdimIdx> m_surfaceCellCorners;

  //!@brief Computational cell (flat index) crossing the surface cell.
  axom::Array<IndexType> m_surfaceCellParents;
  //@}

  double m_contourVal = 0.0;
};


MarchingCubes::MarchingCubes(const conduit::Node& bpMesh,
                             const std::string& coordsetName,
                             const std::string& maskField)
  : m_sd()
  , m_ndim(0)
  , m_coordsetPath("coordsets/" + coordsetName)
  , m_fcnPath()
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
  , m_surfaceMesh(nullptr)
  , m_cellIdField()
  , m_domainIdField()
{
  m_sd.reserve(conduit::blueprint::mesh::number_of_domains(bpMesh));
  for(auto& dom : bpMesh.children())
  {
    m_sd.emplace_back(new MarchingCubesSingleDomain(dom, coordsetName, maskField));
    if(m_ndim == 0)
    {
      m_ndim = m_sd.back()->dimension();
    }
    else
    {
      SLIC_ASSERT(m_ndim == m_sd.back()->dimension());
    }
  }
}

void MarchingCubes::set_function_field(const std::string& fcnField)
{
  m_fcnPath = "fields/" + fcnField;
  for(auto& s : m_sd)
  {
    s->set_function_field(fcnField);
  }
}

/*!
  @brief Set the output surface mesh object.
*/
void MarchingCubes::set_output_mesh(axom::mint::Mesh* surfaceMesh)
{
  m_surfaceMesh = surfaceMesh;

  for(auto& s : m_sd)
  {
    s->set_output_mesh(m_surfaceMesh);
  }
}

void MarchingCubes::compute_iso_surface(double contourVal)
{
  SLIC_ASSERT_MSG(m_surfaceMesh,
                  "You must call set_output_mesh before compute_iso_surface.");
  SLIC_ASSERT_MSG(
    !m_fcnPath.empty(),
    "You must call set_function_field before compute_iso_surface.");

  if(!m_domainIdField.empty() &&
     !m_surfaceMesh->hasField(m_domainIdField, axom::mint::CELL_CENTERED))
  {
    m_surfaceMesh->createField<axom::IndexType>(m_domainIdField, axom::mint::CELL_CENTERED);
  }

  for(int dId = 0; dId < m_sd.size(); ++dId)
  {
    std::shared_ptr<MarchingCubesSingleDomain>& single = m_sd[dId];

    auto nPrev = m_surfaceMesh->getNumberOfCells();
    single->compute_iso_surface(contourVal);
    auto nNew = m_surfaceMesh->getNumberOfCells();

    if(nNew > nPrev && !m_domainIdField.empty())
    {
      auto* domainIdPtr =
        m_surfaceMesh->getFieldPtr<axom::IndexType>(m_domainIdField,
                                        axom::mint::CELL_CENTERED);
      for(int n = nPrev; n < nNew; ++n)
      {
        domainIdPtr[n] = dId;
      }
    }
  }
}

MarchingCubesSingleDomain::MarchingCubesSingleDomain(const conduit::Node& dom,
                                                     const std::string& coordsetName,
                                                     const std::string& maskField)
  : m_dom(nullptr)
  , m_ndim(0)
  , m_cShape()
  , m_logicalOrigin()
  , m_coordsetPath("coordsets/" + coordsetName)
  , m_fcnPath()
  , m_maskPath(maskField.empty() ? std::string() : "fields/" + maskField)
  , m_surfaceMesh(nullptr)
  , m_cellIdField()
  , _contourVal(0.0)
{
  set_domain(dom);
  return;
}

void MarchingCubesSingleDomain::set_domain(const conduit::Node& dom)
{
  SLIC_ASSERT_MSG(
    !conduit::blueprint::mesh::is_multi_domain(dom),
    "MarchingCubesSingleDomain is single-domain only.  Try MarchingCubes.");

  SLIC_ASSERT(dom.has_path(m_coordsetPath));
  SLIC_ASSERT(dom["topologies/mesh/type"].as_string() == "structured");

  if(!m_maskPath.empty())
  {
    SLIC_ASSERT(dom.has_path(m_maskPath + "/values"));
  }

  m_dom = &dom;

  const conduit::Node& dimsNode =
    m_dom->fetch_existing("topologies/mesh/elements/dims");

  m_ndim = dimsNode.number_of_children();

  m_cShape.resize(m_ndim);
  for(int d = 0; d < m_ndim; ++d)
  {
    m_cShape[d] = dimsNode[m_ndim - 1 - d].as_int();
  }

  m_logicalOrigin.resize(m_ndim, 0);
  if(m_dom->has_path("topologies/mesh/elements/origin"))
  {
    const conduit::Node& origins =
      m_dom->fetch_existing("topologies/mesh/elements/origin");
    for(int d = 0; d < m_ndim; ++d)
    {
      m_logicalOrigin[d] = origins[d].as_int();
    }
  }

  SLIC_ASSERT(m_ndim >= 1 && m_ndim <= 3);

  const conduit::Node& coordsValues = dom[m_coordsetPath + "/values"];
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);
  SLIC_ASSERT_MSG(
    !isInterleaved,
    "MarchingCubes currently requires contiguous coordinates layout.");
}

void MarchingCubesSingleDomain::set_function_field(const std::string& fcnField)
{
  m_fcnPath = "fields/" + fcnField;
  SLIC_ASSERT(m_dom->has_path(m_fcnPath));
  SLIC_ASSERT(m_dom->fetch_existing(m_fcnPath + "/association").as_string() ==
              "vertex");
  SLIC_ASSERT(m_dom->has_path(m_fcnPath + "/values"));
}

/*!
  @brief Set the output surface mesh object.
*/
void MarchingCubesSingleDomain::set_output_mesh(axom::mint::Mesh* surfaceMesh)
{
  m_surfaceMesh = surfaceMesh;

  if(m_ndim == 2)
  {
    using SegmentMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
    SegmentMesh* mesh = dynamic_cast<SegmentMesh*>(m_surfaceMesh);
    SLIC_ASSERT_MSG(mesh, "Surface mesh for 2D problem must be a SegmentMesh");
  }
  else
  {
    using TriangleMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
    TriangleMesh* mesh = dynamic_cast<TriangleMesh*>(m_surfaceMesh);
    SLIC_ASSERT_MSG(mesh, "Surface mesh for 3D problem must be a TriangleMesh");
  }
}

void MarchingCubesSingleDomain::compute_iso_surface(double contourVal)
{
  using OutputMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
  OutputMesh* mesh = static_cast<OutputMesh*>(m_surfaceMesh);

  SLIC_ASSERT_MSG(m_surfaceMesh,
                  "You must call set_output_mesh before compute_iso_surface.");
  SLIC_ASSERT_MSG(
    !m_fcnPath.empty(),
    "You must call set_function_field before compute_iso_surface.");

  const int* maskPtr = nullptr;
  if(!m_maskPath.empty())
  {
    auto& maskValues = m_dom->fetch_existing(m_maskPath + "/values");
    maskPtr = maskValues.as_int_ptr();
  }

  if(!m_cellIdField.empty() &&
     !m_surfaceMesh->hasField(m_cellIdField, axom::mint::CELL_CENTERED))
  {
    m_surfaceMesh->createField<axom::IndexType>(m_cellIdField, axom::mint::CELL_CENTERED);
  }

  _contourVal = contourVal;

  if(m_ndim == 2)
  {
    MarchingCubesImpl<2, axom::execution_space<axom::SEQ_EXEC>> impl2d;
    impl2d.initialize(*m_dom, m_coordsetPath, m_fcnPath, m_maskPath);
    impl2d.set_contour_value(contourVal);
    impl2d.mark_crossings();
    impl2d.scan_crossings();
    impl2d.generate_surface();
    impl2d.populate_surface_mesh(mesh, m_cellIdField);
  }
  else
  {
    MarchingCubesImpl<3, axom::execution_space<axom::SEQ_EXEC>> impl3d;
    impl3d.initialize(*m_dom, m_coordsetPath, m_fcnPath, m_maskPath);
    impl3d.set_contour_value(contourVal);
    impl3d.mark_crossings();
    impl3d.scan_crossings();
    impl3d.generate_surface();
    impl3d.populate_surface_mesh(mesh, m_cellIdField);
  }  // END else
}

}  // end namespace quest
}  // end namespace axom
