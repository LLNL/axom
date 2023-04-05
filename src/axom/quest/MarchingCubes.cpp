// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/detail/marching_cubes_lookup.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
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

/*!
  @brief Computations for MarchingCubesSingleDomain, limited to what
  can be run on both hosts and devices.
*/
template<int DIM, typename ExecSpace>
struct MarchingCubesSingleDomainImpl {
  AXOM_HOST void set_domain(const conduit::Node& dom,
                            const std::string& coordsetPath,
                            const std::string& fcnPath,
                            const std::string& maskPath)
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
      m_dom = &dom;

      const conduit::Node& coordValues =
        m_dom->fetch_existing(coordsetPath + "/values");
      bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);

      auto& fcnValues = m_dom->fetch_existing(fcnPath + "/values");
      const double* fcnPtr = fcnValues.as_double_ptr();

      const int* maskPtr = nullptr;
      if(!maskPath.empty())
      {
        auto& maskValues = m_dom->fetch_existing(maskPath + "/values");
        maskPtr = maskValues.as_int_ptr();
      }

      const conduit::Node& dimsNode =
        m_dom->fetch_existing("topologies/mesh/elements/dims");
      for(int d = 0; d < DIM; ++d)
      {
        m_cShape[d] = dimsNode[DIM - 1 - d].as_int();
      }

      const double* coordsPtrs[DIM];
      for(int d=0; d<DIM; ++d) coordsPtrs[d] = coordValues[d].as_double_ptr();

      set_domain(m_cShape, coordsPtrs, fcnPtr, maskPtr, isInterleaved);
    }
  AXOM_HOST_DEVICE void set_domain(const axom::StackArray<axom::IndexType, DIM>& cShape,
                                   const double* coordsPtrs[DIM],
                                   const double* fcnPtr,
                                   const int* maskPtr,
                                   bool isInterleaved)
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
      // TODO: check that memory is compatible with ExecSpace
      const int coordSp = isInterleaved ? DIM : 1;
      m_cShape = cShape;
      m_fcnPtr = fcnPtr;
      m_maskPtr = maskPtr;
      m_checkMask = bool(maskPtr);
      m_computationalCellCount = 1;
      for( int i=0; i<DIM; ++i )
      {
        SLIC_ASSERT(m_cShape[i] > 0);
        m_computationalCellCount *= m_cShape[i];
        m_pShape[i] = 1 + m_cShape[i];
        m_coordsPtrs[i] = coordsPtrs[i];
      }
      for( int i=0; i<DIM; ++i )
      {
        m_coordsViews[i] = axom::ArrayView<const double, DIM>(coordsPtrs[i], m_pShape, coordSp);
      }
      m_fcnView = axom::ArrayView<const double, DIM>(fcnPtr, m_pShape);
      if (m_maskPtr)
      {
        m_maskView = axom::ArrayView<const int, DIM>(maskPtr, m_cShape);
      }
      m_caseIds = axom::Array<std::uint16_t, DIM, axom::MemorySpace::Dynamic>(m_cShape);
      set_offset_pointers();
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type set_offset_pointers()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
      for(int d=0; d<DIM; ++d)
      {
        m_cornerCoordsPtrs[d][0] = &m_coordsViews[d](0 + 1, 0);
        m_cornerCoordsPtrs[d][1] = &m_coordsViews[d](0 + 1, 0 + 1);
        m_cornerCoordsPtrs[d][2] = &m_coordsViews[d](0, 0 + 1);
        m_cornerCoordsPtrs[d][3] = &m_coordsViews[d](0, 0);
      }
      m_cornerFcnPtrs[0] = &m_fcnView(0 + 1, 0);
      m_cornerFcnPtrs[1] = &m_fcnView(0 + 1, 0 + 1);
      m_cornerFcnPtrs[2] = &m_fcnView(0, 0 + 1);
      m_cornerFcnPtrs[3] = &m_fcnView(0, 0);
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type mark_crossings()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
      for(int i = 0; i < m_cShape[0]; ++i)
      {
        for(int j = 0; j < m_cShape[1]; ++j)
        {
          const bool skipZone = m_checkMask && bool(m_maskView(i, j));
          if(!skipZone)
          {
            double nodalValues[4];

            nodalValues[0] = m_fcnView(i + 1, j);
            nodalValues[1] = m_fcnView(i + 1, j + 1);
            nodalValues[2] = m_fcnView(i, j + 1);
            nodalValues[3] = m_fcnView(i, j);

            auto crossingCase = computeIndex(nodalValues);
            m_caseIds(i,j) = crossingCase;

          }  // END if
        }
      }
    }
  //! @brief Scan m_caseIds to count number of cell intersected and dependent data.
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type scan_crossings()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;

      // Reserve memory for crossing info.
      auto* caseFlatIds = m_caseIds.data();
      const auto cellCount = m_caseIds.size();
      axom::IndexType crossingCount = 0;
      for(int i=0; i<cellCount; ++i)
      {
        auto caseId = caseFlatIds[i];
        auto surfaceCellCount = detail::num_segments[caseId];
        crossingCount += bool(surfaceCellCount);
      }
      m_crossings.reserve(crossingCount);

      // Populate crossing info.
      axom::IndexType firstSurfaceCellIdx = 0;
      for(int i=0; i<cellCount; ++i)
      {
        auto caseId = caseFlatIds[i];
        auto surfaceCellCount = detail::num_segments[caseId];
        if(surfaceCellCount != 0)
        {
          m_crossings.push_back({i, caseId, firstSurfaceCellIdx});
          firstSurfaceCellIdx += surfaceCellCount;
        }
      }
      SLIC_ASSERT(m_crossings.size() == crossingCount);
      std::cout << __WHERE << m_crossings.size() << " crossings found." << std::endl;
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type generate_surface()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
      if(m_crossings.empty())
      {
        return;
      }

      // Generate line segments
      // using SegmentMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
      // SegmentMesh* mesh = static_cast<SegmentMesh*>(m_surfaceMesh);
      // SLIC_ASSERT(mesh != NULL);
      // SLIC_ASSERT(mesh->getCellType() == axom::mint::SEGMENT);
      // SLIC_ASSERT(mesh->getDimension() == 2);

      // Determine needed capacity in output mesh and reserve.
      // TODO: For multidomain mesh, capacity should be summed over domains
      // and memory pre-allocated outside this class.
      axom::IndexType surfaceCellCount = m_crossings.back().firstSurfaceCellIdx
        + detail::num_segments[m_crossings.back().caseIdx];
      axom::IndexType surfaceNodeCount = 2*surfaceCellCount;

      m_surfaceNodeCoords.clear();
      m_surfaceCellCorners.clear();
      m_surfaceCellParents.clear();
      m_surfaceNodeCoords.resize(surfaceNodeCount, DIM);
      m_surfaceCellCorners.resize(surfaceCellCount, DIM);
      m_surfaceCellParents.resize(surfaceCellCount);

      for(axom::IndexType iCase=0; iCase<m_crossings.size(); ++iCase)
      {
        const auto& caseInfo = m_crossings[iCase];
        IndexType compCellFlatIdx = caseInfo.compCellFlatIdx;
        IndexType nsegs = detail::num_segments[caseInfo.caseIdx];
        IndexType surfaceCellIdx = caseInfo.firstSurfaceCellIdx;
        IndexType surfaceNodeIdx = caseInfo.firstSurfaceCellIdx*DIM;
        SLIC_ASSERT(nsegs > 0);
        double nodalValues[CELL_CORNER_COUNT];
        double xx[CELL_CORNER_COUNT];
        double yy[CELL_CORNER_COUNT];
        for(int iCorner=0; iCorner<CELL_CORNER_COUNT; ++iCorner)
        {
          xx[iCorner] = m_cornerCoordsPtrs[0][iCorner][compCellFlatIdx];
          yy[iCorner] = m_cornerCoordsPtrs[1][iCorner][compCellFlatIdx];
          nodalValues[iCorner] = m_cornerFcnPtrs[iCorner][compCellFlatIdx];
        }
        for(int iSeg = 0; iSeg < nsegs; ++iSeg)
        {
          const int e1 = detail::cases2D[caseInfo.caseIdx][iSeg * 2];
          const int e2 = detail::cases2D[caseInfo.caseIdx][iSeg * 2 + 1];
          double p[2];

          linear_interp(e1, xx, yy, NULL, nodalValues, p);
          // mesh->appendNode(p[0], p[1]);
          // cell[0] = surfaceNodeCount;
          m_surfaceNodeCoords(surfaceNodeIdx, 0) = p[0];
          m_surfaceNodeCoords(surfaceNodeIdx, 1) = p[1];
          m_surfaceCellCorners(surfaceCellIdx, 0) = surfaceNodeCount;
          ++surfaceNodeIdx;

          linear_interp(e2, xx, yy, NULL, nodalValues, p);
          // mesh->appendNode(p[0], p[1]);
          // cell[1] = surfaceNodeCount;
          m_surfaceNodeCoords(surfaceNodeIdx, 0) = p[0];
          m_surfaceNodeCoords(surfaceNodeIdx, 1) = p[1];
          m_surfaceCellCorners(surfaceCellIdx, 1) = surfaceNodeCount;
          ++surfaceNodeIdx;

          // mesh->appendCell(cell);
          ++surfaceCellIdx;

        }  // END for all segments
      }
      std::cout << __WHERE << m_surfaceNodeCoords.shape() << ' ' << m_surfaceCellCorners.shape() << std::endl;
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type linear_interp(
    int edgeIdx,
    const double* xx,
    const double* yy,
    const double* zz,
    const double* nodeValues,
    double* xyz)
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

      double p1[3];
      p1[0] = xx[n1];
      p1[1] = yy[n1];
      p1[2] = 0.0;

      double p2[3];
      p2[0] = xx[n2];
      p2[1] = yy[n2];
      p2[2] = 0.0;

      if(DIM == 3)
      {
        // set the z--coordinate if in 3-D
        p1[2] = zz[n1];
        p2[2] = zz[n2];
      }

      // STEP 2: check whether the interpolated point is at one of the two corners.
      if(axom::utilities::isNearlyEqual(m_contourVal, f1) ||
         axom::utilities::isNearlyEqual(f1, f2))
      {
        memcpy(xyz, p1, DIM * sizeof(double));
        return;
      }

      if(axom::utilities::isNearlyEqual(m_contourVal, f2))
      {
        memcpy(xyz, p2, DIM * sizeof(double));
        return;
      }

      // STEP 3: point is in between the edge points, interpolate its position
      constexpr double ptiny = 1.0e-80;
      const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
      const double w = (m_contourVal - f1) / df;
      for(int i = 0; i < DIM; ++i)
      {
        xyz[i] = p1[i] + w * (p2[i] - p1[i]);
      }
    }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type set_offset_pointers()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
    }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type mark_crossings()
    {
      std::cout << __WHERE "dim is " << DIM << std::endl;
    }

  void set_contour_value(double contourVal)
    {
      m_contourVal = contourVal;
    }
  int computeIndex(const double* f)
    {
      int index = 0;
      for(int i = 0; i < CELL_CORNER_COUNT; ++i)
      {
        if(f[i] >= m_contourVal)
        {
          const int mask = (1 << i);
          index |= mask;
        }
      }
      return (index);
    }

  /*!
    @brief Info for a computational cell intersecting the surface.
  */
  struct CrossingInfo {
    CrossingInfo(axom::IndexType compCellFlatIdx,
                 std::uint16_t caseIdx,
                 axom::IndexType firstSurfaceCellIdx)
      : compCellFlatIdx (compCellFlatIdx)
      , caseIdx (caseIdx)
      , firstSurfaceCellIdx (firstSurfaceCellIdx) {}

    axom::IndexType compCellFlatIdx; //!< @brief Flat index of computational cell.
    std::uint16_t caseIdx;   //!< @brief Crossing topology index into cases2D or cases3D
    axom::IndexType firstSurfaceCellIdx; //!< @brief First index for generated cells.
  };
  axom::Array<CrossingInfo> m_crossings;

  const conduit::Node *m_dom;
  const std::string m_coordsetPath;
  std::string m_fcnPath;

  //!@brief Number of corners (nodes) on each cell.
  static constexpr std::uint8_t CELL_CORNER_COUNT = (DIM == 3) ? 8 : 4;
  axom::StackArray<axom::IndexType, DIM> m_cShape;  //!< @brief Shape of cell-centered data arrays in m_dom.
  axom::StackArray<axom::IndexType, DIM> m_pShape;  //!< @brief Shape of node-centered data arrays in m_dom.
  axom::IndexType m_computationalCellCount; //!< @brief Cell count in domain.

  const double* m_coordsPtrs[DIM];
  const double* m_fcnPtr;
  const int* m_maskPtr;
  axom::ArrayView<const double, DIM> m_coordsViews[DIM];
  axom::ArrayView<const double, DIM> m_fcnView;
  axom::ArrayView<const int, DIM> m_maskView;
  bool m_checkMask;

  //!@brief Crossing case for each computational mesh cell.
  // TODO: Put this in correct memory space.
  axom::Array<std::uint16_t, DIM, axom::MemorySpace::Dynamic> m_caseIds;

  //!@brief Coordinates of nodes on a cell.  Points inside m_coordsPtrs.
  const double* m_cornerCoordsPtrs[DIM][CELL_CORNER_COUNT];
  //!@brief Function value at nodes on a cell.  Points inside m_fcnPtr.
  const double* m_cornerFcnPtrs[CELL_CORNER_COUNT];

  // Output surface mesh.
  // Migrate away from mint mesh, which don't appear to be device-ported.
  //!@brief Coordinates of generated surface nodes.
  axom::Array<double, 2> m_surfaceNodeCoords;
  //!@brief Corners (index into m_surfaceNodeCoords) of generated surface cells.
  axom::Array<IndexType, 2> m_surfaceCellCorners;
  //!@brief Computational cell (flat index) crossing the surface cell.
  axom::Array<IndexType> m_surfaceCellParents;

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
    m_surfaceMesh->createField<int>(m_domainIdField, axom::mint::CELL_CENTERED);
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
        m_surfaceMesh->getFieldPtr<int>(m_domainIdField,
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
  SLIC_ASSERT_MSG(m_surfaceMesh,
                  "You must call set_output_mesh before compute_iso_surface.");
  SLIC_ASSERT_MSG(
    !m_fcnPath.empty(),
    "You must call set_function_field before compute_iso_surface.");

  const conduit::Node& coordValues =
    m_dom->fetch_existing(m_coordsetPath + "/values");
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  const int coordSp = isInterleaved ? m_ndim : 1;
  const double* xPtr = coordValues["x"].as_double_ptr();
  const double* yPtr = m_ndim >= 2 ? coordValues["y"].as_double_ptr() : nullptr;
  const double* zPtr = m_ndim >= 3 ? coordValues["z"].as_double_ptr() : nullptr;

  auto& fcnValues = m_dom->fetch_existing(m_fcnPath + "/values");
  const double* fcnPtr = fcnValues.as_double_ptr();

  const int* maskPtr = nullptr;
  if(!m_maskPath.empty())
  {
    auto& maskValues = m_dom->fetch_existing(m_maskPath + "/values");
    maskPtr = maskValues.as_int_ptr();
  }

  if(!m_cellIdField.empty() &&
     !m_surfaceMesh->hasField(m_cellIdField, axom::mint::CELL_CENTERED))
  {
    m_surfaceMesh->createField<int>(m_cellIdField, axom::mint::CELL_CENTERED);
  }

  _contourVal = contourVal;

  if(m_ndim == 2)
  {

#if 1
    MarchingCubesSingleDomainImpl<2, axom::execution_space<axom::SEQ_EXEC>> impl2d;
    impl2d.set_domain(*m_dom, m_coordsetPath, m_fcnPath, m_maskPath);
    impl2d.set_contour_value(contourVal);
    impl2d.mark_crossings();
    impl2d.scan_crossings();
    impl2d.generate_surface();
    using SegmentMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
    SegmentMesh* mesh = static_cast<SegmentMesh*>(m_surfaceMesh);
    for( int i=0; i<impl2d.m_surfaceNodeCoords.shape()[0]; ++i )
    {
      mesh->appendNode(impl2d.m_surfaceNodeCoords(i,0), impl2d.m_surfaceNodeCoords(i,1));
    }
    for( int i=0; i<impl2d.m_surfaceCellCorners.shape()[0]; ++i )
    {
      IndexType cell[2] {impl2d.m_surfaceCellCorners(i,0), impl2d.m_surfaceCellCorners(i,1)};
      mesh->appendCell(cell);
    }
#else
    /*
      For now, assume zero offsets and zero ghost width.
      Eventually, we'll have to support index offsets to
      handle data with ghosts.

      By using ArrayView, we're assuming row-major layout.  To support
      column-major layout as well, we have to extend ArrayView to
      compute column-major strides.  We should also use an iteration
      scheme that's cache-efficient for both layouts.
    */
    const axom::StackArray<axom::IndexType, 2> cShape {m_cShape[0], m_cShape[1]};
    const axom::StackArray<axom::IndexType, 2> pShape {1 + m_cShape[0],
                                                       1 + m_cShape[1]};
    axom::ArrayView<const double, 2> xView(xPtr, pShape, coordSp);
    axom::ArrayView<const double, 2> yView(yPtr, pShape, coordSp);
    axom::ArrayView<const double, 2> fcnView(fcnPtr, pShape);
    axom::ArrayView<const int, 2> maskView(maskPtr, cShape);

    // Write as regular nested loops.
    for(int i = 0; i < cShape[0]; ++i)
    {
      for(int j = 0; j < cShape[1]; ++j)
      {
        const bool skipZone = maskPtr && bool(maskView(i, j));
        if(!skipZone)
        {
          double nodalValues[4];
          double xx[4];
          double yy[4];

          nodalValues[0] = fcnView(i + 1, j);
          nodalValues[1] = fcnView(i + 1, j + 1);
          nodalValues[2] = fcnView(i, j + 1);
          nodalValues[3] = fcnView(i, j);

          xx[0] = xView(i + 1, j);
          xx[1] = xView(i + 1, j + 1);
          xx[2] = xView(i, j + 1);
          xx[3] = xView(i, j);

          yy[0] = yView(i + 1, j);
          yy[1] = yView(i + 1, j + 1);
          yy[2] = yView(i, j + 1);
          yy[3] = yView(i, j);

          auto nPrev = m_surfaceMesh->getNumberOfCells();
          this->contourCell2D(xx, yy, nodalValues);
          auto nNew = m_surfaceMesh->getNumberOfCells();

          if(nNew > nPrev && !m_cellIdField.empty())
          {
            int zoneIdx = i + j * cShape[0];  // TODO: Fix for ghost layer size.
            auto* cellIdPtr =
              m_surfaceMesh->getFieldPtr<int>(m_cellIdField,
                                              axom::mint::CELL_CENTERED);
            for(int n = nPrev; n < nNew; ++n)
            {
              cellIdPtr[n] = zoneIdx;
            }
          }
        }  // END if
      }
    }
#endif
  }
  else
  {
    const axom::StackArray<axom::IndexType, 3> cShape {m_cShape[0],
                                                       m_cShape[1],
                                                       m_cShape[2]};
    const axom::StackArray<axom::IndexType, 3> pShape {1 + m_cShape[0],
                                                       1 + m_cShape[1],
                                                       1 + m_cShape[2]};
    axom::ArrayView<const double, 3> xView(xPtr, pShape, coordSp);
    axom::ArrayView<const double, 3> yView(yPtr, pShape, coordSp);
    axom::ArrayView<const double, 3> zView(zPtr, pShape, coordSp);
    axom::ArrayView<const double, 3> fcnView(fcnPtr, pShape);
    axom::ArrayView<const int, 3> maskView(maskPtr, cShape);

    // Write as regular nested loops.
    for(int i = 0; i < cShape[0]; ++i)
    {
      for(int j = 0; j < cShape[1]; ++j)
      {
        for(int k = 0; k < cShape[2]; ++k)
        {
          const bool skipZone = maskPtr && bool(maskView(i, j, k));
          if(!skipZone)
          {
            double nodalValues[8];
            double xx[8];
            double yy[8];
            double zz[8];

            nodalValues[0] = fcnView(i + 1, j, k);
            nodalValues[1] = fcnView(i + 1, j + 1, k);
            nodalValues[2] = fcnView(i, j + 1, k);
            nodalValues[3] = fcnView(i, j, k);
            nodalValues[4] = fcnView(i + 1, j, k + 1);
            nodalValues[5] = fcnView(i + 1, j + 1, k + 1);
            nodalValues[6] = fcnView(i, j + 1, k + 1);
            nodalValues[7] = fcnView(i, j, k + 1);

            xx[0] = xView(i + 1, j, k);
            xx[1] = xView(i + 1, j + 1, k);
            xx[2] = xView(i, j + 1, k);
            xx[3] = xView(i, j, k);
            xx[4] = xView(i + 1, j, k + 1);
            xx[5] = xView(i + 1, j + 1, k + 1);
            xx[6] = xView(i, j + 1, k + 1);
            xx[7] = xView(i, j, k + 1);

            yy[0] = yView(i + 1, j, k);
            yy[1] = yView(i + 1, j + 1, k);
            yy[2] = yView(i, j + 1, k);
            yy[3] = yView(i, j, k);
            yy[4] = yView(i + 1, j, k + 1);
            yy[5] = yView(i + 1, j + 1, k + 1);
            yy[6] = yView(i, j + 1, k + 1);
            yy[7] = yView(i, j, k + 1);

            zz[0] = zView(i + 1, j, k);
            zz[1] = zView(i + 1, j + 1, k);
            zz[2] = zView(i, j + 1, k);
            zz[3] = zView(i, j, k);
            zz[4] = zView(i + 1, j, k + 1);
            zz[5] = zView(i + 1, j + 1, k + 1);
            zz[6] = zView(i, j + 1, k + 1);
            zz[7] = zView(i, j, k + 1);

            auto nPrev = m_surfaceMesh->getNumberOfCells();
            this->contourCell3D(xx, yy, zz, nodalValues);
            auto nNew = m_surfaceMesh->getNumberOfCells();

            if(nNew > nPrev && !m_cellIdField.empty())
            {
              int zoneIdx = i + j * cShape[0] +
                k * cShape[0] * cShape[1];  // TODO: Fix for ghost layer size.
              auto* cellIdPtr =
                m_surfaceMesh->getFieldPtr<int>(m_cellIdField,
                                                axom::mint::CELL_CENTERED);
              for(int n = nPrev; n < nNew; ++n)
              {
                cellIdPtr[n] = zoneIdx;
              }
            }
          }  // END if
        }
      }
    }
  }  // END else
}

//------------------------------------------------------------------------------
void MarchingCubesSingleDomain::linear_interp(int edgeIdx,
                                              const double* xx,
                                              const double* yy,
                                              const double* zz,
                                              const double* nodeValues,
                                              double* xyz)
{
  SLIC_ASSERT(xx != NULL);
  SLIC_ASSERT(yy != NULL);
  SLIC_ASSERT(nodeValues != NULL);

  // STEP 0: get the edge node indices
  // 2 nodes define the edge.  n1 and n2 are the indices of
  // the nodes w.r.t. the square or cubic zone.  There is a
  // agreed-on ordering of these indices in the arrays xx, yy,
  // zz, nodeValues, xyz.
  int n1 = edgeIdx;
  int n2 = (edgeIdx == 3) ? 0 : edgeIdx + 1;

  if(m_ndim == 3)
  {
    SLIC_ASSERT(zz != NULL);

    const int hex_edge_table[] = {
      0, 1, 1, 2, 2, 3, 3, 0,  // base
      4, 5, 5, 6, 6, 7, 7, 4,  // top
      0, 4, 1, 5, 2, 6, 3, 7   // vertical
    };

    n1 = hex_edge_table[edgeIdx * 2];
    n2 = hex_edge_table[edgeIdx * 2 + 1];

  }  // END if 3-D

  // STEP 1: get the fields and coordinates from the two points
  const double f1 = nodeValues[n1];
  const double f2 = nodeValues[n2];

  double p1[3];
  p1[0] = xx[n1];
  p1[1] = yy[n1];
  p1[2] = 0.0;

  double p2[3];
  p2[0] = xx[n2];
  p2[1] = yy[n2];
  p2[2] = 0.0;

  if(m_ndim == 3)
  {
    // set the z--coordinate if in 3-D
    p1[2] = zz[n1];
    p2[2] = zz[n2];
  }

  // STEP 2: check whether the interpolated point is at one of the two corners.
  if(axom::utilities::isNearlyEqual(_contourVal, f1) ||
     axom::utilities::isNearlyEqual(f1, f2))
  {
    memcpy(xyz, p1, m_ndim * sizeof(double));
    return;
  }

  if(axom::utilities::isNearlyEqual(_contourVal, f2))
  {
    memcpy(xyz, p2, m_ndim * sizeof(double));
    return;
  }

  // STEP 3: point is in between the edge points, interpolate its position
  constexpr double ptiny = 1.0e-80;
  const double df = f2 - f1 + ptiny;  //add ptiny to avoid division by zero
  const double w = (_contourVal - f1) / df;
  for(int i = 0; i < m_ndim; ++i)
  {
    xyz[i] = p1[i] + w * (p2[i] - p1[i]);
  }
}

//------------------------------------------------------------------------------
int MarchingCubesSingleDomain::computeIndex(const double* f)
{
  const int numNodes = (m_ndim == 3) ? 8 : 4;

  int index = 0;
  for(int i = 0; i < numNodes; ++i)
  {
    if(f[i] >= _contourVal)
    {
      const int mask = (1 << i);
      index |= mask;
    }
  }
  return (index);
}

//------------------------------------------------------------------------------
void MarchingCubesSingleDomain::contourCell2D(double xx[4],
                                              double yy[4],
                                              double nodeValues[4])
{
  SLIC_ASSERT(xx != NULL);
  SLIC_ASSERT(yy != NULL);
  SLIC_ASSERT(nodeValues != NULL);

  // compute index
  int index = MarchingCubesSingleDomain::computeIndex(nodeValues);
  SLIC_ASSERT((index >= 0) && (index < 16));

  // short-circuit
  if(detail::num_segments[index] == 0)
  {
    return;
  }  // END if

  // Generate line segments
  using SegmentMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
  SegmentMesh* mesh = static_cast<SegmentMesh*>(m_surfaceMesh);
  SLIC_ASSERT(mesh != NULL);
  SLIC_ASSERT(mesh->getCellType() == axom::mint::SEGMENT);
  SLIC_ASSERT(mesh->getDimension() == 2);

  IndexType idx = mesh->getNumberOfNodes();
  IndexType cell[2];
  double p[2];

  const int nsegs = detail::num_segments[index];

  for(int i = 0; i < nsegs; ++i)
  {
    const int e1 = detail::cases2D[index][i * 2];
    const int e2 = detail::cases2D[index][i * 2 + 1];

    MarchingCubesSingleDomain::linear_interp(e1, xx, yy, NULL, nodeValues, p);
    mesh->appendNode(p[0], p[1]);
    cell[0] = idx;
    ++idx;

    MarchingCubesSingleDomain::linear_interp(e2, xx, yy, NULL, nodeValues, p);
    mesh->appendNode(p[0], p[1]);
    cell[1] = idx;
    ++idx;

    mesh->appendCell(cell);

  }  // END for all segments
}

//------------------------------------------------------------------------------
void MarchingCubesSingleDomain::contourCell3D(double xx[8],
                                              double yy[8],
                                              double zz[8],
                                              double nodeValues[8])
{
  SLIC_ASSERT(xx != NULL);
  SLIC_ASSERT(yy != NULL);
  SLIC_ASSERT(zz != NULL);
  SLIC_ASSERT(nodeValues != NULL);

  // compute index
  int index = MarchingCubesSingleDomain::computeIndex(nodeValues);
  SLIC_ASSERT((index >= 0) && (index < 256));

  // short-circuit
  if(detail::num_triangles[index] == 0)
  {
    return;
  }  // END if

  // Generate triangles
  using TriangleMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
  TriangleMesh* mesh = static_cast<TriangleMesh*>(m_surfaceMesh);
  SLIC_ASSERT(mesh != NULL);
  SLIC_ASSERT(mesh->getCellType() == axom::mint::TRIANGLE);
  SLIC_ASSERT(mesh->getDimension() == 3);

  const int numTriangles = detail::num_triangles[index];

  IndexType idx = mesh->getNumberOfNodes();
  IndexType cell[3];
  double p[3];

  for(int i = 0; i < numTriangles; ++i)
  {
    const int e1 = detail::cases3D[index][i * 3];
    const int e2 = detail::cases3D[index][i * 3 + 1];
    const int e3 = detail::cases3D[index][i * 3 + 2];

    MarchingCubesSingleDomain::linear_interp(e1, xx, yy, zz, nodeValues, p);
    mesh->appendNode(p[0], p[1], p[2]);
    cell[0] = idx;
    ++idx;

    MarchingCubesSingleDomain::linear_interp(e2, xx, yy, zz, nodeValues, p);
    mesh->appendNode(p[0], p[1], p[2]);
    cell[1] = idx;
    ++idx;

    MarchingCubesSingleDomain::linear_interp(e3, xx, yy, zz, nodeValues, p);
    mesh->appendNode(p[0], p[1], p[2]);
    cell[2] = idx;
    ++idx;

    mesh->appendCell(cell);

  }  // END for all triangles
}

}  // end namespace quest
}  // end namespace axom
