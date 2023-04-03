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
