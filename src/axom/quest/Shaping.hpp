// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Shaping.hpp
 *
 * \brief Helper class for shaping queries
 */

#ifndef AXOM_QUEST_SHAPING__HPP_
#define AXOM_QUEST_SHAPING__HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"

#include "axom/quest/InOutOctree.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"

#include "mfem.hpp"

#include "fmt/fmt.hpp"
#include "fmt/locale.h"

namespace axom
{
namespace quest
{
namespace shaping
{
using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

using TriVertIndices = primal::Point<axom::IndexType, 3>;
using SpaceTriangle = primal::Triangle<double, 3>;

using Octree3D = quest::InOutOctree<3>;

using GeometricBoundingBox = Octree3D::GeometricBoundingBox;
using SpacePt = Octree3D::SpacePt;
using SpaceVector = Octree3D::SpaceVector;
using GridPt = Octree3D::GridPt;
using BlockIndex = Octree3D::BlockIndex;

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

/** Computes the bounding box of the surface mesh */
GeometricBoundingBox compute_bounds(const mint::Mesh& mesh)
{
  GeometricBoundingBox meshBB;
  SpacePt pt;

  for(int i = 0; i < mesh.getNumberOfNodes(); ++i)
  {
    mesh.getNode(i, pt.data());
    meshBB.addPoint(pt);
  }

  SLIC_ASSERT(meshBB.isValid());

  return meshBB;
}

}  // end namespace shaping

template <int NDIMS>
class SamplingShaper
{
public:
  static constexpr int DIM = NDIMS;
  using InOutOctreeType = quest::InOutOctree<DIM>;

  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpacePt = typename InOutOctreeType::SpacePt;
  using SpaceVector = typename InOutOctreeType::SpaceVector;
  using GridPt = typename InOutOctreeType::GridPt;
  using BlockIndex = typename InOutOctreeType::BlockIndex;

public:
  SamplingShaper(const std::string& shapeName, mint::Mesh* surfaceMesh)
    : m_shapeName(shapeName)
    , m_surfaceMesh(surfaceMesh)
  { }

  ~SamplingShaper() { delete m_octree; }

  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    m_bbox.clear();
    SpacePt pt;

    for(int i = 0; i < m_surfaceMesh->getNumberOfNodes(); ++i)
    {
      m_surfaceMesh->getNode(i, pt.data());
      m_bbox.addPoint(pt);
    }

    SLIC_ASSERT(m_bbox.isValid());

    SLIC_INFO("Mesh bounding box: " << m_bbox);
  }

  void initSpatialIndex()
  {
    // Create octree over mesh's bounding box
    m_octree = new InOutOctreeType(m_bbox, m_surfaceMesh);
    m_octree->generateIndex();
  }

  void sampleInOutField(mfem::DataCollection* dc,
                        shaping::QFunctionCollection& inoutQFuncs,
                        int sampleRes)
  {
    auto* mesh = dc->GetMesh();
    SLIC_ASSERT(mesh != nullptr);
    const int NE = mesh->GetNE();
    const int dim = mesh->Dimension();

    // Generate a Quadrature Function with the geometric positions, if not already available
    if(!inoutQFuncs.Has("positions"))
    {
      shaping::generatePositionsQFunction(mesh, inoutQFuncs, sampleRes);
    }

    // Access the positions QFunc and associated QuadratureSpace
    mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
    mfem::QuadratureSpace* sp = pos_coef->GetSpace();
    const int nq = sp->GetElementIntRule(0).GetNPoints();

    // Sample the in/out field at each point
    // store in QField which we register with the QFunc collection
    const std::string inoutName = fmt::format("inout_{}", m_shapeName);
    const int vdim = 1;
    auto* inout = new mfem::QuadratureFunction(sp, vdim);
    inoutQFuncs.Register(inoutName, inout, true);

    mfem::DenseMatrix m;
    mfem::Vector res;

    axom::utilities::Timer timer(true);
    for(int i = 0; i < NE; ++i)
    {
      pos_coef->GetElementValues(i, m);
      inout->GetElementValues(i, res);

      for(int p = 0; p < nq; ++p)
      {
        const SpacePt pt(m.GetColumn(p), dim);
        const bool in = m_octree->within(pt);
        res(p) = in ? 1. : 0.;

        // SLIC_INFO(fmt::format("[{},{}] Pt: {}, In: {}", i,p,pt, (in? "yes" : "no") ));
      }
    }
    timer.stop();

    SLIC_INFO(fmt::format(std::locale("en_US.UTF-8"),
                          "\t Sampling inout field '{}' took {} seconds (@ "
                          "{:L} queries per second)",
                          inoutName,
                          timer.elapsed(),
                          static_cast<int>((NE * nq) / timer.elapsed())));
  }

  /**
  * Compute volume fractions function for shape on a grid of resolution \a gridRes
  * in region defined by bounding box \a queryBounds
  */
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc,
                                      int AXOM_NOT_USED(sampleRes),
                                      int outputOrder)
  {
    // Step 1 -- generate a QField w/ the spatial coordinates
    mfem::Mesh* mesh = dc->GetMesh();
    const int NE = mesh->GetNE();
    const int dim = mesh->Dimension();

    if(NE < 1)
    {
      SLIC_WARNING("Mesh has no elements!");
      return;
    }

    mfem::L2_FECollection* coll =
      new mfem::L2_FECollection(outputOrder, dim, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
    mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(coll);
    auto volFracName = fmt::format("vol_frac_{}", m_shapeName);
    dc->RegisterField(volFracName, volFrac);

    auto* fe = fes->GetFE(0);
    auto& ir = fe->GetNodes();

    // Assume all elements have the same integration rule
    const int nq = ir.GetNPoints();
    const auto* geomFactors =
      mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

    mfem::DenseTensor pos_coef(dim, nq, NE);

    // Rearrange positions into quadrature function
    {
      for(int i = 0; i < NE; ++i)
      {
        for(int j = 0; j < dim; ++j)
        {
          for(int k = 0; k < nq; ++k)
          {
            pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
          }
        }
      }
    }

    // Step 2 -- sample the in/out field at each point -- store directly in volFrac grid function
    mfem::Vector res(nq);
    mfem::Array<int> dofs;
    for(int i = 0; i < NE; ++i)
    {
      mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const SpacePt pt(m.GetColumn(p), dim);
        const bool in = m_octree->within(pt);
        res(p) = in ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }

private:
  std::string m_shapeName;

  GeometricBoundingBox m_bbox;
  mint::Mesh* m_surfaceMesh {nullptr};
  InOutOctreeType* m_octree {nullptr};
};

/// Helper class to use an MFEM sampling-based shaper
class MFEMShaping
{
public:
  MFEMShaping(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
    : m_shapeSet(shapeSet)
    , m_dc(dc)
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    m_comm = m_dc->GetComm();
#endif
  }

  void setSamplesPerKnotSpan(int nSamples)
  {
    using axom::utilities::clampLower;
    SLIC_WARNING_IF(
      nSamples < 1,
      fmt::format(
        "Samples per knot span must be at least 1. Provided value was {}",
        nSamples));

    m_samplesPerKnotSpan = clampLower(nSamples, 1);
  }

  void setVertexWeldThreshold(double threshold)
  {
    SLIC_WARNING_IF(
      threshold <= 0.,
      fmt::format(
        "Vertex weld threshold should be positive Provided value was {}",
        threshold));

    m_vertexWeldThreshold = threshold;
  }

  sidre::MFEMSidreDataCollection* getDC() { return m_dc; }
  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  /// Loads the shape from file into m_surfaceMesh
  void loadShape(const klee::Shape& shape)
  {
    using axom::utilities::string::endsWith;
    SLIC_ASSERT(shape.getGeometry().getFormat() == "stl" ||
                shape.getGeometry().getFormat() == "c2c");

    const std::string shapeName = shape.getName();
    std::string outMsg = fmt::format(" Loading shape '{}' ", shapeName);
    SLIC_INFO(fmt::format("{:-^80}", outMsg));

    std::string shapePath = m_shapeSet.resolvePath(shape.getGeometry().getPath());
    SLIC_INFO("Reading file: " << shapePath << "...");

    if(endsWith(shapePath, ".stl"))
    {
      quest::internal::read_stl_mesh(shapePath, m_surfaceMesh, m_comm);
    }
#ifdef AXOM_USE_C2C
    else if(endsWith(shapePath, ".contour"))
    {
      quest::internal::read_c2c_mesh(shapePath,
                                     m_samplesPerKnotSpan,
                                     m_vertexWeldThreshold,
                                     m_surfaceMesh,
                                     m_comm);
    }
#endif
    else
    {
      SLIC_ERROR(
        fmt::format("Unsupported filetype for this Axom configuration. "
                    "Provided file was '{}'",
                    shapePath));
    }
  }

  void applyTransforms(const klee::Shape& shape)
  {
    // TODO: Implement this as a set of affine transforms to vertices of mesh
    AXOM_UNUSED_VAR(shape);
  }

  void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape)
  {
    const auto& shapeName = shape.getName();

    SLIC_INFO(fmt::format("{:-^80}", " Generating the octree "));
    switch(shapeDimension)
    {
    case klee::Dimensions::Two:
      m_samplingShaper2D = new SamplingShaper<2>(shapeName, m_surfaceMesh);
      m_samplingShaper2D->computeBounds();
      m_samplingShaper2D->initSpatialIndex();
      m_surfaceMesh = m_samplingShaper2D->getSurfaceMesh();
      break;

    case klee::Dimensions::Three:
      m_samplingShaper3D = new SamplingShaper<3>(shapeName, m_surfaceMesh);
      m_samplingShaper3D->computeBounds();
      m_samplingShaper3D->initSpatialIndex();
      m_surfaceMesh = m_samplingShaper3D->getSurfaceMesh();
      break;

    default:
      SLIC_ERROR(
        "Shaping dimension must be 2 or 3, but requested dimension was "
        << static_cast<int>(shapeDimension));
      break;
    }

    // Check that one of sampling shapers (2D or 3D) is null and the other is not
    SLIC_ASSERT((m_samplingShaper2D == nullptr && m_samplingShaper3D != nullptr) ||
                (m_samplingShaper3D == nullptr && m_samplingShaper2D != nullptr));

    // Output some logging info and dump the mesh
    {
      const int nVerts = m_surfaceMesh->getNumberOfNodes();
      const int nCells = m_surfaceMesh->getNumberOfCells();

      SLIC_INFO(fmt::format(
        "After welding, surface mesh has {} vertices  and {} triangles.",
        nVerts,
        nCells));
      mint::write_vtk(m_surfaceMesh,
                      fmt::format("meldedTriMesh_{}.vtk", shapeName));
    }
  }

  // Handles 2D or 3D shaping, based on the template and associated parameter
  template <typename SamplingShaperType>
  void runShapeQueryImpl(SamplingShaperType* shaper,
                         shaping::VolFracSampling vfSampling,
                         int samplingOrder,
                         int outputOrder)
  {
    // Sample the InOut field at the mesh quadrature points
    switch(vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      shaper->sampleInOutField(m_dc, m_inoutShapeQFuncs, samplingOrder);
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      shaper->computeVolumeFractionsBaseline(m_dc, samplingOrder, outputOrder);
      break;
    }
  }

  void runShapeQuery(shaping::VolFracSampling vfSampling,
                     int samplingOrder,
                     int outputOrder)
  {
    SLIC_INFO(fmt::format("{:-^80}", " Querying the octree "));

    switch(getShapeDimension())
    {
    case klee::Dimensions::Two:
      runShapeQueryImpl(m_samplingShaper2D, vfSampling, samplingOrder, outputOrder);
      break;
    case klee::Dimensions::Three:
      runShapeQueryImpl(m_samplingShaper3D, vfSampling, samplingOrder, outputOrder);
      break;
    }
  }

  void applyReplacementRules(const klee::Shape& shape)
  {
    using axom::utilities::string::splitLastNTokens;

    SLIC_INFO(
      fmt::format("{:=^80}", "Applying replacement rules over the shapes"));

    // Get inout qfunc for shape
    const auto& shapeName = shape.getName();
    auto* shapeQFunc = m_inoutShapeQFuncs.Get(fmt::format("inout_{}", shapeName));
    SLIC_ASSERT_MSG(
      shapeQFunc != nullptr,
      fmt::format("Missing inout samples for shape '{}'", shapeName));

    const auto& materialName = shape.getMaterial();

    // Creata a copy of the inout samples for this shape
    // Replacements will be applied to this and then copied into our shape
    auto* shapeQFuncCopy = new mfem::QuadratureFunction(*shapeQFunc);

    // apply replacement rules to all other materials
    for(auto& mat : m_inoutMaterialQFuncs)
    {
      const std::string otherMatName = splitLastNTokens(mat.first, 2, '_')[1];

      // We'll handle the current shape's material at the end
      if(otherMatName == materialName)
      {
        continue;
      }

      const bool shouldReplace = shape.replaces(otherMatName);
      SLIC_INFO(fmt::format(
        "Should we replace material '{}' with shape '{}' of material '{}'? {}",
        otherMatName,
        shapeName,
        materialName,
        shouldReplace ? "yes" : "no"));

      const std::string& materialQFuncName = mat.first;
      auto* otherMatQFunc = m_inoutMaterialQFuncs.Get(materialQFuncName);
      SLIC_ASSERT_MSG(
        otherMatQFunc != nullptr,
        fmt::format("Missing inout samples for material '{}'", otherMatName));

      quest::shaping::replaceMaterial(shapeQFuncCopy, otherMatQFunc, shouldReplace);
    }

    // Get inout qfunc for material
    const std::string materialQFuncName =
      fmt::format("mat_inout_{}", materialName);
    if(!m_inoutMaterialQFuncs.Has(materialQFuncName))
    {
      // initialize material from shape inout
      auto* mat_inout = new mfem::QuadratureFunction(*shapeQFuncCopy);
      m_inoutMaterialQFuncs.Register(materialQFuncName, mat_inout, true);
    }
    else
    {
      auto* matQFunc = m_inoutMaterialQFuncs.Get(materialQFuncName);
      SLIC_ASSERT_MSG(
        matQFunc != nullptr,
        fmt::format("Missing inout samples for material '{}'", materialName));

      quest::shaping::copyShapeIntoMaterial(shapeQFuncCopy, matQFunc);
    }
  }

  void finalizeShapeQuery()
  {
    delete m_samplingShaper2D;
    m_samplingShaper2D = nullptr;

    delete m_samplingShaper3D;
    m_samplingShaper3D = nullptr;

    delete m_surfaceMesh;
    m_surfaceMesh = nullptr;
  }

  void adjustVolumeFractions(shaping::VolFracSampling vfSampling, int outputOrder)
  {
    for(auto& mat : m_inoutMaterialQFuncs)
    {
      const std::string matName = mat.first;
      SLIC_INFO(
        fmt::format("Generating volume fraction fields for '{}' material",
                    matName));

      // Sample the InOut field at the mesh quadrature points
      switch(vfSampling)
      {
      case shaping::VolFracSampling::SAMPLE_AT_QPTS:
        quest::shaping::computeVolumeFractions(matName,
                                               m_dc,
                                               m_inoutMaterialQFuncs,
                                               outputOrder);
        break;
      case shaping::VolFracSampling::SAMPLE_AT_DOFS:
        /* no-op for now */
        break;
      }
    }
  }

private:
  klee::Dimensions getShapeDimension() const
  {
    const bool has2D = (m_samplingShaper2D != nullptr);
    const bool has3D = (m_samplingShaper3D != nullptr);
    SLIC_ERROR_IF(!(has2D || has3D), "Shape not initialized");
    SLIC_ERROR_IF(has2D && has3D, "Cannot have concurrent 2D and 3D shapes");

    return has2D ? klee::Dimensions::Two : klee::Dimensions::Three;
  }

public:
  const klee::ShapeSet& m_shapeSet;
  sidre::MFEMSidreDataCollection* m_dc;

  // TODO: Use MfemSidreDataCollection QFuncs for this when we upgrade to post mfem@4.3
  shaping::QFunctionCollection m_inoutShapeQFuncs;
  shaping::QFunctionCollection m_inoutMaterialQFuncs;
  shaping::DenseTensorCollection m_inoutDofs;

  mint::Mesh* m_surfaceMesh {nullptr};
  SamplingShaper<2>* m_samplingShaper2D {nullptr};
  SamplingShaper<3>* m_samplingShaper3D {nullptr};

  int m_samplesPerKnotSpan {25};
  double m_vertexWeldThreshold {1e-9};

  MPI_Comm m_comm {MPI_COMM_SELF};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING__HPP_
