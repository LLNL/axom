// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SamplingShaper.hpp
 *
 * \brief Helper class for sampling-based shaping queries
 */

#ifndef AXOM_QUEST_SAMPLING_SHAPER__HPP_
#define AXOM_QUEST_SAMPLING_SHAPER__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "axom/quest/Shaper.hpp"
#include "axom/quest/InOutOctree.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
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
using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

template <int NDIMS>
class InOutSampler
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
  /**
   * \brief Constructor for a SamplingShaper
   *
   * \param shapeName The name of the shape; will be used for the field for the associated samples
   * \param surfaceMesh Pointer to the surface mesh
   *
   * \note Does not take ownership of the surface mesh
   */
  InOutSampler(const std::string& shapeName, mint::Mesh* surfaceMesh)
    : m_shapeName(shapeName)
    , m_surfaceMesh(surfaceMesh)
  { }

  ~InOutSampler() { delete m_octree; }

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
  DISABLE_COPY_AND_ASSIGNMENT(InOutSampler);
  DISABLE_MOVE_AND_ASSIGNMENT(InOutSampler);

  std::string m_shapeName;

  GeometricBoundingBox m_bbox;
  mint::Mesh* m_surfaceMesh {nullptr};
  InOutOctreeType* m_octree {nullptr};
};

}  // end namespace shaping

/*!
 * \brief Concrete class for sample based shaping
 */
class SamplingShaper : public Shaper
{
public:
  SamplingShaper(const klee::ShapeSet& shapeSet,
                 sidre::MFEMSidreDataCollection* dc)
    : Shaper(shapeSet, dc)
  { }

  //@{
  //!  @name Functions to get and set shaping parameters related to sampling; supplements parameters in base class

  void setSamplingType(shaping::VolFracSampling vfSampling)
  {
    m_vfSampling = vfSampling;
  }

  void setQuadratureOrder(int quadratureOrder)
  {
    m_quadratureOrder = quadratureOrder;
  }

  void setVolumeFractionOrder(int volfracOrder)
  {
    m_volfracOrder = volfracOrder;
  }

  //@}

private:
  klee::Dimensions getShapeDimension() const
  {
    const bool has2D = (m_inoutSampler2D != nullptr);
    const bool has3D = (m_inoutSampler3D != nullptr);
    SLIC_ERROR_IF(!(has2D || has3D), "Shape not initialized");
    SLIC_ERROR_IF(has2D && has3D, "Cannot have concurrent 2D and 3D shapes");

    return has2D ? klee::Dimensions::Two : klee::Dimensions::Three;
  }

public:
  //@{
  //!  @name Functions related to the stages for a given shape

  /// Initializes the spatial index for shaping
  void prepareShapeQuery(klee::Dimensions shapeDimension,
                         const klee::Shape& shape) override
  {
    const auto& shapeName = shape.getName();

    SLIC_INFO(fmt::format("{:-^80}", " Generating the octree "));
    switch(shapeDimension)
    {
    case klee::Dimensions::Two:
      m_inoutSampler2D = new shaping::InOutSampler<2>(shapeName, m_surfaceMesh);
      m_inoutSampler2D->computeBounds();
      m_inoutSampler2D->initSpatialIndex();
      m_surfaceMesh = m_inoutSampler2D->getSurfaceMesh();
      break;

    case klee::Dimensions::Three:
      m_inoutSampler3D = new shaping::InOutSampler<3>(shapeName, m_surfaceMesh);
      m_inoutSampler3D->computeBounds();
      m_inoutSampler3D->initSpatialIndex();
      m_surfaceMesh = m_inoutSampler3D->getSurfaceMesh();
      break;

    default:
      SLIC_ERROR(
        "Shaping dimension must be 2 or 3, but requested dimension was "
        << static_cast<int>(shapeDimension));
      break;
    }

    // Check that one of sampling shapers (2D or 3D) is null and the other is not
    SLIC_ASSERT((m_inoutSampler2D == nullptr && m_inoutSampler3D != nullptr) ||
                (m_inoutSampler3D == nullptr && m_inoutSampler2D != nullptr));

    // Output some logging info and dump the mesh
    if(this->getRank() == 0)
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

  void runShapeQuery(const klee::Shape& shape) override
  {
    SLIC_INFO(fmt::format(
      "{:-^80}",
      fmt::format(" Querying the octree for shape '{}'", shape.getName())));

    switch(getShapeDimension())
    {
    case klee::Dimensions::Two:
      runShapeQueryImpl(m_inoutSampler2D);
      break;
    case klee::Dimensions::Three:
      runShapeQueryImpl(m_inoutSampler3D);
      break;
    }
  }

  void applyReplacementRules(const klee::Shape& shape) override
  {
    using axom::utilities::string::splitLastNTokens;

    const auto& shapeName = shape.getName();
    SLIC_INFO(fmt::format(
      "{:-^80}",
      fmt::format("Applying replacement rules over for shape '{}'", shapeName)));

    // Get inout qfunc for this shape
    auto* shapeQFunc = m_inoutShapeQFuncs.Get(fmt::format("inout_{}", shapeName));
    SLIC_ASSERT_MSG(
      shapeQFunc != nullptr,
      fmt::format("Missing inout samples for shape '{}'", shapeName));

    // Create a copy of the inout samples for this shape
    // Replacements will be applied to this and then copied into our shape's material
    auto* shapeQFuncCopy = new mfem::QuadratureFunction(*shapeQFunc);

    // apply replacement rules to all other materials
    const auto& thisMatName = shape.getMaterial();
    for(auto& mat : m_inoutMaterialQFuncs)
    {
      const std::string otherMatName = splitLastNTokens(mat.first, 2, '_')[1];

      // We'll handle the current shape's material at the end
      if(otherMatName == thisMatName)
      {
        continue;
      }

      const bool shouldReplace = shape.replaces(otherMatName);
      SLIC_DEBUG(fmt::format(
        "Should we replace material '{}' with shape '{}' of material '{}'? {}",
        otherMatName,
        shapeName,
        thisMatName,
        shouldReplace ? "yes" : "no"));

      auto* otherMatQFunc = mat.second;
      SLIC_ASSERT_MSG(
        otherMatQFunc != nullptr,
        fmt::format("Missing inout samples for material '{}'", otherMatName));

      quest::shaping::replaceMaterial(shapeQFuncCopy, otherMatQFunc, shouldReplace);
    }

    // Get inout qfunc for the current material
    const std::string materialQFuncName =
      fmt::format("mat_inout_{}", thisMatName);
    if(!m_inoutMaterialQFuncs.Has(materialQFuncName))
    {
      // initialize material from shape inout, the QFunc registry takes ownership
      m_inoutMaterialQFuncs.Register(materialQFuncName, shapeQFuncCopy, true);
    }
    else
    {
      // copy shape data into current material and delete the copy
      auto* matQFunc = m_inoutMaterialQFuncs.Get(materialQFuncName);
      SLIC_ASSERT_MSG(
        matQFunc != nullptr,
        fmt::format("Missing inout samples for material '{}'", thisMatName));

      quest::shaping::copyShapeIntoMaterial(shapeQFuncCopy, matQFunc);

      delete shapeQFuncCopy;
      shapeQFuncCopy = nullptr;
    }
  }

  void finalizeShapeQuery() override
  {
    delete m_inoutSampler2D;
    m_inoutSampler2D = nullptr;

    delete m_inoutSampler3D;
    m_inoutSampler3D = nullptr;

    delete m_surfaceMesh;
    m_surfaceMesh = nullptr;
  }

  //@}

public:
  void adjustVolumeFractions() override
  {
    for(auto& mat : m_inoutMaterialQFuncs)
    {
      const std::string matName = mat.first;
      SLIC_INFO(
        fmt::format("Generating volume fraction fields for '{}' material",
                    matName));

      // Sample the InOut field at the mesh quadrature points
      switch(m_vfSampling)
      {
      case shaping::VolFracSampling::SAMPLE_AT_QPTS:
        quest::shaping::computeVolumeFractions(matName,
                                               m_dc,
                                               m_inoutMaterialQFuncs,
                                               m_volfracOrder);
        break;
      case shaping::VolFracSampling::SAMPLE_AT_DOFS:
        /* no-op for now */
        break;
      }
    }
  }

private:
  // Handles 2D or 3D shaping, based on the template and associated parameter
  template <typename InOutSamplerType>
  void runShapeQueryImpl(InOutSamplerType* shaper)
  {
    // Sample the InOut field at the mesh quadrature points
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      shaper->sampleInOutField(m_dc, m_inoutShapeQFuncs, m_quadratureOrder);
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      shaper->computeVolumeFractionsBaseline(m_dc,
                                             m_quadratureOrder,
                                             m_volfracOrder);
      break;
    }
  }

private:
  // TODO: Use MfemSidreDataCollection QFuncs for this when we upgrade to post mfem@4.3
  shaping::QFunctionCollection m_inoutShapeQFuncs;
  shaping::QFunctionCollection m_inoutMaterialQFuncs;
  shaping::DenseTensorCollection m_inoutDofs;

  shaping::InOutSampler<2>* m_inoutSampler2D {nullptr};
  shaping::InOutSampler<3>* m_inoutSampler3D {nullptr};

  shaping::VolFracSampling m_vfSampling {shaping::VolFracSampling::SAMPLE_AT_QPTS};
  int m_quadratureOrder {5};
  int m_volfracOrder {2};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SAMPLING_SHAPER__HPP_
