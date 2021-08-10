// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IntersectionShaper.hpp
 *
 * \brief Helper class for intersection-based shaping queries
 */

#ifndef AXOM_QUEST_INTERSECTION_SHAPER__HPP_
#define AXOM_QUEST_INTERSECTION_SHAPER__HPP_

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
class IntersectionShaper : public Shaper
{
public:
  IntersectionShaper(const klee::ShapeSet& shapeSet,
                     sidre::MFEMSidreDataCollection* dc)
    : Shaper(shapeSet, dc)
  { }

public:
  //@{
  //!  @name Functions related to the stages for a given shape

  /// Initializes the spatial index for shaping
  void prepareShapeQuery(klee::Dimensions shapeDimension,
                         const klee::Shape& shape) override
  {
    const auto& shapeName = shape.getName();
    AXOM_UNUSED_VAR(shapeDimension);
    AXOM_UNUSED_VAR(shapeName);

    // Implementation here
    //  -- generate the octahedra (probably in a previous step)
    //  -- generate the BVH tree over the octahedra
  }

  void runShapeQuery(const klee::Shape& shape) override
  {
    constexpr int NUM_VERTS_PER_HEX = 8;

    SLIC_INFO(fmt::format("{:-^80}", " Querying the BVH tree "));

    mfem::Mesh* mesh = getDC()->GetMesh();

    // Intersection algorithm only works on linear elements
    SLIC_ASSERT(mesh != nullptr);
    int const NE = mesh->GetNE();
    if(NE > 0)
    {
      SLIC_ASSERT(mesh->GetNodes() == nullptr ||
                  mesh->GetNodes()->FESpace()->GetOrder(0));
    }

    // Create and register a scalar field for this shape's volume fractions
    // The Degrees of Freedom will be in correspondence with the elements
    auto* volFrac = this->newVolFracGridFunction();
    auto volFracName = fmt::format("shape_vol_frac_{}", shape.getName());
    this->getDC()->RegisterField(volFracName, volFrac);

    for(int el = 0; el < NE; ++el)
    {
      // Get the indices of this element's vertices
      mfem::Array<int> verts;
      mesh->GetElementVertices(el, verts);
      SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_HEX);

      // Get the coordinates for the vertices
      double* vertCoords[NUM_VERTS_PER_HEX];
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        vertCoords[j] = mesh->GetVertex(verts[j]);
      }

      // Do something with the vertex coordinates; in this case, just print them
      if(this->isVerbose())
      {
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          SLIC_INFO(
            fmt::format("Element {} -- coords for vertex {} are: {} {} {}",
                        el,
                        j,
                        vertCoords[j][0],
                        vertCoords[j][1],
                        vertCoords[j][2]));
        }
      }

      // run the query for this element
      // HACK: To get a concrete value, let's just set every element's vf to 1 for now
      double vf = 1.0;
      (*volFrac)(el) = vf;
    }
  }

  void applyReplacementRules(const klee::Shape& shape) override
  {
    const auto& shapeName = shape.getName();
    const auto& materialName = shape.getMaterial();
    SLIC_INFO(fmt::format(
      "{:-^80}",
      fmt::format("Applying replacement rules for shape '{}' of material {}",
                  shapeName,
                  materialName)));

    auto shapeVolFracName = fmt::format("shape_vol_frac_{}", shapeName);
    auto materialVolFracName = fmt::format("vol_frac_{}", materialName);

    auto* shapeVolFrac = this->getDC()->GetField(shapeVolFracName);
    SLIC_ASSERT(shapeVolFrac != nullptr);

    // Get or create the volume fraction field for this shape's material
    mfem::GridFunction* matVolFrac = nullptr;
    if(this->getDC()->HasField(materialVolFracName))
    {
      matVolFrac = this->getDC()->GetField(materialVolFracName);
    }
    else
    {
      matVolFrac = newVolFracGridFunction();
      this->getDC()->RegisterField(materialVolFracName, matVolFrac);
    }

    /// Implementation here -- update material volume fractions based on replacement rules

    // HACK: For simplicity at first, let's just set the material vol fraction to that of the shape
    // This will have to be updated
    *matVolFrac = *shapeVolFrac;
  }

  void finalizeShapeQuery() override
  {
    // Implementation here -- destroy BVH tree and other shape-based data structures
    delete m_surfaceMesh;
    m_surfaceMesh = nullptr;
  }

  //@}

public:
  void adjustVolumeFractions() override
  {
    // Implementation here -- not sure if this will require anything for intersection-based shaping
  }

private:
  /// Create and return a new volume fraction grid function for the current mesh
  mfem::GridFunction* newVolFracGridFunction()
  {
    mfem::Mesh* mesh = getDC()->GetMesh();
    SLIC_ASSERT(mesh != nullptr);

    const int vfOrder = 0;
    const int dim = mesh->Dimension();
    mfem::L2_FECollection* coll =
      new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
    mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(coll);

    return volFrac;
  }
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
