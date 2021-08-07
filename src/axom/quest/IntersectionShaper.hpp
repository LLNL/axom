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

  void runShapeQuery() override
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

    for(int i = 0; i < NE; i++)
    {
      // Get the indices of this element's vertices
      mfem::Array<int> verts;
      mesh->GetElementVertices(i, verts);
      SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_HEX);

      // Get the coordinates for the vertices
      double* vertCoords[NUM_VERTS_PER_HEX];
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        vertCoords[j] = mesh->GetVertex(verts[j]);
      }

      // Do something with the vertex coordinates; in this case, just print them
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        SLIC_INFO(
          fmt::format("Element {} -- coords for vertex {} are: {} {} {}",
                      i,
                      j,
                      vertCoords[j][0],
                      vertCoords[j][1],
                      vertCoords[j][2]));
      }
    }
  }

  void applyReplacementRules(const klee::Shape& shape) override
  {
    const auto& shapeName = shape.getName();
    SLIC_INFO(fmt::format(
      "{:-^80}",
      fmt::format("Applying replacement rules over for shape '{}'", shapeName)));

    // Implementation here -- update volume fractions based on replacement rules
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
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
