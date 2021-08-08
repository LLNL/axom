// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file Shaper.hpp
 *
 * \brief Helper class for shaping queries
 */

#ifndef AXOM_QUEST_SHAPER__HPP_
#define AXOM_QUEST_SHAPER__HPP_

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
/**
 * Abstract base class for shaping material volume fractions
 */
class Shaper
{
public:
  Shaper(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
    : m_shapeSet(shapeSet)
    , m_dc(dc)
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    m_comm = m_dc->GetComm();
#endif
  }

  virtual ~Shaper() = default;

public:
  //@{
  //!  @name Functions to get and set shaping parameters

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

  //@}

  sidre::MFEMSidreDataCollection* getDC() { return m_dc; }
  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  virtual bool isValidFormat(const std::string& format) const
  {
    return (format == "stl" || format == "c2c");
  }

public:
  //@{
  //!  @name Functions related to the stages for a given shape

  /// Loads the shape from file into m_surfaceMesh
  virtual void loadShape(const klee::Shape& shape)
  {
    using axom::utilities::string::endsWith;

    SLIC_INFO(fmt::format("{:-^80}",
                          fmt::format(" Loading shape '{}' ", shape.getName())));

    SLIC_ASSERT_MSG(this->isValidFormat(shape.getGeometry().getFormat()),
                    fmt::format("Shape has unsupported format: '{}",
                                shape.getGeometry().getFormat()));

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

  virtual void applyTransforms(const klee::Shape& shape)
  {
    // TODO: Implement this as a set of affine transforms to vertices of mesh
    AXOM_UNUSED_VAR(shape);
  }

  virtual void prepareShapeQuery(klee::Dimensions shapeDimension,
                                 const klee::Shape& shape) = 0;

  virtual void runShapeQuery(const klee::Shape& shape) = 0;

  virtual void applyReplacementRules(const klee::Shape& shape) = 0;

  virtual void finalizeShapeQuery() = 0;

  //@}

public:
  //@{
  //!  @name Functions to generate/adjust volume fractions after all shapes have been applied

  virtual void adjustVolumeFractions() = 0;

  //@}

protected:
  /*!
   * \brief Helper function to get the rank associated with the current process
   *
   * \note This function can be called even in non-mpi configurations
   */
  int getRank() const
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    if(auto* pmesh = static_cast<mfem::ParMesh*>(m_dc->GetMesh()))
    {
      return pmesh->GetMyRank();
    }
    return 0;
#else
    return 0;
#endif
  }

protected:
  const klee::ShapeSet& m_shapeSet;
  sidre::MFEMSidreDataCollection* m_dc;

  mint::Mesh* m_surfaceMesh {nullptr};

  int m_samplesPerKnotSpan {25};
  double m_vertexWeldThreshold {1e-9};

  MPI_Comm m_comm {MPI_COMM_SELF};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPER__HPP_
