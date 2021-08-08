// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Shaper.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "fmt/fmt.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "mfem.hpp"

namespace axom
{
namespace quest
{
Shaper::Shaper(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
  : m_shapeSet(shapeSet)
  , m_dc(dc)
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = m_dc->GetComm();
#endif
}

void Shaper::setSamplesPerKnotSpan(int nSamples)
{
  using axom::utilities::clampLower;
  SLIC_WARNING_IF(
    nSamples < 1,
    fmt::format(
      "Samples per knot span must be at least 1. Provided value was {}",
      nSamples));

  m_samplesPerKnotSpan = clampLower(nSamples, 1);
}

void Shaper::setVertexWeldThreshold(double threshold)
{
  SLIC_WARNING_IF(
    threshold <= 0.,
    fmt::format(
      "Vertex weld threshold should be positive Provided value was {}",
      threshold));

  m_vertexWeldThreshold = threshold;
}

bool Shaper::isValidFormat(const std::string& format) const
{
  return (format == "stl" || format == "c2c");
}

void Shaper::loadShape(const klee::Shape& shape)
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

void Shaper::applyTransforms(const klee::Shape& shape)
{
  // TODO: Implement this as a set of affine transforms to vertices of mesh
  AXOM_UNUSED_VAR(shape);
}

// ----------------------------------------------------------------------------

int Shaper::getRank() const
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

}  // end namespace quest
}  // end namespace axom
