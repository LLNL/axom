// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Shaper.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal/operators/split.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/Shaper.hpp"
#include "axom/quest/DiscreteShape.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace quest
{

// These were needed for linking - but why? They are constexpr.
constexpr int Shaper::DEFAULT_SAMPLES_PER_KNOT_SPAN;
constexpr double Shaper::MINIMUM_PERCENT_ERROR;
constexpr double Shaper::MAXIMUM_PERCENT_ERROR;
constexpr double Shaper::DEFAULT_VERTEX_WELD_THRESHOLD;

#if defined(AXOM_USE_MFEM)
Shaper::Shaper(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
  : m_shapeSet(shapeSet)
  , m_dc(dc)
{
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = m_dc->GetComm();
  #endif
  m_cellCount = m_dc->GetMesh()->GetNE();
}
#endif

Shaper::Shaper(const klee::ShapeSet& shapeSet,
               sidre::Group* bpGrp,
               const std::string& topo)
  : m_shapeSet(shapeSet)
  , m_bpGrp(bpGrp)
  , m_bpTopo(topo.empty() ? bpGrp->getGroup("topologies")->getGroupName(0) : topo)
  , m_bpNode(nullptr)
  , m_comm(MPI_COMM_WORLD)
{
  SLIC_ASSERT(m_bpTopo != sidre::InvalidName);

  auto* topologies = m_bpGrp->getGroup("topologies");
  auto* topoGrp = topologies->getGroup(m_bpTopo);
  auto* coordsetName = topoGrp->getView("coordset");
  std::string coordsName = coordsetName->getString();
#if 0
  std::string coordsName =
    m_bpGrp->getView(axom::fmt::format("topologies/{}/coordset", m_bpTopo))
      ->getString();
#endif
  auto* coordsView =
    m_bpGrp->getView(axom::fmt::format("coordsets/{}/values/x", coordsName));
  m_cellCount = coordsView->getNumElements();
}

void Shaper::setSamplesPerKnotSpan(int nSamples)
{
  using axom::utilities::clampLower;
  SLIC_WARNING_IF(
    nSamples < 1,
    axom::fmt::format(
      "Samples per knot span must be at least 1. Provided value was {}",
      nSamples));

  m_samplesPerKnotSpan = clampLower(nSamples, 1);
}

void Shaper::setVertexWeldThreshold(double threshold)
{
  SLIC_WARNING_IF(
    threshold <= 0.,
    axom::fmt::format(
      "Vertex weld threshold should be positive Provided value was {}",
      threshold));

  m_vertexWeldThreshold = threshold;
}

void Shaper::setPercentError(double percent)
{
  using axom::utilities::clampVal;
  SLIC_WARNING_IF(
    percent <= MINIMUM_PERCENT_ERROR,
    axom::fmt::format("Percent error must be greater than {}. Provided value "
                      "was {}. Dynamic refinement will not be used.",
                      MINIMUM_PERCENT_ERROR,
                      percent));
  SLIC_WARNING_IF(percent > MAXIMUM_PERCENT_ERROR,
                  axom::fmt::format(
                    "Percent error must be less than {}. Provided value was {}",
                    MAXIMUM_PERCENT_ERROR,
                    percent));
  if(percent <= MINIMUM_PERCENT_ERROR)
  {
    m_refinementType = DiscreteShape::RefinementUniformSegments;
  }
  m_percentError =
    clampVal(percent, MINIMUM_PERCENT_ERROR, MAXIMUM_PERCENT_ERROR);
}

void Shaper::setRefinementType(Shaper::RefinementType t)
{
  m_refinementType = t;
}

bool Shaper::isValidFormat(const std::string& format) const
{
  return (format == "stl" || format == "proe" || format == "c2c" ||
          format == "blueprint-tets" || format == "tet3D" ||
          format == "hex3D" || format == "plane3D" || format == "sphere3D" ||
          format == "vor3D" || format == "none");
}

void Shaper::loadShape(const klee::Shape& shape)
{
  AXOM_ANNOTATE_SCOPE("loadShape");

  // Do not save the revolved volume in the default shaper.
  double revolved = 0.;
  loadShapeInternal(shape, m_percentError, revolved);
}

void Shaper::loadShapeInternal(const klee::Shape& shape,
                               double percentError,
                               double& revolvedVolume)
{
  internal::ScopedLogLevelChanger logLevelChanger(
    this->isVerbose() ? slic::message::Debug : slic::message::Warning);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(" Loading shape '{}' ", shape.getName())));

  SLIC_ASSERT_MSG(this->isValidFormat(shape.getGeometry().getFormat()),
                  axom::fmt::format("Shape has unsupported format: '{}",
                                    shape.getGeometry().getFormat()));

  // Code for discretizing shapes has been factored into DiscreteShape class.
  DiscreteShape discreteShape(shape, m_dataStore.getRoot(), m_shapeSet.getPath());
  discreteShape.setVertexWeldThreshold(m_vertexWeldThreshold);
  discreteShape.setRefinementType(m_refinementType);
  if(percentError > 0)
  {
    discreteShape.setPercentError(percentError);
  }
  m_surfaceMesh = discreteShape.createMeshRepresentation();
  revolvedVolume = discreteShape.getRevolvedVolume();
}

// ----------------------------------------------------------------------------

int Shaper::getRank() const
{
#if defined(AXOM_USE_MPI)
  int rank = -1;
  MPI_Comm_rank(m_comm, &rank);
  return rank;
#endif
  return 0;
}

double Shaper::allReduceSum(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_SUM, m_comm);
  return global;
#endif
  return val;
}

}  // end namespace quest
}  // end namespace axom
