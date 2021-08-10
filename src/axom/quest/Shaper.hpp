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
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/mint.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"

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
  Shaper(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc);

  virtual ~Shaper() = default;

public:
  //@{
  //!  @name Functions to get and set shaping parameters

  void setSamplesPerKnotSpan(int nSamples);
  void setVertexWeldThreshold(double threshold);
  void setVerbosity(bool isVerbose) { m_verboseOutput = isVerbose; }

  //@}

  bool isVerbose() const { return m_verboseOutput; }

  sidre::MFEMSidreDataCollection* getDC() { return m_dc; }
  mint::Mesh* getSurfaceMesh() const { return m_surfaceMesh; }

  /*!
   * \brief Predicate to determine if the specified format is valid
   *
   * \param format A string listing the format to check
   */
  virtual bool isValidFormat(const std::string& format) const;

public:
  //@{
  //!  @name Functions related to the stages for a given shape

  /// Loads the shape from file into m_surfaceMesh
  virtual void loadShape(const klee::Shape& shape);

  virtual void applyTransforms(const klee::Shape& shape);

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
  int getRank() const;

protected:
  const klee::ShapeSet& m_shapeSet;
  sidre::MFEMSidreDataCollection* m_dc;

  mint::Mesh* m_surfaceMesh {nullptr};

  int m_samplesPerKnotSpan {25};
  double m_vertexWeldThreshold {1e-9};
  bool m_verboseOutput {false};

  MPI_Comm m_comm {MPI_COMM_SELF};
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPER__HPP_
