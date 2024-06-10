// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Shaper.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include "axom/fmt.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "mfem.hpp"

namespace axom
{
namespace quest
{
namespace internal
{
/*!
 * \brief Implementation of a GeometryOperatorVisitor for processing klee shape operators
 *
 * This class extracts the matrix form of supported operators and marks the operator as unvalid otherwise
 * To use, check the \a isValid() function after visiting and then call the \a getMatrix() function.
 */
class AffineMatrixVisitor : public klee::GeometryOperatorVisitor
{
public:
  AffineMatrixVisitor() : m_matrix(4, 4) { }

  void visit(const klee::Translation& translation) override
  {
    m_matrix = translation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Rotation& rotation) override
  {
    m_matrix = rotation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Scale& scale) override
  {
    m_matrix = scale.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::UnitConverter& converter) override
  {
    m_matrix = converter.toMatrix();
    m_isValid = true;
  }

  void visit(const klee::CompositeOperator&) override
  {
    SLIC_WARNING("CompositeOperator not supported for Shaper query");
    m_isValid = false;
  }
  void visit(const klee::SliceOperator&) override
  {
    SLIC_WARNING("SliceOperator not yet supported for Shaper query");
    m_isValid = false;
  }

  const numerics::Matrix<double>& getMatrix() const { return m_matrix; }
  bool isValid() const { return m_isValid; }

private:
  bool m_isValid {false};
  numerics::Matrix<double> m_matrix;
};

}  // end namespace internal

// These were needed for linking - but why? They are constexpr.
constexpr int Shaper::DEFAULT_SAMPLES_PER_KNOT_SPAN;
constexpr double Shaper::MINIMUM_PERCENT_ERROR;
constexpr double Shaper::MAXIMUM_PERCENT_ERROR;
constexpr double Shaper::DEFAULT_VERTEX_WELD_THRESHOLD;

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
    m_refinementType = RefinementUniformSegments;
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
          format == "none");
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
  using axom::utilities::string::endsWith;

  internal::ScopedLogLevelChanger logLevelChanger(
    this->isVerbose() ? slic::message::Debug : slic::message::Warning);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(" Loading shape '{}' ", shape.getName())));

  const std::string& file_format = shape.getGeometry().getFormat();
  SLIC_ASSERT_MSG(
    this->isValidFormat(file_format),
    axom::fmt::format("Shape has unsupported format: '{}", file_format));

  if(!shape.getGeometry().hasGeometry())
  {
    SLIC_DEBUG(
      axom::fmt::format("Current shape '{}' of material '{}' has no geometry",
                        shape.getName(),
                        shape.getMaterial()));
    return;
  }

  std::string shapePath = m_shapeSet.resolvePath(shape.getGeometry().getPath());
  SLIC_INFO("Reading file: " << shapePath << "...");

  // Initialize revolved volume.
  revolvedVolume = 0.;

  if(endsWith(shapePath, ".stl"))
  {
    SLIC_ASSERT_MSG(
      file_format == "stl",
      axom::fmt::format(" '{}' format requires .stl file type", file_format));

    quest::internal::read_stl_mesh(shapePath, m_surfaceMesh, m_comm);
    // Transform the coordinates of the linearized mesh.
    applyTransforms(shape);
  }
  else if(endsWith(shapePath, ".proe"))
  {
    SLIC_ASSERT_MSG(
      file_format == "proe",
      axom::fmt::format(" '{}' format requires .proe file type", file_format));

    quest::internal::read_pro_e_mesh(shapePath, m_surfaceMesh, m_comm);
  }
#ifdef AXOM_USE_C2C
  else if(endsWith(shapePath, ".contour"))
  {
    SLIC_ASSERT_MSG(
      file_format == "c2c",
      axom::fmt::format(" '{}' format requires .contour file type", file_format));

    // Get the transforms that are being applied to the mesh. Get them
    // as a single concatenated matrix.
    auto transform = getTransforms(shape);

    // Pass in the transform so any transformations can figure into
    // computing the revolved volume.
    if(m_refinementType == RefinementDynamic &&
       percentError > MINIMUM_PERCENT_ERROR)
    {
      quest::internal::read_c2c_mesh_non_uniform(shapePath,
                                                 transform,
                                                 percentError,
                                                 m_vertexWeldThreshold,
                                                 m_surfaceMesh,
                                                 revolvedVolume,  // output arg
                                                 m_comm);
    }
    else
    {
      quest::internal::read_c2c_mesh_uniform(shapePath,
                                             transform,
                                             m_samplesPerKnotSpan,
                                             m_vertexWeldThreshold,
                                             m_surfaceMesh,
                                             revolvedVolume,  // output arg
                                             m_comm);
    }

    // Transform the coordinates of the linearized mesh.
    applyTransforms(transform);
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("Unsupported filetype for this Axom configuration. "
                        "Provided file was '{}', with format '{}'",
                        shapePath,
                        file_format));
  }
}

numerics::Matrix<double> Shaper::getTransforms(const klee::Shape& shape) const
{
  const auto identity4x4 = numerics::Matrix<double>::identity(4);
  numerics::Matrix<double> transformation(identity4x4);
  auto& geometryOperator = shape.getGeometry().getGeometryOperator();
  auto composite =
    std::dynamic_pointer_cast<const klee::CompositeOperator>(geometryOperator);
  if(composite)
  {
    // Concatenate the transformations
    for(auto op : composite->getOperators())
    {
      // Use visitor pattern to extract the affine matrix from supported operators
      internal::AffineMatrixVisitor visitor;
      op->accept(visitor);
      if(!visitor.isValid())
      {
        continue;
      }
      const auto& matrix = visitor.getMatrix();
      numerics::Matrix<double> res(identity4x4);
      numerics::matrix_multiply(matrix, transformation, res);
      transformation = res;
    }
  }
  return transformation;
}

void Shaper::applyTransforms(const klee::Shape& shape)
{
  // Concatenate the transformations
  numerics::Matrix<double> transformation = getTransforms(shape);
  // Transform the surface mesh coordinates.
  applyTransforms(transformation);
}

void Shaper::applyTransforms(const numerics::Matrix<double>& transformation)
{
  AXOM_ANNOTATE_SCOPE("applyTransforms");

  // Apply transformation to coordinates of each vertex in mesh
  if(!transformation.isIdentity())
  {
    const int spaceDim = m_surfaceMesh->getDimension();
    const int numSurfaceVertices = m_surfaceMesh->getNumberOfNodes();
    double* x = m_surfaceMesh->getCoordinateArray(mint::X_COORDINATE);
    double* y = m_surfaceMesh->getCoordinateArray(mint::Y_COORDINATE);
    double* z = spaceDim > 2
      ? m_surfaceMesh->getCoordinateArray(mint::Z_COORDINATE)
      : nullptr;

    double xformed[4];
    for(int i = 0; i < numSurfaceVertices; ++i)
    {
      double coords[4] = {x[i], y[i], (z == nullptr ? 0. : z[i]), 1.};
      numerics::matrix_vector_multiply(transformation, coords, xformed);
      x[i] = xformed[0];
      y[i] = xformed[1];
      if(z != nullptr)
      {
        z[i] = xformed[2];
      }
    }
  }
}

// ----------------------------------------------------------------------------

int Shaper::getRank() const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(auto* pmesh = static_cast<mfem::ParMesh*>(m_dc->GetMesh()))
  {
    return pmesh->GetMyRank();
  }
#endif
  return 0;
}

double Shaper::allReduceSum(double val) const
{
  if(m_dc != nullptr && m_dc->GetNumProcs() > 1)
  {
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    double global;
    MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_SUM, m_comm);
    return global;
#endif
  }
  return val;
}

}  // end namespace quest
}  // end namespace axom
