// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

bool Shaper::isValidFormat(const std::string& format) const
{
  return (format == "stl" || format == "c2c");
}

void Shaper::loadShape(const klee::Shape& shape)
{
  using axom::utilities::string::endsWith;

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(" Loading shape '{}' ", shape.getName())));

  SLIC_ASSERT_MSG(this->isValidFormat(shape.getGeometry().getFormat()),
                  axom::fmt::format("Shape has unsupported format: '{}",
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
      axom::fmt::format("Unsupported filetype for this Axom configuration. "
                        "Provided file was '{}'",
                        shapePath));
  }
}

void Shaper::applyTransforms(const klee::Shape& shape)
{
  auto& geometryOperator = shape.getGeometry().getGeometryOperator();
  auto composite =
    std::dynamic_pointer_cast<const klee::CompositeOperator>(geometryOperator);
  if(composite)
  {
    // Get surface mesh coordinates
    const int spaceDim = m_surfaceMesh->getDimension();
    const int numSurfaceVertices = m_surfaceMesh->getNumberOfNodes();
    double* x = m_surfaceMesh->getCoordinateArray(mint::X_COORDINATE);
    double* y = m_surfaceMesh->getCoordinateArray(mint::Y_COORDINATE);
    double* z = spaceDim > 2
      ? m_surfaceMesh->getCoordinateArray(mint::Z_COORDINATE)
      : nullptr;

    // Loop through operators and apply transformations to vertices
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

      // Apply transformation to coordinates of each vertex in mesh
      double xformed[4];
      for(int i = 0; i < numSurfaceVertices; ++i)
      {
        double coords[4] = {x[i], y[i], (z == nullptr ? 0. : z[i]), 1.};
        numerics::matrix_vector_multiply(matrix, coords, xformed);
        x[i] = xformed[0];
        y[i] = xformed[1];
        if(z != nullptr)
        {
          z[i] = xformed[2];
        }
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
  return 0;
#else
  return 0;
#endif
}

double Shaper::allReduceSum(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_SUM, m_comm);
  return global;
#else
  return val;
#endif
}

}  // end namespace quest
}  // end namespace axom
