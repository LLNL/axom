// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_DISCRETE_SHAPE_HPP
#define AXOM_KLEE_DISCRETE_SHAPE_HPP

#include <string>
#include <vector>

#include "axom/klee/Shape.hpp"
#include "axom/mint/mesh/Mesh.hpp"

#if defined(AXOM_USE_MPI)
  #include "mpi.h"
#endif

namespace axom
{
namespace quest
{
/*!
  @brief Post-processed klee::Shape, with the geometry discretized and transformed
  according to the Shape's operators.

  TODO: Move this class into internal namespace.
*/
class DiscreteShape
{
public:
  /// Refinement type.
  using RefinementType = enum { RefinementUniformSegments, RefinementDynamic };

  using Point3D = axom::primal::Point<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;
  using TetType = axom::primal::Tetrahedron<double, 3>;
  using OctType = axom::primal::Octahedron<double, 3>;
  using HexType = axom::primal::Hexahedron<double, 3>;

  static constexpr int DEFAULT_SAMPLES_PER_KNOT_SPAN {25};
  static constexpr double MINIMUM_PERCENT_ERROR {0.};
  static constexpr double MAXIMUM_PERCENT_ERROR {100.};
  static constexpr double DEFAULT_VERTEX_WELD_THRESHOLD {1e-9};

  /*!
    @brief Constructor.

    @param shape The Klee specifications for the shape.
    @param prefixPath Path prefix for shape files specified
      with a relative path.
    @param parentGroup Group under which to put the discrete mesh.
      If null, don't use sidre.

    Refinement type is set to DiscreteShape::RefinementUniformSegments
    and percent erro is set to 0.  See setPercentError() and
    setRefinementType().
  */
  DiscreteShape(const axom::klee::Shape& shape,
                const std::string& prefixPath = {},
                axom::sidre::Group* parentGroup = nullptr);

  virtual ~DiscreteShape() { clearInternalData(); }

  //@{
  //! @name Functions to get and set shaping parameters

  void setPrefixPath(const std::string& prefixPath)
  {
    m_prefixPath = prefixPath;
  }

  /*!
    @brief Set the refinement type.
    Refinement type is used for shaping with C2C contours.
  */
  void setRefinementType(RefinementType refinementType)
  {
    m_refinementType = refinementType;
  }

  void setSamplesPerKnotSpan(int nSamples);
  void setVertexWeldThreshold(double threshold);
  /*!
    @brief Set the percentage error tolerance.

    If percent <= MINIMUM_PERCENT_ERROR, the refinement type
    will change to DiscreteShape::RefinementUniformSegments.
  */
  void setPercentError(double percent);

  //@}

#if defined(AXOM_USE_MPI)
  /**
   @brief Set the MPI communicator used when reading C2C files.
   */
  void setMPICommunicator(MPI_Comm comm) { m_comm = comm; }
#endif

  /*!
    Get the name of this shape.
    \return the shape's name
  */
  const axom::klee::Shape& getShape() const { return m_shape; }

  /*!
    \brief Get the discrete mesh representation.

    If the discrete mesh isn't generated yet (for analytical shapes)
    generate it.
  */
  std::shared_ptr<mint::Mesh> createMeshRepresentation();

  //!@brief Get the discrete mesh representation.
  std::shared_ptr<mint::Mesh> getMeshRepresentation() const
  {
    return m_meshRep;
  }

  /*!
     \brief Get the revolved volume for volumes of revolution.
  */
  double getRevolvedVolume() const { return m_revolvedVolume; }

private:
  const axom::klee::Shape& m_shape;

  /*!
    \brief Discrete mesh representation.

    This is either an internally generated mesh representing a
    discretized analytical shape or a modifiable copy of the
    discrete input geometry for geometries specified as a
    discrete mesh.
  */
  std::shared_ptr<axom::mint::Mesh> m_meshRep;

  //!@brief Internal DataStore for working space
  axom::sidre::DataStore m_dataStore;

  //!@brief Sidre store for m_meshRep.
  axom::sidre::Group* m_sidreGroup {nullptr};

  //!@brief Prefix for disc files with relative path.
  std::string m_prefixPath;

  //@{
  //!@name Various parameters for discretization of analytical shapes.
  RefinementType m_refinementType;
  double m_percentError {MINIMUM_PERCENT_ERROR};
  int m_samplesPerKnotSpan {DEFAULT_SAMPLES_PER_KNOT_SPAN};
  double m_vertexWeldThreshold {DEFAULT_VERTEX_WELD_THRESHOLD};
  double m_revolvedVolume {0.0};
  //@}

#if defined(AXOM_USE_MPI)
  MPI_Comm m_comm {MPI_COMM_WORLD};
#endif

  void applyTransforms();
  numerics::Matrix<double> getTransforms() const;

  //! @brief Returns the full geometry path.
  std::string resolvePath() const;

  /*!
    @brief Set the parent group for this object to store data.
  */
  void setParentGroup(axom::sidre::Group* parentGroup);

  //!@brief Return a 3x3 matrix that rotate coordinates from the x-axis to the given direction.
  numerics::Matrix<double> vorAxisRotMatrix(const Vector3D& dir);

  void clearInternalData();
};

}  // namespace quest
}  // namespace axom

#endif
