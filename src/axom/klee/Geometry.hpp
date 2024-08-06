// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRY_HPP
#define AXOM_KLEE_GEOMETRY_HPP

#include <memory>
#include <string>

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"

namespace axom
{
namespace klee
{
class GeometryOperator;

/**
 * Properties of a geometric object which can be transformed by operators
 */
struct TransformableGeometryProperties
{
  Dimensions dimensions;
  LengthUnit units;
};

/**
 * Compare transformable properties for equality.
 * \param lhs the left-hand-side operand
 * \param rhs the right-hand-side operand
 * \return true if and only if all properties are equal
 */
bool operator==(const TransformableGeometryProperties &lhs,
                const TransformableGeometryProperties &rhs);

/**
 * Compare transformable properties for inequality.
 * \param lhs the left-hand-side operand
 * \param rhs the right-hand-side operand
 * \return false if and only if all properties are equal
 */
inline bool operator!=(const TransformableGeometryProperties &lhs,
                       const TransformableGeometryProperties &rhs)
{
  return !(lhs == rhs);
}

/**
 * Represents the geometry specified in a Shape.
 */
class Geometry
{
public:
  /**
   * Create a new Geometry object based on a file representation.
   *
   * \param startProperties the transformable properties before any
   * operators are applied
   * \param format the format of the file
   * \param path the path of the file
   * \param operator_ a possibly null operator to apply to the geometry.
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           std::string format,
           std::string path,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a new Geometry object based on a blueprint tetrahedral mesh.
   *
   * \param startProperties the transformable properties before any
   * operators are applied
   * \param meshGroup a simplex geometry in blueprint format.
   *   The elements should be segments, triangles or tetrahedra.
   * \param topology The \c meshGroup topology to use.
   * \param operator_ a possibly null operator to apply to the geometry.
   *
   * \internal TODO: Is this the simplex requirement overly restrictive?
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           axom::sidre::Group *meshGroup,
           const std::string& topology,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Create a new sphere Geometry object.
   *
   * \param startProperties the transformable properties before any
   * operators are applied
   * \param sphere Analytical sphere specifications
   * \param generationCount Number of generations to use for
   *        discretizing the sphere.
   * \param operator_ a possibly null operator to apply to the geometry.
   *
   * \internal TODO: Is this the simplex requirement overly restrictive?
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           const axom::primal::Sphere<double, 3>& sphere,
           axom::IndexType generationCount,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Get the format in which the geometry is specified.
   *
   * Values are:
   * - "c2c" = C2C file
   * - "proe" = ProE file
   * - "memory-blueprint" = Blueprint tetrahedral mesh in memory
   * - "memory-xy" = discretized 2D, non-negative function as an
   *   array of points in memory
   * - "sphere3D" = 3D sphere, as \c primal::Sphere<double,3>
   * - "cone3D" = 3D cone, as \c primal::Cone<double,3>
   *   "cylinder3D" = 3D cylinder, as \c primal::Cylinder<double,3>
   *
   * \return the format of the shape
   */
  const std::string &getFormat() const { return m_format; }

  /**
   * Get the path at which to find the specification of the geometry,
   * for geometries stored in files.
   *
   * \return the path to the geometry file
   */
  const std::string &getPath() const { return m_path; }

  /**
   * @brief Return the blueprint mesh, for formats that are specified
   * by a blueprint mesh or have been converted to a blueprint mesh.
   */
  axom::sidre::Group* getBlueprintMesh() const;

  /**
   * @brief Return the blueprint mesh topology, for formats that are specified
   * by a blueprint mesh or have been converted to a blueprint mesh.
   */
  const std::string& getBlueprintTopology() const;

  /*! @brief Predicate that returns true when the shape has an associated geometry

    A false means that this is set up to determine volume fractions without
    computing on the geometry.

    TODO: We should just create a new format to represent getting
    volume fractions without geometries.  Or move this logic into
    Shape, because it's confusing to have a Geometry that has no
    geometry.
  */
  bool hasGeometry() const;

  /**
   * Get a GeometryOperator to apply to this geometry. Can be null.
   *
   * \return a potentially null operator to apply to the geometry
   */
  std::shared_ptr<GeometryOperator const> const &getGeometryOperator() const
  {
    return m_operator;
  }

  /**
   * Get the initial transformable properties of this geometry
   *
   * \return the initial transformable properties of this geometry
   */
  const TransformableGeometryProperties &getStartProperties() const
  {
    return m_startProperties;
  }

  /**
   * Get the final transformable properties of this geometry after
   * operators are applied
   *
   * \return the initial transformable properties of this geometry
   */
  TransformableGeometryProperties getEndProperties() const;

  /**
   @brief Return the number of generations for discretization of
   analytical curves.

   This number is unused for geometries that are specified in discrete
   form.
  */
  axom::IndexType getGenerationCount() const
  {
    return m_generationCount;
  }

  /**
   @brief Return the sphere geometry, when the Geometry
   represents an alalytical sphere.
  */
  const axom::primal::Sphere<double, 3>& getSphere() const
  {
    return m_sphere;
  }

private:
  TransformableGeometryProperties m_startProperties;

  //!@brief Geometry file format.
  std::string m_format;

  //!@brief Geometry file path, if it's file-based.
  std::string m_path;

  //!@brief Geometry blueprint simplex mesh, when/if it's in memory.
  axom::sidre::Group* m_meshGroup;

  //!@brief Topology of the blueprint simplex mesh, if it's in memory.
  std::string m_topology;

  //!@brief The analytical sphere, if used.
  axom::primal::Sphere<double, 3> m_sphere;

#if 0
  //!@brief The analytical cylinder, if used.
  axom::primal::Sphere<double, 3> m_cylinder;

  //!@brief The analytical cone, if used.
  axom::primal::Sphere<double, 3> m_cone;

  //!@brief The discrete 2D curve, if used.
  axom::Array<primal::Point<double, 2> m_discreteFunction;
#endif

  //!@brief Generations of refinement for discretizing analytical shapes
  // and surfaces of revolutions.
  axom::IndexType m_generationCount = 0;

  std::shared_ptr<const GeometryOperator> m_operator;
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEE_GEOMETRY_HPP
