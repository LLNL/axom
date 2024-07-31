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
#include "axom/mint/mesh/UnstructuredMesh.hpp"

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
  //!@brief Type for geometry represented by a simplex mesh.
  using SimplexMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;

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
   * Create a new Geometry object based on a bluieprint tetrahedra mesh.
   *
   * \param startProperties the transformable properties before any
   * operators are applied
   * \param simplexMesh a simplex geometry in blueprint format.
   *   The elements should be segments, triangles or tetrahedra.
   * \param topology The \c simplexMesh topology to use.
   * \param operator_ a possibly null operator to apply to the geometry.
   *
   * \internal TODO: Is this the simplex requirement overly restrictive?
   */
  Geometry(const TransformableGeometryProperties &startProperties,
           axom::sidre::Group *simplexMesh,
           const std::string& topology,
           std::shared_ptr<GeometryOperator const> operator_);

  /**
   * Get the format in which the geometry is specified.
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

  axom::sidre::Group* getBlueprintMesh() const { return m_simplexMesh; }

  const std::string& getBlueprintTopology() const { return m_topology; }

  /// Predicate that returns true when the shape has an associated geometry
  /* This seems a poorly named function and probably should be renamed
     to isFromFile.  Especially true in light of my adding memor-based
     geometry.  Maybe this method was orginally from the Shape class.
     But note that Kenny intentionally renamed it from hasPath(),
     to "better match how it will be used.".  Ask Kenny.
  */
  bool hasGeometry() const { return !m_path.empty(); }

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

private:
  TransformableGeometryProperties m_startProperties;

  //!@brief Geometry file format, if it's file-based.
  std::string m_format;

  //!@brief Geometry file path, if it's file-based.
  std::string m_path;

  //!@brief Geometry blueprint simplex mesh, if it's in memory.
  axom::sidre::Group* m_simplexMesh;
  // std::shared_ptr<SimplexMesh> m_simplexMesh;

  //!@brief Topology of the blueprint simplex mesh, if it's in memory.
  std::string m_topology;

  std::shared_ptr<const GeometryOperator> m_operator;
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEE_GEOMETRY_HPP
