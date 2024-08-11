// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"

#include "axom/klee/GeometryOperators.hpp"

#include <utility>

namespace axom
{
namespace klee
{
bool operator==(const TransformableGeometryProperties& lhs,
                const TransformableGeometryProperties& rhs)
{
  return lhs.dimensions == rhs.dimensions && lhs.units == rhs.units;
}

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   std::string format,
                   std::string path,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format(std::move(format))
  , m_path(std::move(path))
  , m_levelOfRefinement(0)
  , m_operator(std::move(operator_))
{ }

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   axom::sidre::Group* meshGroup,
                   const std::string& topology,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("memory-blueprint")
  , m_path()
  , m_meshGroup(meshGroup)
  , m_topology(topology)
  , m_levelOfRefinement(0)
  , m_operator(std::move(operator_))
{ }

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::primal::Hexahedron<double, 3>& hex,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("hex3D")
  , m_path()
  , m_meshGroup(nullptr)
  , m_topology()
  , m_hex(hex)
  , m_levelOfRefinement(0)
  , m_operator(std::move(operator_))
{ }

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const Sphere3D& sphere,
                   axom::IndexType levelOfRefinement,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("sphere3D")
  , m_path()
  , m_meshGroup(nullptr)
  , m_topology()
  , m_levelOfRefinement(levelOfRefinement)
  , m_sphere(sphere)
  , m_operator(std::move(operator_))
{ }

Geometry::Geometry(const TransformableGeometryProperties& startProperties,
                   const axom::Array<double, 2>& discreteFunction,
                   const Point3D& vorBase,
                   const Vector3D& vorDirection,
                   axom::IndexType levelOfRefinement,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("vor3D")
  , m_path()
  , m_meshGroup(nullptr)
  , m_topology()
  , m_levelOfRefinement(levelOfRefinement)
  , m_sphere()
  , m_discreteFunction(discreteFunction)
  , m_vorBase(vorBase)
  , m_vorDirection(vorDirection)
  , m_operator(std::move(operator_))
{ }

bool Geometry::hasGeometry() const
{
  bool isInMemory = m_format == "memory-blueprint" || m_format == "sphere3D" ||
    m_format == "cone3D" || m_format == "cylinder3D";
  if(isInMemory)
  {
    return true;
  }
  return !m_path.empty();
}

TransformableGeometryProperties Geometry::getEndProperties() const
{
  if(m_operator)
  {
    return m_operator->getEndProperties();
  }
  return m_startProperties;
}

axom::sidre::Group* Geometry::getBlueprintMesh() const
{
  SLIC_ASSERT_MSG(
    m_meshGroup,
    axom::fmt::format(
      "The Geometry format '{}' is not specified "
      "as a blueprint mesh and/or has not been converted into one.",
      m_format));
  return m_meshGroup;
}

const std::string& Geometry::getBlueprintTopology() const
{
  SLIC_ASSERT_MSG(
    m_meshGroup,
    axom::fmt::format(
      "The Geometry format '{}' is not specified "
      "as a blueprint mesh and/or has not been converted into one.",
      m_format));
  return m_topology;
}

}  // namespace klee
}  // namespace axom
